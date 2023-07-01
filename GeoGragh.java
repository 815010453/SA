package lyd.SA;

import org.checkerframework.checker.units.qual.A;
import org.geotools.data.FeatureWriter;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.referencing.CRS;
import org.geotools.referencing.GeodeticCalculator;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.*;
import org.opengis.feature.Property;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static org.geotools.data.Transaction.AUTO_COMMIT;

public class GeoGragh {
    private final HashMap<Integer, GeoVertex> id_vertex;  // 点号对应的节点字典 {id1: vertex1, id2: vertex2, ...}
    private final HashMap<Integer, GeoEdge> id_edge; // 边号对应的边字典 {id1: edge1, id2: edge2}
    private final HashMap<ArrayList<GeoVertex>, GeoEdge> vertices_edge;  // 点号对应的边号字典 {(vertex1, vertex2): edge1, (vertex3, vertex4): edge2,...}
    private final HashMap<Integer, GeoSegment> id_segment;
    private final String name; // 图名
    private final HashMap<ArrayList<Double>, GeoVertex> vertexCoord; // 点的坐标对应的点号
    private String proj; //坐标系字符串
    private final HashMap<GeoVertex, ArrayList<GeoVertex>> subVertex; //辅助字典1
    private final HashMap<GeoVertex, GeoVertex> auxVertex; //辅助字典2
    private final ArrayList<ArrayList<GeoVertex>> roadDict; // 当前图中路的路径点几何
    private final HashSet<GeoSegment> circleRoad; //环路
    private ArrayList<GeoVertex> node; //顶点

    public GeoGragh(String name) {
        this.id_vertex = new HashMap<>();
        this.id_edge = new HashMap<>();
        this.vertices_edge = new HashMap<>();
        this.name = name;
        this.vertexCoord = new HashMap<>();
        this.proj = "";
        this.roadDict = new ArrayList<>();
        this.subVertex = new HashMap<>();
        this.auxVertex = new HashMap<>();
        this.circleRoad = new HashSet<>();
        this.node = new ArrayList<>();
        this.id_segment = new HashMap<>();
    }

    //通过边文件路径名以及图名构建该图
    public GeoGragh(String name, String path) throws Exception {
        //点与边的id号均从1开始
        this(name);
        File file = new File(path);
        // 坐标系WKT
        FileReader fileReader = new FileReader(path.replace(".shp", ".prj"));
        this.proj = new Scanner(fileReader).nextLine();
        System.out.println(this.proj);
        fileReader.close();
        //创建要素类
        FileDataStore store = FileDataStoreFinder.getDataStore(file);
        SimpleFeatureSource featureSource = store.getFeatureSource();
        SimpleFeatureIterator items = featureSource.getFeatures().features();
        int vertexId = 1;
        int edgeId = 1;
        //遍历要素迭代器
        while (items.hasNext()) {
            //获取单个要素
            SimpleFeature sf = items.next();
            Collection<Property> properties = sf.getProperties();
            ArrayList<GeoEdge> geoEdges = new ArrayList<>();
            //属性
            HashMap<String, String> att = new HashMap<>();
            //循环遍历多边
            for (Property property : properties) {
                String shpFieldName = String.valueOf(property.getName());
                String shpFieldVal = String.valueOf(property.getValue());
                if ("the_geom".equals(shpFieldName)) {
                    //坐标数组，需根据要素类型分类讨论
                    ArrayList<Coordinate[]> coordinates = new ArrayList<>();
                    //几何要素类型
                    if (property.getValue() instanceof MultiLineString) {
                        //多线
                        int numGeometries = ((MultiLineString) (property.getValue())).getNumGeometries();
                        for (int i = 0; i < numGeometries; i++) {
                            Coordinate[] coord = ((MultiLineString) (property.getValue())).getGeometryN(i).getCoordinates();
                            coordinates.add(coord);
                        }
                    } else if (property.getValue() instanceof MultiPolygon) {
                        int numGeometries = ((MultiPolygon) (property.getValue())).getNumGeometries();
                        for (int i = 0; i < numGeometries; i++) {
                            Coordinate[] coord = ((MultiPolygon) (property.getValue())).getGeometryN(i).getCoordinates();
                            coordinates.add(coord);
                        }
                    } else if (property.getValue() instanceof LineString) {
                        Coordinate[] coord = ((LineString) (property.getValue())).getGeometryN(0).getCoordinates();
                        coordinates.add(coord);

                    } else if (property.getValue() instanceof Polygon) {
                        int numGeometries = ((Polygon) (property.getValue())).getNumGeometries();
                        for (int i = 0; i < numGeometries; i++) {
                            Coordinate[] coord = ((Polygon) (property.getValue())).getGeometryN(i).getCoordinates();
                            coordinates.add(coord);
                        }
                    }
                    for (Coordinate[] cs : coordinates) {
                        //前置节点与当前节点
                        GeoVertex preVertex = null;
                        GeoVertex nowVertex;
                        for (Coordinate c : cs) {
                            ArrayList<Double> temp = new ArrayList<>();
                            temp.add(c.getX());
                            temp.add(c.getY());
                            double[] tempCoord = {c.getX(), c.getY()};
                            //当前点是否在图中
                            if (preVertex == null) {
                                // 如果前置节点不存在，则查找前置节点
                                preVertex = this.vertexCoord.containsKey(temp) ? this.vertexCoord.get(temp) : new GeoVertex(vertexId++, new HashMap<>(), tempCoord);
                                continue;
                            } else {
                                nowVertex = this.vertexCoord.containsKey(temp) ? this.vertexCoord.get(temp) : new GeoVertex(vertexId++, new HashMap<>(), tempCoord);
                            }
                            // 通过前置节点和现在的节点构建该边，要确保该边不在图中，如果在图中则跳过。
                            if (findEdgeVertices(preVertex, nowVertex) == null) {
                                GeoEdge nowEdge = new GeoEdge(edgeId++, preVertex, nowVertex, new HashMap<>());
                                geoEdges.add(nowEdge);
                                this.addEdge(nowEdge);
                            }
                            preVertex = nowVertex;
                        }
                    }
                } else {
                    //获取属性
                    att.put(shpFieldName, shpFieldVal);
                }
            }
            for (GeoEdge e : geoEdges) {
                e.setEdgeAttribute(att);
            }
        }
        items.close();
        store.dispose();
        //构建路段(两个node之间的边的集合) node的度大于2
        for (GeoVertex v : id_vertex.values()) {
            if (v.getConVertex().size() > 2) node.add(v);
        }
        for(GeoVertex v : id_vertex.values()){
            if (v.getConVertex().size()==0) System.out.println("WTFFFF");
        }
        //图中需要删除的环路点
        ArrayList<GeoVertex> removeVertices = new ArrayList<>();
        System.out.println(node.size());
        int segId = 1;
        for (GeoVertex v : node) {
            ArrayList<GeoVertex> conVertices = v.getConVertex();
            //有可能这个segment重新回到了当前node（环路）
            HashSet<GeoVertex> replication = new HashSet<>();
            //每一个t可以找出一个segment
            for (GeoVertex t : conVertices) {
                if (replication.contains(t)) continue;
                ArrayList<GeoVertex> vertices = new ArrayList<>();
                GeoVertex pre;
                GeoVertex next = t;
                ArrayList<GeoEdge> edges = new ArrayList<>();
                GeoEdge e = findEdgeVertices(v, t);
                edges.add(e);
                vertices.add(v);
                pre = v;
                while (next.getConVertex().size() == 2) {
                    int index;
                    if (next.getConVertex().get(0) == pre) index = 1;
                    else index = 0;
                    pre = next;
                    next = pre.getConVertex().get(index);
                    e = findEdgeVertices(pre, next);
                    edges.add(e);
                    vertices.add(pre);
                }
                vertices.add(next);
                GeoSegment seg = new GeoSegment(segId++, edges, vertices);
                id_segment.put(seg.getId(), seg);
                if (next == v) {
                    //找到一条环路
                    replication.add(pre);
                    this.circleRoad.add(seg);
                    roadDict.add(vertices);
                    //在图中构建环路
                    ArrayList<GeoVertex> rvertices = new ArrayList<>(vertices);
                    rvertices.remove(v);
                    rvertices.remove(v);
                    removeVertices.addAll(rvertices);
                } else v.addSegment(seg);
            }
        }
        for (GeoVertex v : removeVertices) {
            this.removeVertex(v);
        }
        //更新node
        node.clear();
        for (GeoVertex v : id_vertex.values()) {
            if (v.getConVertex().size() > 2) node.add(v);
        }
        System.out.println(node.size());
        System.out.println(roadDict);
    }

    public ArrayList<ArrayList<GeoVertex>> getRoadDict() {
        return roadDict;
    }

    @Override
    public String toString() {
        return "图名：" + name + "坐标系字符串：" + proj;
    }

    public void addVertex(GeoVertex vertex) {
        this.id_vertex.put(vertex.getId(), vertex);
        ArrayList<Double> t = new ArrayList<>();
        t.add(vertex.getCoord()[0]);
        t.add(vertex.getCoord()[1]);
        this.vertexCoord.put(t, vertex);
    }


    public void removeVertex(GeoVertex vertex) {
        int id = vertex.getId();
        id_vertex.remove(id);
        ArrayList<GeoVertex> conVertices = new ArrayList<>(vertex.getConVertex());
        for (GeoVertex v : conVertices) {
            v.remove_conVertex(vertex);
            GeoEdge edge = this.findEdgeVertices(vertex, v);
            ArrayList<Double> t = new ArrayList<>();
            t.add(vertex.getCoord()[0]);
            t.add(vertex.getCoord()[1]);
            vertexCoord.remove(t);
            id_edge.remove(edge.getId());
        }
    }

    public void addEdge(GeoEdge edge) {
        GeoVertex[] vertices = edge.getGeoVertices();
        // 先在图中添加点
        for (GeoVertex i : vertices) {
            this.addVertex(i);
        }
        this.id_edge.put(edge.getId(), edge);
        this.vertices_edge.put(new ArrayList<>() {{
            add(vertices[0]);
            add(vertices[1]);
        }}, edge);
        this.vertices_edge.put(new ArrayList<>() {{
            add(vertices[1]);
            add(vertices[0]);
        }}, edge);
        vertices[0].add_conVertex(vertices[1]);
        vertices[0].addEdge(edge);
        vertices[1].addEdge(edge);
    }

    public void removeEdge(GeoEdge edge) {
        GeoVertex[] vertices = edge.getGeoVertices();
        this.id_edge.remove(edge.getId());
        vertices[0].remove_conVertex(vertices[1]);
        vertices[0].removeEdge(edge);
        vertices[1].removeEdge(edge);
    }

    /*
    public void removeEdge(GeoEdge edge) {
        GeoVertex[] vertices = edge.getGeoVertices();
        id_edge.remove(edge.getId());
        vertices[0].remove_conVertex(vertices[1]);
        vertices_edge.remove(new ArrayList<GeoVertex>() {{
            add(vertices[0]);
            add(vertices[1]);
        }});
        vertices_edge.remove(new ArrayList<GeoVertex>() {{
            add(vertices[1]);
            add(vertices[0]);
        }});
    }
    */
    public GeoVertex findVertexId(int id) {
        if (id_vertex.containsKey(id)) {
            return this.id_vertex.get(id);
        }
        return null;
    }

    public GeoEdge findEdgeId(int id) {
        if (id_edge.containsKey(id)) {
            return this.id_edge.get(id);
        }
        return null;
    }

    // 通过两节点找边
    public GeoEdge findEdgeVertices(GeoVertex vertex1, GeoVertex vertex2) {
        ArrayList<GeoVertex> temp1 = new ArrayList<>() {{
            add(vertex1);
            add(vertex2);
        }};
        ArrayList<GeoVertex> temp2 = new ArrayList<>() {{
            add(vertex2);
            add(vertex1);
        }};
        if (vertices_edge.containsKey(temp1)) {
            return vertices_edge.get(temp1);
        }
        if (vertices_edge.containsKey(temp2)) {
            return vertices_edge.get(temp2);
        }
        return null;

    }

    public int disconnectVertex(GeoVertex vertex, ArrayList<GeoVertex> vertices, int newId) {
        GeoVertex newVertex = new GeoVertex(newId, vertex.getNodeAttributes(), vertex.getCoord());
        id_vertex.put(newId, newVertex);
        if (!auxVertex.containsKey(vertex)) {
            auxVertex.put(newVertex, vertex);
            ArrayList<GeoVertex> temp = new ArrayList<>();
            temp.add(newVertex);
            subVertex.put(vertex, temp);
        } else {
            auxVertex.put(newVertex, auxVertex.get(vertex));
            GeoVertex key = auxVertex.get(newVertex);
            subVertex.get(key).add(newVertex);
        }
        for (GeoVertex v : vertices) {
            vertex.remove_conVertex(v);
            newVertex.add_conVertex(v);
        }
        newId++;
        return newId;
    }

    public void reconstructEdgeMinDeltaAngle(double k) {
        int newId = id_vertex.size() + 1;
        int judgeId = newId;
        int vertexId = 1;
        while (id_vertex.containsKey(vertexId)) {
            GeoVertex vertex = id_vertex.get(vertexId);
            ArrayList<GeoVertex> conVertices = new ArrayList<>(vertex.getConVertex());
            int count = conVertices.size();
            vertexId++;
            if (count > 2 || (count == 2 && vertexId > judgeId)) {
                HashMap<GeoVertex[], Double> calDict = new HashMap<>();
                ArrayList<GeoVertex> oldConVertices = new ArrayList<>(conVertices);
                while (!conVertices.isEmpty()) {
                    //分别计算不同组合下的夹角值
                    GeoVertex vertex1 = conVertices.get(0);
                    conVertices.remove(0);
                    for (GeoVertex v2 : conVertices) {
                        calDict.put(new GeoVertex[]{vertex1, v2}, calculateAngle(vertex1.getCoord(), v2.getCoord(), vertex.getCoord()));
                    }
                }
                List<Map.Entry<GeoVertex[], Double>> infoIds = new ArrayList<>(calDict.entrySet());
                //对value进行升序排列
                infoIds.sort((o1, o2) -> Double.compare(o1.getValue() - o2.getValue(), 0));
                //如果最小的比临界值小，则应该是这两个节点进行连接
                if (infoIds.get(0).getValue() <= k) {
                    GeoVertex[] linkedVertices = infoIds.get(0).getKey();
                    oldConVertices.remove(linkedVertices[0]);
                    oldConVertices.remove(linkedVertices[1]);
                }
                // 最小的变化角比临界值大，则该点应该作为起始节点
                else oldConVertices.remove(0);
                if (!oldConVertices.isEmpty()) newId = this.disconnectVertex(vertex, oldConVertices, newId);
            }
        }
        this.roadTrace();
    }

    public void reconstructEdgeSA(double t, double alpha) {
        int newId = id_vertex.size() + 1;
        int judgeId = newId;
        int vertexId = 1;
        HashMap<GeoVertex, GeoVertex> auxVertex = new HashMap<>();
        HashMap<GeoVertex, HashMap<ArrayList<GeoVertex>, Double>> angleDict = new HashMap<>();
        if (t <= 0.01) return;
        System.out.println("---------------正在进行初始化---------------");
        long start = System.currentTimeMillis();
        while (id_vertex.containsKey(vertexId)) {
            GeoVertex vertex = id_vertex.get(vertexId++);
            ArrayList<GeoVertex> conVertices = new ArrayList<>(vertex.getConVertex());
            int count = conVertices.size();
            if (count > 2 || (count == 2 && vertexId > judgeId)) {
                ArrayList<GeoVertex> oldConVertices = new ArrayList<>(conVertices);
                oldConVertices.remove(0);
                newId = disconnectVertex(vertex, oldConVertices, newId);
            }
        }
        for (Map.Entry<GeoVertex, ArrayList<GeoVertex>> entry : subVertex.entrySet()) {
            ArrayList<GeoVertex> chosen = new ArrayList<>(entry.getValue());
            chosen.add(entry.getKey());
            angleDict.put(entry.getKey(), new HashMap<>());
            while (!chosen.isEmpty()) {
                GeoVertex vertex1 = chosen.get(0);
                chosen.remove(0);
                for (GeoVertex vertex2 : chosen) {
                    GeoVertex conVertex1 = vertex1.getConVertex().get(0);
                    GeoVertex conVertex2 = vertex2.getConVertex().get(0);
                    double angle = Math.PI - calculateAngle(conVertex1.getCoord(), conVertex2.getCoord(), entry.getKey().getCoord());
                    angleDict.get(entry.getKey()).put(new ArrayList<>() {{
                        add(vertex1);
                        add(vertex2);
                    }}, angle / Math.PI);
                }
            }
        }
        System.out.println("初始化完成，共耗时：" + (System.currentTimeMillis() - start) + "ms");
        System.out.println("---------------正在进行退火---------------");
        HashMap<GeoVertex, ArrayList<GeoVertex>> resDict = new HashMap<>(); // 最终相邻关系表
        start = System.currentTimeMillis();
        long mid = start;
        //总损失函数
        double totalCost = Double.MAX_VALUE;
        while (t > 0.001) {
            //long mid = System.currentTimeMillis();
            System.out.println("当前温度：" + String.format("%.4f", t));
            double angleCost = 0; //角度损失
            //double disCost = 0.0; //孤立损失
            // 相连关系表
            HashMap<GeoVertex, ArrayList<GeoVertex>> conDict = new HashMap<>();
            for (GeoVertex v : id_vertex.values()) {
                conDict.put(v, new ArrayList<>(v.getConVertex()));
            }
            for (Map.Entry<GeoVertex, HashMap<ArrayList<GeoVertex>, Double>> entry : angleDict.entrySet()) {
                //当前节点上已经选择相连的节点集合
                HashSet<GeoVertex> chosenVertex = new HashSet<>();
                // 决定不相连的节点组合
                HashSet<ArrayList<GeoVertex>> disconnectVertex = new HashSet<>();
                //nor = {[v1,v2] = 0.2} 表示如果v1和v2合并的话角度收益为0.2
                HashMap<ArrayList<GeoVertex>, Double> nor = new HashMap<>(entry.getValue());
                while (!nor.isEmpty()) {
                    nor = normalize(nor);
                    double chosenNum = Math.random();
                    double nowNum = 0.0;
                    for (Map.Entry<ArrayList<GeoVertex>, Double> val : nor.entrySet()) {
                        nowNum += val.getValue();
                        if (chosenNum <= nowNum) {
                            //选中!
                            //判断是否要合并
                            double chosenNum2 = Math.random();
                            double nowAngle = entry.getValue().get(val.getKey());
                            if (chosenNum2 - 0.5 > Math.log(nowAngle + 1) || nowAngle <= 0.8) {
                                //这两点不合并，也就是它们的相邻点所组成的边不相连
                                disconnectVertex.add(val.getKey());
                                //计算孤立损失
                                //disCost = disCost + nowAngle * 2;
                            } else {
                                //要合并
                                GeoVertex vertex1 = val.getKey().get(0);
                                GeoVertex vertex2 = val.getKey().get(1);
                                chosenVertex.add(vertex1);
                                chosenVertex.add(vertex2);
                                //计算角度损失
                                angleCost = angleCost + nowAngle;
                                //合并两个节点，相当于是在删除vertex2，但要保留vertex2的相邻关系
                                GeoVertex conVertex2 = conDict.get(vertex2).get(0);
                                conDict.get(vertex1).add(conVertex2);
                                conDict.get(conVertex2).remove(vertex2);
                                conDict.get(conVertex2).add(vertex1);
                                conDict.remove(vertex2);
                            }
                            nor.clear();
                            for (Map.Entry<ArrayList<GeoVertex>, Double> temp : entry.getValue().entrySet()) {
                                if (!chosenVertex.contains(temp.getKey().get(0)) && !chosenVertex.contains(temp.getKey().get(1)) && !disconnectVertex.contains(temp.getKey())) {
                                    nor.put(temp.getKey(), temp.getValue());
                                }
                            }
                            break;
                        }
                    }
                }
            }
            //计算路径损失
            double roadCost = calRoadCost(conDict); //路径损失
            double nowCost = (Math.pow(roadCost, 2.0) + angleCost); //+ disCost
            if (totalCost > nowCost) { // || Math.random() < Math.exp((nowCost - totalCost) / t)
                totalCost = nowCost;
                resDict = new HashMap<>(conDict);
                System.out.println("接受当前目标函数值：" + totalCost);
                System.out.println("" + (System.currentTimeMillis() - mid) + "ms");
                mid = System.currentTimeMillis();
            }
            t *= alpha;
        }
        System.out.println("退火完成，共耗时：" + (System.currentTimeMillis() - start) + "ms");
        System.out.println("---------------正在重建相邻关系---------------");
        start = System.currentTimeMillis();
        // 根据最终相邻关系表重建相邻关系
        ArrayList<GeoVertex> delVertices = new ArrayList<>();
        for (Map.Entry<Integer, GeoVertex> entry : id_vertex.entrySet()) {
            if (!resDict.containsKey(entry.getValue())) delVertices.add(entry.getValue());
        }
        for (GeoVertex v : delVertices) id_vertex.remove(v.getId());
        for (Map.Entry<GeoVertex, ArrayList<GeoVertex>> entry : resDict.entrySet()) {
            ArrayList<GeoVertex> conVertices = new ArrayList<>();
            for (GeoVertex v : entry.getValue()) {
                conVertices.add(id_vertex.get(v.getId()));
            }
            id_vertex.get(entry.getKey().getId()).setConVertex(conVertices);
        }
        roadTrace();
        System.out.println("重建相邻关系完成，共耗时：" + (System.currentTimeMillis() - start) + "ms");
    }

    public void joinByMaxAngelOfEdgeAndSegment() {
        //如果这两条边以及它们对应的路段的夹角最大的话，则一定相连(可选)
        //需要重新构建node节点
        long start = System.currentTimeMillis();
        System.out.println("---------------开始根据边以及路段重建node---------------");
        //构建哈希表方便查询两个路段（segment）和边（edge）之间的角度值
        HashMap<GeoVertex, HashMap<ArrayList<GeoEdge>, Double>> angleEdgeVertex = new HashMap<>();
        HashMap<GeoVertex, HashMap<ArrayList<GeoSegment>, Double>> angleSegmentVertex = new HashMap<>();
        for (GeoVertex v : node) {
            HashMap<ArrayList<GeoEdge>, Double> angleEdge = new HashMap<>();
            HashMap<ArrayList<GeoSegment>, Double> angleSegment = new HashMap<>();
            ArrayList<GeoSegment> segments1 = new ArrayList<>(v.getSegments());
            for (GeoSegment calSeg1 : segments1) {
                ArrayList<GeoVertex> vs1 = calSeg1.getVertices();
                ArrayList<GeoSegment> segments2 = new ArrayList<>(segments1);
                segments2.remove(calSeg1);
                GeoVertex v1 = vs1.get(0);
                GeoVertex v2 = vs1.get(vs1.size() - 1);
                assert v1 == v : "wtf1 " + v + " " + v1;
                assert v1 != v2 : "wtf2 " + v1 + " " + v2;
                for (GeoSegment calSeg2 : segments2) {
                    ArrayList<GeoVertex> vs2 = calSeg2.getVertices();
                    GeoVertex v3 = vs2.get(0);
                    GeoVertex v4 = vs2.get(vs2.size() - 1);
                    assert v3 == v : "wtf3 " + v + " " + v3;
                    assert v3 != v4 : "wtf4 " + v3 + " " + v4 + " " + vs2;
                    double angle = calculateAngle(v4.getCoord(), v2.getCoord(), v.getCoord());
                    angleSegment.put(new ArrayList<>() {{
                        add(calSeg1);
                        add(calSeg2);
                    }}, angle);
                }
            }
            ArrayList<GeoEdge> edges1 = new ArrayList<>(v.getEdges());
            for (GeoEdge e1 : edges1) {
                GeoVertex v1;
                ArrayList<GeoEdge> edges2 = new ArrayList<>(edges1);
                edges2.remove(e1);
                if (v == e1.getGeoVertices()[0]) {
                    v1 = e1.getGeoVertices()[1];
                } else {
                    v1 = e1.getGeoVertices()[0];
                    assert v == e1.getGeoVertices()[1] : "wtf5 " + v + " " + Arrays.toString(e1.getGeoVertices());
                }
                for (GeoEdge e2 : edges2) {
                    GeoVertex v2;
                    if (v == e2.getGeoVertices()[0]) {
                        v2 = e2.getGeoVertices()[1];
                    } else {
                        v2 = e2.getGeoVertices()[0];
                        assert v == e2.getGeoVertices()[1] : "wtf6 " + v + " " + Arrays.toString(e2.getGeoVertices());
                    }
                    double angle = calculateAngle(v1.getCoord(), v2.getCoord(), v.getCoord());
                    angleEdge.put(new ArrayList<>() {{
                        add(e1);
                        add(e2);
                    }}, angle);
                }

            }
            angleEdgeVertex.put(v, angleEdge);
            angleSegmentVertex.put(v, angleSegment);
        }
        int k = 1;
        boolean flag = true;
        while (flag) {
            flag = false;
            System.out.println("第" + k++ + "次");
            for (GeoVertex vertex : node) {
                List<Map.Entry<ArrayList<GeoEdge>, Double>> edgeList = new ArrayList<>(angleEdgeVertex.get(vertex).entrySet());
                List<Map.Entry<ArrayList<GeoSegment>, Double>> segmentList = new ArrayList<>(angleSegmentVertex.get(vertex).entrySet());
                //按value由小到大排序
                edgeList.sort(Map.Entry.comparingByValue());
                segmentList.sort(Map.Entry.comparingByValue());
                ArrayList<GeoEdge> targetEdge = edgeList.get(0).getKey();
                ArrayList<GeoSegment> targetSegment = segmentList.get(0).getKey();
                boolean judge = (targetSegment.get(0).getEdges().contains(targetEdge.get(0)) && targetSegment.get(1).getEdges().contains(targetEdge.get(1))) ||
                        (targetSegment.get(0).getEdges().contains(targetEdge.get(1)) && targetSegment.get(1).getEdges().contains(targetEdge.get(0)));
                if (judge) {
                    flag = true;
                    //这两条边以及它们对应的路段的夹角最大,将他们相连
                    GeoVertex conVertex1 = targetEdge.get(0).getGeoVertices()[0] == vertex ? targetEdge.get(0).getGeoVertices()[1] : targetEdge.get(0).getGeoVertices()[0];
                    GeoVertex conVertex2 = targetEdge.get(1).getGeoVertices()[0] == vertex ? targetEdge.get(1).getGeoVertices()[1] : targetEdge.get(1).getGeoVertices()[0];
                    assert targetEdge.get(0).getGeoVertices()[0] == vertex || targetEdge.get(0).getGeoVertices()[1] == vertex : "wtf1";
                    assert targetEdge.get(1).getGeoVertices()[0] == vertex || targetEdge.get(1).getGeoVertices()[1] == vertex : "wtf2";
                    GeoSegment segment1 = targetSegment.get(0);
                    GeoSegment segment2 = targetSegment.get(1);
                    GeoVertex newVertex = new GeoVertex(id_vertex.size() + 1, new HashMap<>(), vertex.getCoord());
                    //修改边对应的点
                    GeoEdge edge1 = targetEdge.get(0);
                    GeoEdge edge2 = targetEdge.get(1);
                    //删除旧边，添加新边
                    this.removeEdge(edge1);
                    this.removeEdge(edge2);
                    edge1.setGeoVertices(new GeoVertex[]{newVertex, conVertex1});
                    edge2.setGeoVertices(new GeoVertex[]{newVertex, conVertex2});
                    this.addEdge(edge1);
                    this.addEdge(edge2);

                    ArrayList<GeoEdge> newEdges1 = new ArrayList<>();
                    ArrayList<GeoEdge> newEdges2 = new ArrayList<>();
                    ArrayList<GeoVertex> newVertices1 = new ArrayList<>();
                    ArrayList<GeoVertex> newVertices2 = new ArrayList<>();

                    GeoVertex vertex1 = segment1.getVertices().get(segment1.getVertices().size() - 1);
                    GeoVertex vertex2 = segment2.getVertices().get(segment2.getVertices().size() - 1);

                    assert segment1.getVertices().get(0) == segment2.getVertices().get(0) : "wfffffff";
                    for (int i = segment1.getEdges().size() - 1; i >= 0; i--) {
                        newEdges1.add(segment1.getEdges().get(i));
                    }
                    newEdges1.addAll(segment2.getEdges());
                    for (int i = segment2.getEdges().size() - 1; i >= 0; i--) {
                        newEdges2.add(segment2.getEdges().get(i));
                    }
                    newEdges2.addAll(segment1.getEdges());

                    for (int i = segment1.getVertices().size() - 1; i >= 0; i--) {
                        newVertices1.add(segment1.getVertices().get(i));
                    }
                    newVertices1.remove(vertex);
                    newVertices1.addAll(segment2.getVertices());
                    for (int i = segment2.getVertices().size() - 1; i >= 0; i--) {
                        newVertices2.add(segment2.getVertices().get(i));
                    }
                    newVertices2.remove(vertex);
                    newVertices2.addAll(segment1.getVertices());
                    GeoSegment newSegment1 = new GeoSegment(id_segment.size() + 1, newEdges1, newVertices1);
                    id_segment.put(newSegment1.getId(), newSegment1);
                    GeoSegment newSegment2 = new GeoSegment(id_segment.size() + 1, newEdges2, newVertices2);
                    id_segment.put(newSegment2.getId(), newSegment2);
                    vertex.removeSegment(segment1);
                    vertex.removeSegment(segment2);
                    //更新segment角度表以及edge角度表
                    ArrayList<ArrayList<GeoSegment>> tarSeg = new ArrayList<>(angleSegmentVertex.get(vertex).keySet());
                    ArrayList<ArrayList<GeoSegment>> rmSeg = new ArrayList<>();
                    for (ArrayList<GeoSegment> key : tarSeg) {
                        if (key.contains(segment1) || key.contains(segment2)) rmSeg.add(key);
                    }
                    for (ArrayList<GeoSegment> key : rmSeg) {
                        angleSegmentVertex.get(vertex).remove(key);
                    }

                    ArrayList<ArrayList<GeoEdge>> tarEdge = new ArrayList<>(angleEdgeVertex.get(vertex).keySet());
                    ArrayList<ArrayList<GeoEdge>> rmEdge = new ArrayList<>();
                    for (ArrayList<GeoEdge> key : tarEdge) {
                        if (key.contains(edge1) || key.contains(edge2)) rmEdge.add(key);
                    }
                    for (ArrayList<GeoEdge> key : rmEdge) {
                        angleEdgeVertex.get(vertex).remove(key);
                    }

                    if (newSegment1.getVertices().get(0) == newSegment1.getVertices().get(newSegment1.getVertices().size() - 1)) {
                        // 出现了环路
                        ArrayList<GeoVertex> removeVertices1 = new ArrayList<>(newSegment1.getVertices());
                        circleRoad.add(newSegment1);
                        roadDict.add(newSegment1.getVertices());
                        removeVertices1.remove(removeVertices1.size() - 1);
                        removeVertices1.remove(0);
                        for (GeoVertex rrv : removeVertices1) this.removeVertex(rrv);
                    } else {
                        if (vertex1.getConVertex().size() > 2) {
                            ArrayList<ArrayList<GeoSegment>> tarSeg1 = new ArrayList<>(angleSegmentVertex.get(vertex1).keySet());
                            ArrayList<ArrayList<GeoSegment>> rmSeg1 = new ArrayList<>();
                            for (ArrayList<GeoSegment> key : tarSeg1) {
                                if (key.contains(segment1)) rmSeg1.add(key);
                            }
                            for (ArrayList<GeoSegment> key : rmSeg1) {
                                angleSegmentVertex.get(vertex).remove(key);
                            }
                            vertex1.removeSegment(segment1);
                            vertex1.addSegment(newSegment1);
                            for (GeoSegment key : vertex1.getSegments()) {
                                GeoVertex tarS = key.getVertices().get(key.getVertices().size() - 1);
                                double angle = calculateAngle(tarS.getCoord(), vertex2.getCoord(), vertex1.getCoord());
                                angleSegmentVertex.get(vertex1).put(new ArrayList<>() {{
                                    add(key);
                                    add(newSegment1);
                                }}, angle);
                            }
                        }
                    }
                    if (newSegment2.getVertices().get(0) == newSegment2.getVertices().get(newSegment2.getVertices().size() - 1)) {
                        // 出现了环路
                        ArrayList<GeoVertex> removeVertices2 = new ArrayList<>(newSegment2.getVertices());
                        circleRoad.add(newSegment2);
                        roadDict.add(newSegment2.getVertices());
                        removeVertices2.remove(removeVertices2.size() - 1);
                        removeVertices2.remove(0);
                        for (GeoVertex rrv : removeVertices2) this.removeVertex(rrv);
                    } else {
                        if (vertex2.getConVertex().size() > 2) {
                            ArrayList<ArrayList<GeoSegment>> tarSeg2 = new ArrayList<>(angleSegmentVertex.get(vertex2).keySet());
                            ArrayList<ArrayList<GeoSegment>> rmSeg2 = new ArrayList<>();
                            for (ArrayList<GeoSegment> key : tarSeg2) {
                                if (key.contains(segment2)) rmSeg2.add(key);
                            }
                            for (ArrayList<GeoSegment> key : rmSeg2) {
                                angleSegmentVertex.get(vertex).remove(key);
                            }
                            vertex2.removeSegment(segment2);
                            vertex2.addSegment(newSegment2);
                            for (GeoSegment key : vertex2.getSegments()) {
                                GeoVertex tarS = key.getVertices().get(key.getVertices().size() - 1);
                                double angle = calculateAngle(tarS.getCoord(), vertex1.getCoord(), vertex2.getCoord());
                                angleSegmentVertex.get(vertex2).put(new ArrayList<>() {{
                                    add(key);
                                    add(newSegment2);
                                }}, angle);
                            }
                        }
                    }
                }
            }
            //重新构建路段(两个node之间的边的集合) node的度大于2
            node.clear();
            for (GeoVertex vv : id_vertex.values()) {
                if (vv.getConVertex().size() > 2) node.add(vv);
            }
            System.out.println(node.size());
        }
        ArrayList<GeoVertex> rmVertices = new ArrayList<>();
        for(GeoVertex v: id_vertex.values()){
            if(v.getConVertex().size()==0) rmVertices.add(v);
        }
        for(GeoVertex v: rmVertices){
            this.removeVertex(v);
        }
        System.out.println("---------------重建node完成---------------");
        System.out.println("重建node完成，共耗时：" + (System.currentTimeMillis() - start) + "ms");

    }

    public void strokeBuildingSA(double t, double alpha) {
        /*使用模拟退火算法重建图的顶点(node)之间的新的相邻边关系
        Args:
        t（double）：模拟退火算法的初始温度。
        alpha（double）：模拟退火算法的温度衰减率。
        Returns:
        构建图中的道路二维数组。每个数组是节点集合。*/
        //如果这两条边以及它们对应的路段的夹角最大的话，则相连，重建node
        long start = System.currentTimeMillis();
        System.out.println("---------------正在进行初始化---------------");
        joinByMaxAngelOfEdgeAndSegment();
        //----TODO----
        int newId = Collections.max(id_vertex.keySet()) + 1;
        if (t <= 0.01) return;
        //构建哈希表方便查询两个路段（segment）和边（edge）之间的角度值
        HashMap<GeoVertex, HashMap<ArrayList<GeoEdge>, Double>> angleEdgeVertex = new HashMap<>();
        HashMap<GeoVertex, HashMap<ArrayList<GeoSegment>, Double>> angleSegmentVertex = new HashMap<>();
        for (GeoVertex v : node) {
            HashMap<ArrayList<GeoEdge>, Double> angleEdge = new HashMap<>();
            HashMap<ArrayList<GeoSegment>, Double> angleSegment = new HashMap<>();
            ArrayList<GeoSegment> segments1 = new ArrayList<>(v.getSegments());
            for (GeoSegment calSeg1 : segments1) {
                ArrayList<GeoVertex> vs1 = calSeg1.getVertices();
                ArrayList<GeoSegment> segments2 = new ArrayList<>(segments1);
                segments2.remove(calSeg1);
                GeoVertex v1 = vs1.get(0);
                GeoVertex v2 = vs1.get(vs1.size() - 1);
                assert v1 == v : "wtf1 " + v + " " + v1;
                assert v1 != v2 : "wtf2 " + v1 + " " + v2;
                for (GeoSegment calSeg2 : segments2) {
                    ArrayList<GeoVertex> vs2 = calSeg2.getVertices();
                    GeoVertex v3 = vs2.get(0);
                    GeoVertex v4 = vs2.get(vs2.size() - 1);
                    assert v3 == v : "wtf3 " + v + " " + v3;
                    assert v3 != v4 : "wtf4 " + v3 + " " + v4 + " " + vs2;
                    double angle = calculateAngle(v4.getCoord(), v2.getCoord(), v.getCoord());
                    angleSegment.put(new ArrayList<>() {{
                        add(calSeg1);
                        add(calSeg2);
                    }}, angle);
                }
            }
            ArrayList<GeoEdge> edges1 = new ArrayList<>(v.getEdges());
            for (GeoEdge e1 : edges1) {
                GeoVertex v1;
                ArrayList<GeoEdge> edges2 = new ArrayList<>(edges1);
                edges2.remove(e1);
                if (v == e1.getGeoVertices()[0]) {
                    v1 = e1.getGeoVertices()[1];
                } else {
                    v1 = e1.getGeoVertices()[0];
                    assert v == e1.getGeoVertices()[1] : "wtf5 " + v + " " + Arrays.toString(e1.getGeoVertices());
                }
                for (GeoEdge e2 : edges2) {
                    GeoVertex v2;
                    if (v == e2.getGeoVertices()[0]) {
                        v2 = e2.getGeoVertices()[1];
                    } else {
                        v2 = e2.getGeoVertices()[0];
                        assert v == e2.getGeoVertices()[1] : "wtf6 " + v + " " + Arrays.toString(e2.getGeoVertices());
                    }
                    double angle = calculateAngle(v1.getCoord(), v2.getCoord(), v.getCoord());
                    angleEdge.put(new ArrayList<>() {{
                        add(e1);
                        add(e2);
                    }}, angle);
                }

            }
            angleEdgeVertex.put(v, angleEdge);
            angleSegmentVertex.put(v, angleSegment);
        }
        System.out.println("---------------初始化完成---------------");
        System.out.println("初始化完成，共耗时：" + (System.currentTimeMillis() - start) + "ms");
        System.out.println("---------------正在进行退火---------------");
        HashMap<GeoVertex, ArrayList<GeoVertex>> resDict = new HashMap<>(); // 最终相邻关系表
        start = System.currentTimeMillis();
        long mid = start;
        //总损失函数
        double totalCost = Double.MAX_VALUE;
        /*while (t > 0.001) {
            //long mid = System.currentTimeMillis();
            System.out.println("当前温度：" + String.format("%.4f", t));
            double angleCost = 0; //角度损失
            //double disCost = 0.0; //孤立损失
            // 相连关系表
            HashMap<GeoVertex, ArrayList<GeoVertex>> conDict = new HashMap<>();
            for (GeoVertex v : id_vertex.values()) {
                conDict.put(v, new ArrayList<>(v.getConVertex()));
            }
            for (Map.Entry<GeoVertex, HashMap<ArrayList<GeoVertex>, Double>> entry : angleDict.entrySet()) {
                //当前节点上已经选择相连的节点集合
                HashSet<GeoVertex> chosenVertex = new HashSet<>();
                // 决定不相连的节点组合
                HashSet<ArrayList<GeoVertex>> disconnectVertex = new HashSet<>();
                //nor = {[v1,v2] = 0.2} 表示如果v1和v2合并的话角度收益为0.2
                HashMap<ArrayList<GeoVertex>, Double> nor = new HashMap<>(entry.getValue());
                while (!nor.isEmpty()) {
                    nor = normalize(nor);
                    double chosenNum = Math.random();
                    double nowNum = 0.0;
                    for (Map.Entry<ArrayList<GeoVertex>, Double> val : nor.entrySet()) {
                        nowNum += val.getValue();
                        if (chosenNum <= nowNum) {
                            //选中!
                            //判断是否要合并
                            double chosenNum2 = Math.random();
                            double nowAngle = entry.getValue().get(val.getKey());
                            if (chosenNum2 - 0.5 > Math.log(nowAngle + 1) || nowAngle <= 0.8) {
                                //这两点不合并，也就是它们的相邻点所组成的边不相连
                                disconnectVertex.add(val.getKey());
                                //计算孤立损失
                                //disCost = disCost + nowAngle * 2;
                            } else {
                                //要合并
                                GeoVertex vertex1 = val.getKey().get(0);
                                GeoVertex vertex2 = val.getKey().get(1);
                                chosenVertex.add(vertex1);
                                chosenVertex.add(vertex2);
                                //计算角度损失
                                angleCost = angleCost + nowAngle;
                                //合并两个节点，相当于是在删除vertex2，但要保留vertex2的相邻关系
                                GeoVertex conVertex2 = conDict.get(vertex2).get(0);
                                conDict.get(vertex1).add(conVertex2);
                                conDict.get(conVertex2).remove(vertex2);
                                conDict.get(conVertex2).add(vertex1);
                                conDict.remove(vertex2);
                            }
                            nor.clear();
                            for (Map.Entry<ArrayList<GeoVertex>, Double> temp : entry.getValue().entrySet()) {
                                if (!chosenVertex.contains(temp.getKey().get(0)) && !chosenVertex.contains(temp.getKey().get(1)) && !disconnectVertex.contains(temp.getKey())) {
                                    nor.put(temp.getKey(), temp.getValue());
                                }
                            }
                            break;
                        }
                    }
                }
            }
            //计算路径损失
            double roadCost = calRoadCost(conDict); //路径损失
            double nowCost = (Math.pow(roadCost, 2.0) + angleCost); //+ disCost
            if (totalCost > nowCost) { // || Math.random() < Math.exp((nowCost - totalCost) / t)
                totalCost = nowCost;
                resDict = new HashMap<>(conDict);
                System.out.println("接受当前目标函数值：" + totalCost);
                System.out.println("" + (System.currentTimeMillis() - mid) + "ms");
                mid = System.currentTimeMillis();
            }
            t *= alpha;
        }
        System.out.println("退火完成，共耗时：" + (System.currentTimeMillis() - start) + "ms");
        System.out.println("---------------正在重建相邻关系---------------");
        start = System.currentTimeMillis();
        // 根据最终相邻关系表重建相邻关系
        ArrayList<GeoVertex> delVertices = new ArrayList<>();
        for (Map.Entry<Integer, GeoVertex> entry : id_vertex.entrySet()) {
            if (!resDict.containsKey(entry.getValue())) delVertices.add(entry.getValue());
        }
        for (GeoVertex v : delVertices) id_vertex.remove(v.getId());
        for (Map.Entry<GeoVertex, ArrayList<GeoVertex>> entry : resDict.entrySet()) {
            ArrayList<GeoVertex> conVertices = new ArrayList<>();
            for (GeoVertex v : entry.getValue()) {
                conVertices.add(id_vertex.get(v.getId()));
            }
            id_vertex.get(entry.getKey().getId()).setConVertex(conVertices);
        }
        roadTrace();
        System.out.println("重建相邻关系完成，共耗时：" + (System.currentTimeMillis() - start) + "ms");*/
    }

    public void roadTrace() {
        int road_count = roadDict.size();
        ArrayList<GeoVertex> nodeVertex = new ArrayList<>();
        HashMap<GeoVertex, Boolean> visitedVertex = new HashMap<>();
        for (GeoVertex vertex : this.id_vertex.values()) {
            visitedVertex.put(vertex, false);
            if (vertex.getConVertex().size() == 1) {
                nodeVertex.add(vertex);
            } else if (vertex.getConVertex().size() != 2) {
                visitedVertex.put(vertex, true);
            }
        }
        //端点
        for (GeoVertex vertex : nodeVertex) {
            if (!visitedVertex.get(vertex)) {
                visitedVertex.put(vertex, true);
                roadDict.add(new ArrayList<>());
                GeoVertex nowVertex = vertex;
                GeoVertex nextVertex = vertex.getConVertex().get(0);
                roadDict.get(road_count).add(nowVertex);
                while (nextVertex.getConVertex().size() != 1) {
                    GeoVertex preVertex = nowVertex;
                    nowVertex = nextVertex;
                    ArrayList<GeoVertex> vertices = new ArrayList<>(nowVertex.getConVertex());
                    if (vertices.get(0) == preVertex) nextVertex = vertices.get(1);
                    else nextVertex = vertices.get(0);
                    visitedVertex.put(nowVertex, true);
                    roadDict.get(road_count).add(nowVertex);
                }
                roadDict.get(road_count).add(nextVertex);
                visitedVertex.put(nextVertex, true);
                road_count++;
            }
        }
        //环路
        for (GeoVertex vertex : this.id_vertex.values()) {
            if (visitedVertex.get(vertex)) continue;
            visitedVertex.put(vertex, true);
            roadDict.add(new ArrayList<>());
            GeoVertex nextVertex = vertex.getConVertex().get(0);
            roadDict.get(road_count).add(vertex);
            GeoVertex nowVertex = vertex;
            while (!visitedVertex.get(nextVertex)) {
                GeoVertex preVertex = nowVertex;
                nowVertex = nextVertex;
                ArrayList<GeoVertex> vertices = new ArrayList<>(nowVertex.getConVertex());
                if (vertices.get(0) == preVertex) nextVertex = vertices.get(1);
                else nextVertex = vertices.get(0);
                visitedVertex.put(nowVertex, true);
                roadDict.get(road_count).add(nowVertex);
            }
            roadDict.get(road_count).add(nextVertex);
            road_count += 1;
        }
    }

    public double calRoadCost(HashMap<GeoVertex, ArrayList<GeoVertex>> conDict) {
        int road_count = 1;
        ArrayList<GeoVertex> nodeVertex = new ArrayList<>();
        HashMap<GeoVertex, Boolean> visitedVertex = new HashMap<>();
        for (Map.Entry<GeoVertex, ArrayList<GeoVertex>> entry : conDict.entrySet()) {
            visitedVertex.put(entry.getKey(), false);
            if (entry.getValue().size() == 1) {
                nodeVertex.add(entry.getKey());
            } else if (entry.getValue().size() != 2) {
                visitedVertex.put(entry.getKey(), true);
            }
        }
        //端点
        for (GeoVertex vertex : nodeVertex) {
            if (!visitedVertex.get(vertex)) {
                visitedVertex.put(vertex, true);
                GeoVertex nowVertex = vertex;
                GeoVertex nextVertex = conDict.get(vertex).get(0);
                while (conDict.get(nextVertex).size() != 1) {
                    GeoVertex preVertex = nowVertex;
                    nowVertex = nextVertex;
                    ArrayList<GeoVertex> vertices = new ArrayList<>(conDict.get(nowVertex));
                    if (vertices.get(0) == preVertex) nextVertex = vertices.get(1);
                    else nextVertex = vertices.get(0);
                    visitedVertex.put(nowVertex, true);
                }
                visitedVertex.put(nextVertex, true);
                road_count++;
            }
        }
        //环路
        for (GeoVertex vertex : conDict.keySet()) {
            if (visitedVertex.get(vertex)) continue;
            visitedVertex.put(vertex, true);
            GeoVertex nextVertex = conDict.get(vertex).get(0);
            GeoVertex nowVertex = vertex;
            while (!visitedVertex.get(nextVertex)) {
                GeoVertex preVertex = nowVertex;
                nowVertex = nextVertex;
                ArrayList<GeoVertex> vertices = new ArrayList<>(conDict.get(nowVertex));
                if (vertices.get(0) == preVertex) nextVertex = vertices.get(1);
                else nextVertex = vertices.get(0);
                visitedVertex.put(nowVertex, true);
            }
            road_count += 1;
        }
        return road_count;
    }

    public void outputGeographRoads(String outputPath) throws Exception {
        File file = new File(outputPath);
        if (!file.createNewFile()) {
            return;
        }
        if (roadDict == null) return;
        GeometryFactory geometryFactory = new GeometryFactory();
        Map<String, Serializable> params = new HashMap<>();
        ShapefileDataStore ds;
        params.put(ShapefileDataStoreFactory.URLP.key, file.toURI().toURL());
        params.put("create spatial index", Boolean.TRUE);
        ds = (ShapefileDataStore) new ShapefileDataStoreFactory().createNewDataStore(params);
        // 定义图形信息和属性信息
        SimpleFeatureTypeBuilder tb = new SimpleFeatureTypeBuilder();
        CoordinateReferenceSystem crsTarget = CRS.parseWKT(this.proj);
        tb.setCRS(crsTarget);
        tb.setName(this.name);
        //几何类型
        tb.add("the_geom", LineString.class);
        Map.Entry<Integer, GeoEdge> attSchema = id_edge.entrySet().iterator().next();
        //字段名
        for (Map.Entry<String, String> item : attSchema.getValue().getEdgeAttribute().entrySet()) {
            tb.add(item.getKey(), String.class);
        }
        tb.add("id", String.class);
        ds.createSchema(tb.buildFeatureType());
        ds.setCharset(StandardCharsets.UTF_8);
        //写入shapefile
        FeatureWriter<SimpleFeatureType, SimpleFeature> writer = ds.getFeatureWriter(ds.getTypeNames()[0], AUTO_COMMIT);
        int nowId = 1;
        for (ArrayList<GeoVertex> item : roadDict) {
            //获取当前写对象
            SimpleFeature feature = writer.next();
            //几何坐标，LineString是通过Coordinate[]构建的，同理MultiLineString是通过LineString[]构建
            Coordinate[] coordinates = new Coordinate[item.size()];
            //获取线段的属性集合，方便后面求众数、平均数等。
            HashMap<String, ArrayList<String>> att = new HashMap<>();
            int i = 0;
            GeoVertex preVertex = null;
            for (GeoVertex v : item) {
                coordinates[i++] = new Coordinate(v.getCoord()[0], v.getCoord()[1]);
                if (preVertex != null) {
                    v = auxVertex.getOrDefault(v, v);
                    preVertex = auxVertex.getOrDefault(preVertex, preVertex);
                    GeoEdge nowEdge = this.findEdgeVertices(preVertex, v);
                    assert nowEdge != null: " " + preVertex + " " + v;
                    HashMap<String, String> edgeAttribute = nowEdge.getEdgeAttribute();
                    for (Map.Entry<String, String> entry : edgeAttribute.entrySet()) {
                        if (att.containsKey(entry.getKey())) att.get(entry.getKey()).add(entry.getValue());
                        else att.put(entry.getKey(), new ArrayList<>() {{
                            add(entry.getValue());
                        }});
                    }
                }
                preVertex = v;
            }
            HashMap<String, String> attr = new HashMap<>();
            LineString lineString = geometryFactory.createLineString(coordinates);
            //设置几何坐标
            feature.setAttribute("the_geom", lineString);
            //获取出现次数最多的字符串作为属性字段的值（求众数）
            for (Map.Entry<String, ArrayList<String>> entry : att.entrySet()) {
                HashMap<String, Integer> temp = new HashMap<>();
                for (String s : entry.getValue()) {
                    if (temp.containsKey(s)) {
                        int i1 = temp.get(s);
                        temp.put(s, ++i1);
                    } else temp.put(s, 1);
                }
                int max = 0;
                String s = null;
                for (Map.Entry<String, Integer> t : temp.entrySet()) {
                    String tempStr = t.getKey().strip();
                    if (max < t.getValue() && !"".equals(tempStr)) {
                        s = tempStr;
                        max = t.getValue();
                    }
                }
                attr.put(entry.getKey(), s);
            }
            //设置属性
            for (Map.Entry<String, String> entry : attr.entrySet()) {
                feature.setAttribute(entry.getKey(), entry.getValue());
            }
            feature.setAttribute("id", nowId++);
        }
        writer.write();
        writer.close();
        ds.dispose();
        FileWriter fw = new FileWriter(outputPath.replace(".shp", ".prj"));
        fw.write(this.proj);
        fw.close();
    }

    public static double getDistance(double lat1, double lon1, double lat2, double lon2) {
        // 84坐标系构造GeodeticCalculator
        GeodeticCalculator geodeticCalculator = new GeodeticCalculator(DefaultGeographicCRS.WGS84);
        // 起点经纬度
        geodeticCalculator.setStartingGeographicPoint(lon1, lat1);
        // 末点经纬度
        geodeticCalculator.setDestinationGeographicPoint(lon2, lat2);
        // 计算距离，单位：米
        return geodeticCalculator.getOrthodromicDistance();
    }

    //归一化
    public HashMap<ArrayList<GeoVertex>, Double> normalize(HashMap<ArrayList<GeoVertex>, Double> hashMap) {
        double sum = 0.0;
        HashMap<ArrayList<GeoVertex>, Double> res = new HashMap<>(hashMap);
        res.replaceAll((e, v) -> Math.exp(res.get(e)) - 1);
        for (Map.Entry<ArrayList<GeoVertex>, Double> entry : res.entrySet()) {
            sum += entry.getValue();
        }
        double finalSum = sum;
        res.replaceAll((e, v) -> res.get(e) / finalSum);
        return res;
    }

    public static double calculateAngle(double[] coord1, double[] coord2, double[] mid) {
        if (coord1.length != 2 || coord2.length != 2 || mid.length != 2) return 0.0;
        double[] a = new double[]{coord1[0] - mid[0], coord1[1] - mid[1]};
        double[] b = new double[]{coord2[0] - mid[0], coord2[1] - mid[1]};
        double dot_product = 0.0;
        for (int i = 0; i < a.length; i++) dot_product += (a[i] * b[i]);
        double a_norm = 0.0;
        for (double x : a) a_norm += (x * x);
        a_norm = Math.sqrt(a_norm);
        double b_norm = 0.0;
        for (double y : b) b_norm += (y * y);
        b_norm = Math.sqrt(b_norm);
        double cos_theta = dot_product / (a_norm * b_norm);
        cos_theta = Double.parseDouble(String.format("%.8f", cos_theta));
        double theta = Math.acos(cos_theta);
        //if (theta > Math.PI) theta = 2 * Math.PI - theta;

        return Math.PI - theta;
    }
}
