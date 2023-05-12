package lyd.SA;

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
    private final HashMap<List<GeoVertex>, GeoEdge> vertices_edge;  // 点号对应的边号字典 {(vertex1, vertex2): edge1, (vertex3, vertex4): edge2,...}
    private String name; // 图名
    private final HashMap<List<Double>, GeoVertex> vertexCoord; // 点的坐标对应的点号
    private String proj; //坐标系字符串

    private HashMap<Integer, ArrayList<GeoVertex>> roadDict; // 当前图中路的路径点几何

    public GeoGragh(String name) {
        this.id_vertex = new HashMap<>();
        this.id_edge = new HashMap<>();
        this.vertices_edge = new HashMap<>();
        this.name = name;
        this.vertexCoord = new HashMap<>();
        this.proj = "";
        this.roadDict = null;
    }

    //通过边文件路径名以及图名构建该图
    public GeoGragh(String name, String path) throws Exception {
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
                            // 通过前置节点和现在的节点构建该边
                            GeoEdge nowEdge = new GeoEdge(edgeId++, preVertex, nowVertex, new HashMap<>());
                            geoEdges.add(nowEdge);
                            this.addEdge(nowEdge);
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
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getProj() {
        return proj;
    }

    public void setProj(String proj) {
        this.proj = proj;
    }

    public HashMap<Integer, ArrayList<GeoVertex>> getRoadDict() {
        return roadDict;
    }

    @Override
    public String toString() {
        return "图名：" + name + "坐标系字符串：" + proj;
    }

    public void addVertex(GeoVertex vertex) {
        if (this.id_vertex.containsKey(vertex.getId())) return;
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
            vertices_edge.remove(new ArrayList<GeoVertex>() {{
                add(vertex);
                add(v);
            }});
            vertices_edge.remove(new ArrayList<GeoVertex>() {{
                add(v);
                add(vertex);
            }});
            ArrayList<Double> t = new ArrayList<>();
            t.add(vertex.getCoord()[0]);
            t.add(vertex.getCoord()[1]);
            vertexCoord.remove(t);
            this.id_edge.remove(edge.getId());
        }
    }

    public void addEdge(GeoEdge edge) {
        if (this.id_edge.containsKey(edge.getId())) return;
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
        vertices[0].add_conVertex(vertices[1]);
    }

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

    public int disconnectVertex(GeoVertex vertex, ArrayList<GeoVertex> vertices, int newId, HashMap<GeoVertex, ArrayList<GeoVertex>> subVertex, HashMap<GeoVertex, GeoVertex> auxVertex) {
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
            GeoEdge nowEdge = findEdgeVertices(v, vertex);
            vertices_edge.put(new ArrayList<>() {{
                add(newVertex);
                add(v);
            }}, nowEdge);
            vertex.remove_conVertex(v);
            newVertex.add_conVertex(v);
        }
        newId++;
        return newId;
    }

    public void connectVertex(GeoVertex vertex1, GeoVertex vertex2) {
        ArrayList<GeoVertex> conVertices2 = vertex2.getConVertex();
        for (GeoVertex v : conVertices2) {
            vertex1.add_conVertex(v);
            vertex2.remove_conVertex(v);
        }
    }

    public void reconstructEdgeMinDeltaAngle(double k) {
        int newId = id_vertex.size() + 1;
        int judgeId = newId;
        int vertexId = 1;
        HashMap<GeoVertex, ArrayList<GeoVertex>> subVertex = new HashMap<>();
        HashMap<GeoVertex, GeoVertex> auxVertex = new HashMap<>();
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
                if (!oldConVertices.isEmpty())
                    newId = this.disconnectVertex(vertex, oldConVertices, newId, subVertex, auxVertex);
            }
        }
        roadDict = this.roadTrace();
    }

    public void reconstructEdgeSA(double t, double alpha) {
        int newId = id_vertex.size() + 1;
        int judgeId = newId;
        int vertexId = 1;
        HashMap<GeoVertex, ArrayList<GeoVertex>> subVertex = new HashMap<>();
        HashMap<GeoVertex, GeoVertex> auxVertex = new HashMap<>();
        HashMap<GeoVertex, HashMap<ArrayList<GeoVertex>, Double>> angleDict = new HashMap<>();
        if (t <= 1) return;
        System.out.println("---------------正在进行初始化---------------");
        long start = System.currentTimeMillis();
        while (id_vertex.containsKey(vertexId)) {
            GeoVertex vertex = id_vertex.get(vertexId++);
            ArrayList<GeoVertex> conVertices = new ArrayList<>(vertex.getConVertex());
            int count = conVertices.size();
            if (count > 2 || (count == 2 && vertexId > judgeId)) {
                ArrayList<GeoVertex> oldConVertices = new ArrayList<>(conVertices);
                oldConVertices.remove(0);
                newId = disconnectVertex(vertex, oldConVertices, newId, subVertex, auxVertex);
            }
        }
        System.out.println(id_edge);
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
        HashMap<GeoVertex, GeoVertex[]> reflectVertex = new HashMap<>();// 合并的两个节点 key为合并后的节点，value是删除的节点和它的相邻点
        start = System.currentTimeMillis();
        //总损失函数
        double totalCost = 999999999;
        while (t > 1) {
            System.out.println("当前温度：" + String.format("%.3f", t));
            long mid = System.currentTimeMillis();
            double angleCost = 0.0; //角度损失
            double disCost = 0.0; //孤立损失
            // 相连关系表
            HashMap<GeoVertex, ArrayList<GeoVertex>> conDict = new HashMap<>();
            HashMap<GeoVertex, GeoVertex[]> tempReflectVertex = new HashMap<>();
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
                    //指数归一化 exp(x)
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
                            if (chosenNum2 - 0.2 > Math.log(nowAngle + 1)) {
                                //这两点不合并，也就是它们的相邻点所组成的边不相连
                                disconnectVertex.add(val.getKey());
                                //计算孤立损失
                                disCost = disCost + nowAngle * 2;
                            } else {
                                //要合并
                                GeoVertex vertex1 = val.getKey().get(0);
                                GeoVertex vertex2 = val.getKey().get(1);
                                chosenVertex.add(vertex1);
                                chosenVertex.add(vertex2);
                                //计算角度损失
                                angleCost = angleCost + nowAngle * nowAngle;
                                //合并两个节点，相当于是在删除vertex2，但要保留vertex2的相邻关系
                                GeoVertex conVertex2 = conDict.get(vertex2).get(0);
                                conDict.get(vertex1).add(conVertex2);
                                conDict.get(conVertex2).remove(vertex2);
                                conDict.get(conVertex2).add(vertex1);
                                conDict.remove(vertex2);
                                tempReflectVertex.put(vertex1, new GeoVertex[]{vertex2, conVertex2});
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
            double nowCost = roadCost + disCost + angleCost;
            if (totalCost > nowCost) {//|| Math.random() < Math.exp((totalCost - nowCost) / t)
                totalCost = nowCost;
                resDict = new HashMap<>(conDict);
                reflectVertex = new HashMap<>(tempReflectVertex);
                System.out.println("接受当前目标函数值：" + totalCost);
            }
            t *= alpha;
            System.out.println("单次循环耗时：" + (System.currentTimeMillis() - mid) + "ms");
        }
        System.out.println("退火完成，共耗时：" + (System.currentTimeMillis() - start) + "ms");
        System.out.println("---------------正在重建相邻关系---------------");
        start = System.currentTimeMillis();
        // 根据最终相邻关系表重建相邻关系
        for (Map.Entry<GeoVertex, GeoVertex[]> entry : reflectVertex.entrySet()) {
            GeoEdge e = findEdgeVertices(entry.getValue()[0], entry.getValue()[1]);
            if (e == null) System.out.println("WTFFFFF!");
            vertices_edge.put(new ArrayList<>() {{
                add(entry.getKey());
                add(entry.getValue()[1]);
            }}, e);
        }
        for (GeoVertex[] v : reflectVertex.values()) removeVertex(v[0]);
        for (Map.Entry<GeoVertex, ArrayList<GeoVertex>> entry : resDict.entrySet()) {
            ArrayList<GeoVertex> conVertices = new ArrayList<>();
            for (GeoVertex v : entry.getValue()) {
                conVertices.add(id_vertex.get(v.getId()));
            }
            id_vertex.get(entry.getKey().getId()).setConVertex(conVertices);
        }
        roadDict = roadTrace();
        System.out.println("重建相邻关系完成，共耗时：" + (System.currentTimeMillis() - start) + "ms");
    }

    public HashMap<Integer, ArrayList<GeoVertex>> roadTrace() {
        HashMap<Integer, ArrayList<GeoVertex>> roadDict = new HashMap<>();
        int road_count = 1;
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
                roadDict.put(road_count, new ArrayList<>());
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
            roadDict.put(road_count, new ArrayList<>());
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
        return roadDict;
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
        return Math.sqrt(road_count);
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
        ds.createSchema(tb.buildFeatureType());
        ds.setCharset(StandardCharsets.UTF_8);
        //写入shapefile
        FeatureWriter<SimpleFeatureType, SimpleFeature> writer = ds.getFeatureWriter(ds.getTypeNames()[0], AUTO_COMMIT);
        for (Map.Entry<Integer, ArrayList<GeoVertex>> item : roadDict.entrySet()) {
            //获取当前写对象
            SimpleFeature feature = writer.next();
            //几何坐标，LineString是通过Coordinate[]构建的，同理MultiLineString是通过LineString[]构建
            Coordinate[] coordinates = new Coordinate[item.getValue().size()];
            //获取线段的属性集合，方便后面求众数、平均数等。
            HashMap<String, ArrayList<String>> att = new HashMap<>();
            int i = 0;
            GeoVertex preVertex = null;
            for (GeoVertex v : item.getValue()) {
                coordinates[i++] = new Coordinate(v.getCoord()[0], v.getCoord()[1]);
                if (preVertex != null) {
                    GeoEdge nowEdge = this.findEdgeVertices(preVertex, v);
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
        cos_theta = Double.parseDouble(String.format("%.10f", cos_theta));
        double theta = Math.acos(cos_theta);
        if (theta > Math.PI) theta = 2 * Math.PI - theta;

        return Math.PI - theta;
    }
}
