package lyd.SA;


import java.util.*;
import java.util.logging.Logger;



/**
 * Prompts the user for a shapefile and displays the contents on the screen in a map frame.
 *
 * <p>This is the GeoTools Quickstart application used in documentation and tutorials. *
 */
public class Main {

    private static final Logger LOGGER = org.geotools.util.logging.Logging.getLogger(Main.class);

    /**
     * GeoTools Quickstart demo application. Prompts the user for a shapefile and displays its
     * contents on the screen in a map frame
     */
    public static void main(String[] args) throws Exception {
        //File file = JFileDataStoreChooser.showOpenFile(".shp", null);

        String path = "C:\\Users\\81501\\Desktop\\SA\\src\\main\\resources\\data\\cd_road\\cd_road_processed.shp";
        long start = System.currentTimeMillis();
        GeoGragh gragh = new GeoGragh("cd_road", path);
        System.out.println("数据构图完成，共耗时" + (System.currentTimeMillis() - start) + "ms");
        start = System.currentTimeMillis();
        gragh.reconstructEdgeSA(100, 0.8);
        System.out.println("模拟退火算法重构完成，共耗时" + (System.currentTimeMillis() - start) + "ms");
        HashMap<Integer, ArrayList<GeoVertex>> myDict = gragh.getRoadDict();
        System.out.println(""+myDict.entrySet().size());
        start = System.currentTimeMillis();
        gragh.outputGeographRoads("C:\\Users\\81501\\Desktop\\SA\\src\\main\\resources\\out\\cd_road_sa.shp");
        System.out.println("数据写入完成，共耗时" + (System.currentTimeMillis() - start) + "ms");

        /*
        GeoGragh gragh = new GeoGragh("cd_road");
        ArrayList<GeoVertex> vertices = new ArrayList<>();

        GeoVertex vertex;
        vertex = new GeoVertex(1, new HashMap<>(), new double[]{0.0, 1.0});
        vertices.add(vertex);
        vertex = new GeoVertex(2, new HashMap<>(), new double[]{2.0, 1.0});
        vertices.add(vertex);
        vertex = new GeoVertex(3, new HashMap<>(), new double[]{1.0, 2.0});
        vertices.add(vertex);
        vertex = new GeoVertex(4, new HashMap<>(), new double[]{0.0, 0.0});
        vertices.add(vertex);
        vertex = new GeoVertex(5, new HashMap<>(), new double[]{2.0, 0.0});
        vertices.add(vertex);
        vertex = new GeoVertex(6, new HashMap<>(), new double[]{0.0, 2.0});
        vertices.add(vertex);
        vertex = new GeoVertex(7, new HashMap<>(), new double[]{2.0, 2.0});
        vertices.add(vertex);

        GeoEdge edge;
        edge = new GeoEdge(1, vertices.get(0), vertices.get(1), new HashMap<>());
        gragh.addEdge(edge);
        edge = new GeoEdge(2, vertices.get(0), vertices.get(2), new HashMap<>());
        gragh.addEdge(edge);
        edge = new GeoEdge(3, vertices.get(0), vertices.get(3), new HashMap<>());
        gragh.addEdge(edge);
        edge = new GeoEdge(4, vertices.get(1), vertices.get(2), new HashMap<>());
        gragh.addEdge(edge);
        edge = new GeoEdge(5, vertices.get(1), vertices.get(4), new HashMap<>());
        gragh.addEdge(edge);
        edge = new GeoEdge(6, vertices.get(3), vertices.get(4), new HashMap<>());
        gragh.addEdge(edge);
        edge = new GeoEdge(7, vertices.get(0), vertices.get(5), new HashMap<>());
        gragh.addEdge(edge);
        edge = new GeoEdge(8, vertices.get(1), vertices.get(6), new HashMap<>());
        gragh.addEdge(edge);

        gragh.reconstructEdgeSA(100, 0.99);
        HashMap<Integer, ArrayList<GeoVertex>> myDict = gragh.getRoadDict();
        for (Map.Entry<Integer, ArrayList<GeoVertex>> entry : myDict.entrySet()) {
            for (GeoVertex v : entry.getValue()) {
                System.out.print(Arrays.toString(v.getCoord()));
            }
            System.out.println();
        }
        System.out.println("" + myDict);
        /*
        // Create a map content and add our shapefile to it
        MapContent map = new MapContent();
        map.setTitle("Using cached features");
        Style style = SLD.createSimpleStyle(featureSource.getSchema());
        Layer layer = new FeatureLayer(cachedSource, style);
        map.addLayer(layer);

        // Now display the map
        JMapFrame.showMap(map);
        */
    }
}