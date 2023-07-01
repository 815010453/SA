package lyd.SA;

import java.util.ArrayList;

public class GeoSegment {
    private final ArrayList<GeoEdge> edges;
    private final ArrayList<GeoVertex> vertices;
    private final int id;

    public GeoSegment(int id, ArrayList<GeoEdge> edges,ArrayList<GeoVertex> vertices) {
        this.edges = edges;
        this.id = id;
        this.vertices = vertices;
    }

    public ArrayList<GeoEdge> getEdges() {
        return edges;
    }

    public ArrayList<GeoVertex> getVertices(){
        return vertices;
    }

    public int getId() {
        return id;
    }
    public int containVertex(GeoVertex vertex){
        return vertices.indexOf(vertex);
    }

    public int containEdge(GeoEdge edge){
        return edges.indexOf(edge);
    }

    @Override
    public String toString() {
        return "" + id;
    }

    @Override
    public int hashCode() {
        return this.id;
    }
}
