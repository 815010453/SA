package lyd.SA;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class GeoVertex {

    private final int id; // 唯一标识符 id
    private final double[] coord; // 点坐标[x, y] 单位为m
    private ArrayList<GeoVertex> conVertex; // 相邻点[vertex1, vertex2, ...]
    private final HashMap<String, String> nodeAttributes; // 节点属性{'x':12315,'y':21546,...}

    private final ArrayList<GeoEdge> edges;
    private final ArrayList<GeoSegment> segments;

    public void setConVertex(ArrayList<GeoVertex> conVertex) {
        this.conVertex = conVertex;
    }

    public GeoVertex(int v_id, HashMap<String, String> att, double[] coord) {
        this.conVertex = new ArrayList<>();
        this.id = v_id;
        this.coord = coord;
        this.nodeAttributes = att;
        this.edges = new ArrayList<>();
        this.segments = new ArrayList<>();
    }

    public ArrayList<GeoVertex> getConVertex() {
        return conVertex;
    }

    public int getId() {
        return id;
    }

    public double[] getCoord() {
        return coord;
    }

    public void addEdge(GeoEdge edge) {
        if(!this.edges.contains(edge)){
            this.edges.add(edge);
        }
    }

    public void removeEdge(GeoEdge edge) {
        this.edges.remove(edge);
    }

    public void addSegment(GeoSegment segment) {
        if(!this.segments.contains(segment)) this.segments.add(segment);
    }

    public ArrayList<GeoSegment> getSegments() {
        return segments;
    }

    public ArrayList<GeoEdge> getEdges() {
        return edges;
    }

    public void removeSegment(GeoSegment seg) {
        segments.remove(seg);
    }

    public HashMap<String, String> getNodeAttributes() {
        return nodeAttributes;
    }

    public void add_conVertex(GeoVertex vertex) {
        if (vertex != this && !conVertex.contains(vertex)) {
            this.conVertex.add(vertex);
            vertex.add_conVertex(this);
        }
    }

    public void remove_conVertex(GeoVertex vertex) {
        if (this.conVertex.contains(vertex)) {
            this.conVertex.remove(vertex);
            vertex.remove_conVertex(this);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        GeoVertex geoVertex = (GeoVertex) o;
        return id == geoVertex.id;
    }
    public void clearSegment(){
        this.segments.clear();
    }
    @Override
    public int hashCode() {
        return this.id;
    }

    @Override
    public String toString() {
        return "" + id;
    }
}
