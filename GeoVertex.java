package lyd.SA;

import java.util.ArrayList;
import java.util.HashMap;

public class GeoVertex {

    private int id; // 唯一标识符 id
    private double[] coord; // 点坐标[x, y] 单位为m
    private ArrayList<GeoVertex> conVertex; // 相邻点[vertex1, vertex2, ...]
    private HashMap<String, String> nodeAttributes; // 节点属性{'x':12315,'y':21546,...}

    public void setConVertex(ArrayList<GeoVertex> conVertex) {
        this.conVertex = conVertex;
    }

    public void setId(int id) {
        this.id = id;
    }

    public void setCoord(double[] coord) {
        this.coord = coord;
    }

    public void setNodeAttributes(HashMap<String, String> nodeAttributes) {
        this.nodeAttributes = nodeAttributes;
    }

    public GeoVertex(int v_id, HashMap<String, String> att, double[] coord) {
        this.conVertex = new ArrayList<>();
        this.id = v_id;
        this.coord = coord;
        this.nodeAttributes = att;
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

    @Override
    public int hashCode() {
        return this.id;
    }

    @Override
    public String toString() {
        return "" + id;
    }
}
