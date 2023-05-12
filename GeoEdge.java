package lyd.SA;

import java.util.HashMap;

public class GeoEdge {
    private int id; // 唯一标识符 id
    private GeoVertex[] geoVertices; // 组成该边的点
    private HashMap<String, String> edgeAttribute; // 该边的属性 {'fclass': 'highway', 'name': '公路', ...}

    public GeoEdge(int id, GeoVertex vertex_a, GeoVertex vertex_b, HashMap<String, String> att) {
        this.id = id;
        this.geoVertices = new GeoVertex[]{vertex_a, vertex_b};
        this.edgeAttribute = att;
    }

    public int getId() {
        return id;
    }

    public GeoVertex[] getGeoVertices() {
        return geoVertices;
    }


    public HashMap<String, String> getEdgeAttribute() {
        return edgeAttribute;
    }

    public void setId(int id) {
        this.id = id;
    }

    public void setGeoVertices(GeoVertex[] geoVertices) {
        this.geoVertices = geoVertices;
    }

    public void setEdgeAttribute(HashMap<String, String> edgeAttribute) {
        this.edgeAttribute = edgeAttribute;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        GeoEdge geoEdge = (GeoEdge) o;
        return id == geoEdge.id;
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
