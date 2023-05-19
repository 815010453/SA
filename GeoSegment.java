package lyd.SA;

import java.util.ArrayList;

public class GeoSegment {
    private final GeoVertex start;
    private final GeoVertex vertex2;
    private final ArrayList<GeoEdge> edges;
    private final int id;

    public GeoSegment(int id, ArrayList<GeoEdge> edges, GeoVertex start, GeoVertex vertex2) {
        this.edges = edges;
        this.id = id;
        this.start = start;
        this.vertex2 = vertex2;
    }

    public GeoVertex getStart() {
        return start;
    }

    public GeoVertex getVertex2() {
        return vertex2;
    }

    public ArrayList<GeoEdge> getEdges() {
        return edges;
    }

    public int getId() {
        return id;
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
