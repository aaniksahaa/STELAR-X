package stelarx.taxon;

import java.util.HashMap;
import java.util.Map;

/**
 * Maps taxon names <-> contiguous integer IDs [0, n).
 * Two-phase: register all names first, then lock.
 */
public class TaxonRegistry {

    private final Map<String, Integer> nameToId = new HashMap<>();
    private String[] idToName;
    private boolean locked = false;

    /** Register a name; returns its ID (new or existing). Thread-safe during pass 1. */
    public synchronized int register(String name) {
        if (locked) throw new IllegalStateException("Registry is locked");
        return nameToId.computeIfAbsent(name, k -> nameToId.size());
    }

    /** Lock registry and build the reverse array. Call after all registrations. */
    public void lock() {
        locked = true;
        idToName = new String[nameToId.size()];
        for (var e : nameToId.entrySet()) idToName[e.getValue()] = e.getKey();
    }

    public int getId(String name) {
        Integer id = nameToId.get(name);
        if (id == null) throw new IllegalArgumentException("Unknown taxon: " + name);
        return id;
    }

    /** Return a registered ID, or -1 when the name is outside this registry. */
    public int findId(String name) {
        Integer id = nameToId.get(name);
        return id == null ? -1 : id;
    }

    public String getName(int id) { return idToName[id]; }
    public int size()             { return nameToId.size(); }
    public boolean isLocked()     { return locked; }
}
