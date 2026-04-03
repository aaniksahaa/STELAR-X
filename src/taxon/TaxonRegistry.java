package taxon;

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
        Integer existing = nameToId.get(name);
        if (existing != null) return existing;
        int id = nameToId.size();
        nameToId.put(name, id);
        return id;
    }

    /** Lock registry and build the reverse array. Call after all registrations. */
    public void lock() {
        locked = true;
        idToName = new String[nameToId.size()];
        for (Map.Entry<String, Integer> e : nameToId.entrySet()) {
            idToName[e.getValue()] = e.getKey();
        }
    }

    public int getId(String name) {
        Integer id = nameToId.get(name);
        if (id == null) throw new IllegalArgumentException("Unknown taxon: " + name);
        return id;
    }

    public String getName(int id) { return idToName[id]; }
    public int size()             { return nameToId.size(); }
    public boolean isLocked()     { return locked; }
}
