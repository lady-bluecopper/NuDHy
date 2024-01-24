package samosa.utils.io;

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

/**
 *
 * @author giulia
 */
public class Transformer {
    
    Object2IntOpenHashMap<String> outer2innerNodeIds;
    Int2ObjectOpenHashMap<String> inner2outerNodeIds;
    
    public Transformer() {
        this.outer2innerNodeIds = new Object2IntOpenHashMap();
        this.outer2innerNodeIds.defaultReturnValue(-1);
        this.inner2outerNodeIds = new Int2ObjectOpenHashMap();
        this.inner2outerNodeIds.defaultReturnValue("");
    }
    
    /**
     * Associates an inner id to outer.
     * 
     * @param outer node id
     */
    public void saveMapping(String outer) {
        this.outer2innerNodeIds.putIfAbsent(outer, this.outer2innerNodeIds.size());
        this.inner2outerNodeIds.putIfAbsent(this.outer2innerNodeIds.getInt(outer), outer);
    }
    
    /**
     * 
     * @param outer node id
     * @return inner id associated with outer
     */
    public int getInnerId(String outer) {
        if (!this.outer2innerNodeIds.containsKey(outer)) {
            System.out.println(outer + " not in map of size " + this.outer2innerNodeIds.size());
        }
        return this.outer2innerNodeIds.getInt(outer);
    }
    
    /**
     * 
     * @param inner node id
     * @return outer id associated with inner
     */
    public String getOuterId(int inner) {
        if (!this.inner2outerNodeIds.containsKey(inner)) {
            System.out.println(inner + " not in map of size " + this.inner2outerNodeIds.size());
        }
        return this.inner2outerNodeIds.get(inner);
    }
    
    public void printMapping() {
        System.out.println("Inner 2 Outer");
        System.out.println(this.inner2outerNodeIds.toString());
        System.out.println("Outer 2 Inner");
        System.out.println(this.outer2innerNodeIds.toString());
    }
}
