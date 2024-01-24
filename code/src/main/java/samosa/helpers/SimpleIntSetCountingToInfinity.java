package samosa.helpers;

import java.util.Arrays;

public class SimpleIntSetCountingToInfinity {
    // id of the hyperedge head/tail we are currently filling
    // (we start from 1 beacause map_key_presence is initialized with zeros)
    // this value does not correspond to the id of the hyperedge in the hypergraph
    protected int threshold;
    // max number of hyperedges we can create
    protected final int maximum_allowed_threshold;
    // for each node, id of the last hyperedge head/tail where we inserted the node
    protected int[] map_key_presence;
    // size of this set
    protected int size;

    public SimpleIntSetCountingToInfinity(int max_key_value) {
        this.threshold = 1;
        this.maximum_allowed_threshold = Integer.MAX_VALUE - 2;
        this.map_key_presence = new int[max_key_value + 1];
        this.size = 0;
    }

    public SimpleIntSetCountingToInfinity(final SimpleIntSetCountingToInfinity set) {
        this.threshold = set.threshold;
        this.maximum_allowed_threshold = set.maximum_allowed_threshold;
        this.map_key_presence = Arrays.copyOf(set.map_key_presence, set.map_key_presence.length);
        this.size = set.size;
    }

    /**
     * 
     * @return number of elements in the head/tail we are filling
     */
    public int size() {
        return this.size;
    }

    /**
     * Adds a new node to the head/tail with id = threshold, if the node is not
     * already present.
     * @param key node to add
     * @return true if key has been added to the head/tail; false otherwise
     */
    public boolean informativeAdd(int key) {
        if (!this.contains(key)) {
            // update the id associated to the node key
            this.map_key_presence[key] = this.threshold;
            // the size of the current head/tail increases by 1
            this.size++;
            return true;
        }
        return false;
    }

    /**
     * 
     * @param key element to search
     * @return true if the element is stored in the head/tail with id = threshold; 
     * false otherwise
     */
    public boolean contains(int key) {
        return (this.map_key_presence[key] >= this.threshold);
    }

    /**
     * Prepares the data structure for the population of the next head/tail.
     */
    public void clear() {
        // we increase the counter by 1 because we are ready to start
        // populating a new head/tail
        this.threshold++;
        // restart from 1 if we reached the maximum id
        if (this.threshold > this.maximum_allowed_threshold) {
            //
            this.threshold = 1;
            // restart the array storing the id of the last head/tail including
            // each node
            Arrays.fill(this.map_key_presence, 0);
        }
        this.size = 0;
    }

}
