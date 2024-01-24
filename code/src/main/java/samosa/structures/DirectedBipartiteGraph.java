package samosa.structures;

import java.io.IOException;
import java.util.Arrays;
import samosa.utils.io.Writer;

/**
 *
 * @author giulia
 */
public class DirectedBipartiteGraph {

    // out-neighbors of the left nodes
    private final int[][] left_nodes_out_adj;
    // in-neighbors of the left nodes
    private final int[][] left_nodes_in_adj;
    // out-neighbors of the right nodes
    private final int[][] right_nodes_out_adj;
    // in-neighbors of the right nodes
    private final int[][] right_nodes_in_adj;
    // number of edges from L to R nodes
    private final int numLeft2RightEdges;
    // number of edges from R to L nodes
    private final int numRight2LeftEdges;

    public DirectedBipartiteGraph(
            int[][] left_P,
            int[][] left_M,
            int[][] right_P,
            int[][] right_M,
            int size_D0,
            int size_D1) {

        this.left_nodes_out_adj = left_P;
        this.left_nodes_in_adj = left_M;
        this.right_nodes_out_adj = right_P;
        this.right_nodes_in_adj = right_M;
        this.numLeft2RightEdges = size_D0;
        this.numRight2LeftEdges = size_D1;
    }

    /**
     * 
     * @return out-neighbors of the left nodes
     */
    public int[][] getLeftP() {
        return left_nodes_out_adj;
    }

    /**
     * 
     * @return in-neighbors of the left nodes
     */
    public int[][] getLeftM() {
        return left_nodes_in_adj;
    }

    /**
     * 
     * @return out-neighbors of the right nodes
     */
    public int[][] getRightP() {
        return right_nodes_out_adj;
    }

    /**
     * 
     * @return in-neighbors of the right nodes
     */
    public int[][] getRightM() {
        return right_nodes_in_adj;
    }

    /**
     * 
     * @return number of edges from L to R nodes
     */
    public int numLeft2RightEdges() {
        return numLeft2RightEdges;
    }

    /**
     * 
     * @return number of edges from R to L nodes
     */
    public int numRight2LeftEdges() {
        return numRight2LeftEdges;
    }

    /**
     * 
     * @return total number of edges
     */
    public int getNumEdges() {
        return numLeft2RightEdges + numRight2LeftEdges;
    }

    /**
     * Returns the adjacency matrix that will be converted into a transactional
     * database for frequent itemset mining.
     *
     * @param pos kind of dataset of itemsets
     * @return adjacency matrix associated with pos
     */
    public int[][] getTransactions(int pos) {

        if (pos == 0) {
            return getRightM();
        }
        if (pos == 1) {
            return getRightP();
        }
        if (pos == 2) {
            return getLeftM();
        }
        return getLeftP();
    }
    
    /**
     * Write the dataset of itemsets corresponding to pos to disk.
     * @param pos kind of dataset of itemsets
     * @param fileName output file name
     * @throws IOException 
     */
    public void saveTransactions(int pos, String fileName) throws IOException {
        int[][] trans = getTransactions(pos);
        Writer.writeTransactions(trans, fileName);
    }

    /**
     * 
     * @return a deep copy of this DirBipGraph
     */
    public DirectedBipartiteGraph copy() {
        int[][] newLP = new int[left_nodes_out_adj.length][];
        for (int i = 0; i < newLP.length; i++) {
            newLP[i] = Arrays.copyOf(left_nodes_out_adj[i], left_nodes_out_adj[i].length);
        }
        int[][] newLM = new int[left_nodes_in_adj.length][];
        for (int i = 0; i < newLM.length; i++) {
            newLM[i] = Arrays.copyOf(left_nodes_in_adj[i], left_nodes_in_adj[i].length);
        }
        int[][] newRP = new int[right_nodes_out_adj.length][];
        for (int i = 0; i < newRP.length; i++) {
            newRP[i] = Arrays.copyOf(right_nodes_out_adj[i], right_nodes_out_adj[i].length);
        }
        int[][] newRM = new int[right_nodes_in_adj.length][];
        for (int i = 0; i < newRM.length; i++) {
            newRM[i] = Arrays.copyOf(right_nodes_in_adj[i], right_nodes_in_adj[i].length);
        }
        return new DirectedBipartiteGraph(newLP, newLM, newRP, newRM, numLeft2RightEdges, numRight2LeftEdges);
    }

    /**
     * 
     * @param other a DirBipGraph
     * @return true if this is identical to other (no isomorphic but identical);
     * false otherwise
     */
    public boolean equals(DirectedBipartiteGraph other) {
        for (int i = 0; i < left_nodes_out_adj.length; i++) {
            if (!Arrays.equals(left_nodes_out_adj[i], other.left_nodes_out_adj[i])) {
                return false;
            }
        }
        for (int i = 0; i < left_nodes_in_adj.length; i++) {
            if (!Arrays.equals(left_nodes_in_adj[i], other.left_nodes_in_adj[i])) {
                return false;
            }
        }
        for (int i = 0; i < right_nodes_out_adj.length; i++) {
            if (!Arrays.equals(right_nodes_out_adj[i], other.right_nodes_out_adj[i])) {
                return false;
            }
        }
        for (int i = 0; i < right_nodes_in_adj.length; i++) {
            if (!Arrays.equals(right_nodes_in_adj[i], other.right_nodes_in_adj[i])) {
                return false;
            }
        }
        return true;
    }
}
