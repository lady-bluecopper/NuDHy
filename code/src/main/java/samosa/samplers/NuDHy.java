package samosa.samplers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import samosa.helpers.ArrayHandler;
import samosa.structures.InnerOuterIDTransformer;
import samosa.utils.Timer;
import samosa.utils.Utils;

public abstract class NuDHy extends InnerOuterIDTransformer implements Sampler {

    // out-neighbors left nodes
    protected int[][] left_P = null;
    // in-neighbors left nodes
    protected int[][] left_M = null;
    // out-neighbors right nodes
    protected int[][] right_P = null;
    // in-neighbors right nodes
    protected int[][] right_M = null;
    // number of edges from left to right nodes
    public int size_D0 = -1;
    // number of edges from right to left nodes
    public int size_D1 = -1;

    public NuDHy() {
        super();
    }

    public NuDHy(final NuDHy obj) {
        super(obj);
        this.left_P = new int[obj.left_P.length][];
        Utils.initializeArrayOfArrays(obj.left_P, left_P);
        this.left_M = new int[obj.left_M.length][];
        Utils.initializeArrayOfArrays(obj.left_M, left_M);
        this.right_P = new int[obj.right_P.length][];
        Utils.initializeArrayOfArrays(obj.right_P, right_P);
        this.right_M = new int[obj.right_M.length][];
        Utils.initializeArrayOfArrays(obj.right_M, right_M);
        this.size_D0 = obj.size_D0;
        this.size_D1 = obj.size_D1;
    }

    /**
     * Determines whether S_with_outer_ids appears in the head of some hyperedge of the
     * hypergraph. The array S_with_outer_ids contains the outer ids of the nodes.
     * @param S_with_outer_ids
     * @return true if S_with_outer_ids appears in the head of some 
     * hyperedge; false otherwise
     */
    protected boolean isSinAtLeastOneHyperedgeHead(final String[] S_with_outer_ids) {
        return this.isSinAtLeastOneHyperedge(S_with_outer_ids, true);
    }

    /**
     * Determines whether S_with_outer_ids appears in the tail of some hyperedge of the
     * hypergraph. The array S_with_outer_ids contains the outer ids of the nodes.
     * @param S_with_outer_ids
     * @return true if S_with_outer_ids appears in the tail of some 
     * hyperedge; false otherwise
     */
    protected boolean isSinAtLeastOneHyperedgeTail(final String[] S_with_outer_ids) {
        return this.isSinAtLeastOneHyperedge(S_with_outer_ids, false);
    }

    /**
     * Determines whether S_with_outer_ids appears in some hyperedge of the
     * hypergraph. 
     * The array S_with_outer_ids contains the outer ids of the nodes.
     * The node group is searched either in the heads or in the tails of the
     * hyperedges.
     * @param S_with_outer_ids node group to search (outer ids)
     * @param head_or_tail whether the search must be performed in the heads
     * or in the tails of the hyperedges
     * @return true if S_with_outer_ids appears in the head/tail of some 
     * hyperedge; false otherwise
     */
    public boolean isSinAtLeastOneHyperedge(final String[] S_with_outer_ids, 
            final boolean head_or_tail) {
        final int[] S = new int[S_with_outer_ids.length];
        for (int i = 0; i < S_with_outer_ids.length; i++) {
            if (!this.map__outer_inner_id.containsKey(S_with_outer_ids[i])) {
                return false;
            }
            S[i] = this.map__outer_inner_id.getInt(S_with_outer_ids[i]);
        }
        Arrays.sort(S);
        boolean S_is_contained_in_the_current_hyperedge = false;
        
        int[][] left;
        int[][] right;
        if (head_or_tail) {
            left = this.left_P;
            right = this.right_M;
        } else {
            left = this.left_M;
            right = this.right_P;
        }
        int i;
        int best_s = -1;
        // select the vertex belonging to the minimum number of hyperedges.
        int best_s_value = this.right_M.length;
        int s;
        int s_value;
        for (i = 0; i < S.length; i++) {
            s = S[i];
            s_value = left[s].length;
            if (s_value < best_s_value) {
                best_s_value = s_value;
                best_s = s;
            }
        }
        // iterate over all the hyperedges of interest.
        ArrayHandler ah = new ArrayHandler();
        int h;
        for (i = 0; i < left[best_s].length; i++) {
            h = left[best_s][i];
            S_is_contained_in_the_current_hyperedge = ah.isAcontainedInB(S, right[h]);
            if (S_is_contained_in_the_current_hyperedge) {
                return true;
            }
        }
        return false;
    }

    /**
     * Read hypergraph and initialize relevant structures.
     * @param input_file_name path to input hypergrph
     * @throws Exception 
     */
    public void parseInputFile(final String input_file_name) throws Exception {
        //
        int[] max_id_l__max_id_r = super.parseInputFileForCreatingOuterInnerBijection(input_file_name);
        int max_id_l = max_id_l__max_id_r[0];
        int max_id_r = max_id_l__max_id_r[1];
        //
        IntArrayList[] tmp_left_out = new IntArrayList[max_id_l + 1];
        for (int i = 0; i < tmp_left_out.length; i++) {
            tmp_left_out[i] = new IntArrayList();
        }
        IntArrayList[] tmp_left_in = new IntArrayList[max_id_l + 1];
        for (int i = 0; i < tmp_left_in.length; i++) {
            tmp_left_in[i] = new IntArrayList();
        }
        IntArrayList[] tmp_right_out = new IntArrayList[max_id_r + 1];
        for (int i = 0; i < tmp_right_out.length; i++) {
            tmp_right_out[i] = new IntArrayList();
        }
        IntArrayList[] tmp_right_in = new IntArrayList[max_id_r + 1];
        for (int i = 0; i < tmp_right_in.length; i++) {
            tmp_right_in[i] = new IntArrayList();
        }
        BufferedReader br = new BufferedReader(new FileReader(input_file_name));
        String line;
        int c_r_id = -1;
        int c_l_id = -1;
        while ((line = br.readLine()) != null) {
            c_r_id++;
            String[] s_0__s_1 = line.split("\t");
            String[] s_0 = s_0__s_1[0].split(",");
            String[] s_1 = null;
            if (s_0__s_1.length == 2) {
                s_1 = s_0__s_1[1].split(",");
            }
            for (String c_l_outher_id : s_0) {
                c_l_id = this.map__outer_inner_id.getInt(c_l_outher_id);
                tmp_left_out[c_l_id].add(c_r_id);
                tmp_right_in[c_r_id].add(c_l_id);
            }
            if (s_1 != null) {
                for (String c_l_outher_id : s_1) {
                    c_l_id = this.map__outer_inner_id.getInt(c_l_outher_id);
                    tmp_left_in[c_l_id].add(c_r_id);
                    tmp_right_out[c_r_id].add(c_l_id);
                }
            }
        }
        br.close();
        //
        this.left_P = new int[max_id_l + 1][];
        this.size_D0 = Utils.convertTempLeftRightInDefinitiveLeftRight(tmp_left_out, left_P);
        tmp_left_out = null;
        this.left_M = new int[max_id_l + 1][];
        Utils.convertTempLeftRightInDefinitiveLeftRight(tmp_left_in, left_M);
        tmp_left_in = null;
        this.right_P = new int[max_id_r + 1][];
        this.size_D1 = Utils.convertTempLeftRightInDefinitiveLeftRight(tmp_right_out, right_P);
        tmp_right_out = null;
        this.right_M = new int[max_id_r + 1][];
        Utils.convertTempLeftRightInDefinitiveLeftRight(tmp_right_in, right_M);
        tmp_right_in = null;
        //
        System.gc();
    }

    /**
     * Read hypergraph and initialize relevant structures.
     * @param br BufferedReader to read the file and get the node ids
     * @param br2 BufferedReader to read the adjacency matrix
     * @throws Exception 
     */
    public void parseInputFile(BufferedReader br, BufferedReader br2) throws Exception {

        ObjectOpenHashSet<String> set__outher_ids = new ObjectOpenHashSet();
        String line;
        int max_id_r = -1;
        String[][] array_s_0__s_1 = new String[2][];
        // read the first time
        while ((line = br.readLine()) != null) {
            max_id_r++;
            String[] s_0__s_1 = line.split("\t");
            String[] s_0 = s_0__s_1[0].split(",");
            String[] s_1 = s_0__s_1[1].split(",");
            array_s_0__s_1[0] = s_0;
            array_s_0__s_1[1] = s_1;
            for (String[] s_x : array_s_0__s_1) {
                for (String id_as_string : s_x) {
                    set__outher_ids.add(id_as_string);
                }
            }
        }
        br.close();

        int max_id_l = set__outher_ids.size() - 1;
        int c_id_l = 0;
        String[] sorted__outher_ids = new String[set__outher_ids.size()];
        set__outher_ids.toArray(sorted__outher_ids);
        Arrays.sort(sorted__outher_ids);
        set__outher_ids = null;
        this.map__inner_outer_id = sorted__outher_ids;
        this.map__outer_inner_id = new Object2IntOpenHashMap(3 * this.map__inner_outer_id.length);
        this.map__outer_inner_id.defaultReturnValue(-1);
        for (c_id_l = 0; c_id_l < this.map__inner_outer_id.length; c_id_l++) {
            this.map__outer_inner_id.put(this.map__inner_outer_id[c_id_l], c_id_l);
        }
        // we first put the in- and out-neighbors into ArrayList and
        // then we transform the lists into arrays
        IntArrayList[] tmp_left_out = new IntArrayList[max_id_l + 1];
        for (int i = 0; i < tmp_left_out.length; i++) {
            tmp_left_out[i] = new IntArrayList();
        }
        IntArrayList[] tmp_left_in = new IntArrayList[max_id_l + 1];
        for (int i = 0; i < tmp_left_in.length; i++) {
            tmp_left_in[i] = new IntArrayList();
        }
        IntArrayList[] tmp_right_out = new IntArrayList[max_id_r + 1];
        for (int i = 0; i < tmp_right_out.length; i++) {
            tmp_right_out[i] = new IntArrayList();
        }
        IntArrayList[] tmp_right_in = new IntArrayList[max_id_r + 1];
        for (int i = 0; i < tmp_right_in.length; i++) {
            tmp_right_in[i] = new IntArrayList();
        }
        int c_r_id = -1;
        int c_l_id = -1;
        // read the second time
        while ((line = br2.readLine()) != null) {
            c_r_id++;
            String[] s_0__s_1 = line.split("\t");
            String[] s_0 = s_0__s_1[0].split(",");
            String[] s_1 = s_0__s_1[1].split(",");
            for (String c_l_outher_id : s_0) {
                c_l_id = this.map__outer_inner_id.getInt(c_l_outher_id);
                tmp_left_out[c_l_id].add(c_r_id);
                tmp_right_in[c_r_id].add(c_l_id);
            }
            for (String c_l_outher_id : s_1) {
                c_l_id = this.map__outer_inner_id.getInt(c_l_outher_id);
                tmp_left_in[c_l_id].add(c_r_id);
                tmp_right_out[c_r_id].add(c_l_id);
            }
        }
        br2.close();

        this.left_P = new int[max_id_l + 1][];
        this.size_D0 = Utils.convertTempLeftRightInDefinitiveLeftRight(tmp_left_out, left_P);
        tmp_left_out = null;
        this.left_M = new int[max_id_l + 1][];
        Utils.convertTempLeftRightInDefinitiveLeftRight(tmp_left_in, left_M);
        tmp_left_in = null;
        this.right_P = new int[max_id_r + 1][];
        this.size_D1 = Utils.convertTempLeftRightInDefinitiveLeftRight(tmp_right_out, right_P);
        tmp_right_out = null;
        this.right_M = new int[max_id_r + 1][];
        Utils.convertTempLeftRightInDefinitiveLeftRight(tmp_right_in, right_M);
        tmp_right_in = null;

        System.gc();
    }

    /**
     * Save hypergraph to disk.
     * @param input_file_complete_name path to original hypergraph
     * @param max_number_of_iterations number of steps performed in the Markov graph
     * @param random_seed seed used to generate the sample
     * @param execTimeINmsec time spent to generate the sample
     * @param number_of_rejections number of self-loops performed in the Markov graph
     * @param output_directory directory to store the random hypergraph
     * @throws Exception 
     */
    public void printOnFile(final String input_file_complete_name, 
            final int max_number_of_iterations,
            final int random_seed, 
            final long execTimeINmsec,
            final int number_of_rejections, 
            final String output_directory) throws Exception {

        final int index = this.getClass().toString().lastIndexOf(".");
        // create output file name including all the relevant hyper-parameters
        String output_file_name = input_file_complete_name.replace(".tsv", "")
                .substring(input_file_complete_name.lastIndexOf('/') + 1);
        String final_out_file_name = output_directory + "/" 
                + output_file_name + "__algorithm_"
                + this.getClass().toString().substring(index + 1)
                + "__iterations_" + max_number_of_iterations
                + "__randomSeed_" + random_seed
                + "__excTimemsec_" + execTimeINmsec 
                + "__numRejs_" + number_of_rejections + "__.tsv";
        System.out.println(final_out_file_name);

        
        StringBuilder sb_head = new StringBuilder();
        StringBuilder sb_tail = new StringBuilder();
        BufferedWriter bw = new BufferedWriter(new FileWriter(final_out_file_name));
        String out_head, out_tail;
        for (int c_node_id = 0; c_node_id < this.right_M.length; c_node_id++) {
            // string representing the head of the hyperedge
            out_head = toStringHeadOrTailOfHyperedge(this.right_M[c_node_id], sb_head);
            // string representing the tail of the hyperedge
            out_tail = toStringHeadOrTailOfHyperedge(this.right_P[c_node_id], sb_tail);
            bw.write(out_head);
            bw.write("\t");
            bw.write(out_tail);
            bw.write("\n");
        }
        bw.close();
    }

    /**
     * Fix the input bias if there are no edges with a certain direction.
     * @param bias input bias
     * @param size_D0 number of left to right edges
     * @param size_D1 number of right to left edges
     * @return 0 if size_D0 is 0; 1 if size_D1 is 0; bias otherwise
     */
    protected double rectifyBias(double bias, final int size_D0, final int size_D1) {

        double new_bias = bias;

        if ((size_D0 != 0) && (size_D1 != 0)) {
            if (bias <= Double.MIN_VALUE || bias >= 1. - Double.MIN_VALUE) {
                return -1;
            }
        }
        if ((size_D0 != 0) && (size_D1 == 0)) {
            // no right to left edges; we can swap only left to right edges 
            new_bias = 1.;
        }
        if ((size_D0 == 0) && (size_D1 != 0)) {
            // no left to right edges; we can swap only right to left edges 
            new_bias = 0.;
        }
        return new_bias;
    }

    /**
     * Starting from the original hypergraph, move in the Markov graph and store
     * on disk the current state after numSwaps steps.
     * @param fileName path to original hypergraph
     * @param numSwaps number of steps in the Markov graph before returning the sample
     * @param seed seed for reproducibility
     * @param timer object to store the running times
     * @param save whether the sample should be saved on disk
     * @throws java.lang.Exception
     */
    public abstract void sample(String fileName, int numSwaps, long seed, Timer timer, boolean save) throws Exception;
    
    /**
     * Moves in the Markov graph for max_number_of_iterations steps.
     *
     * @param max_number_of_iterations number of iterations
     * @param seed seed for reproducibility
     * @param timer Timer to store running times
     * @return number of self-loops performed
     */
    public abstract int runMC(final int max_number_of_iterations, int seed, Timer timer);

    /**
     * Starting from the original hypergraph, move in the Markov graph and store
     * on disk the current state after k steps up to max_number_of_iterations steps.
     * @param max_number_of_iterations maximum number of steps in the Markov graph
     * @param seed seed for reproducibility
     * @param file_name path to original hypergraph
     * @return number of self-loops performed
     * @throws Exception 
     */
    public abstract int runMC_multiFlush(final int max_number_of_iterations, int seed, String file_name)
            throws Exception;

    /**
     *
     * @return total number of edges in the bipartite graph
     */
    public int getTotalEdges() {
        return size_D0 + size_D1;
    }

    /**
     * Initialize a new NuDHy sampler of type name.
     * @param name sampler type
     * @return an instance of sampler of type name
     * @throws ClassNotFoundException
     */
    public static NuDHy getSampler(String name) throws ClassNotFoundException {
        if (name.equalsIgnoreCase(NuDHy_Degs.class.getName())
                || name.equalsIgnoreCase("NuDHy_A")
                || name.equalsIgnoreCase("NuDHy_Degs")) {
            return new NuDHy_Degs();
        }
        if (name.equalsIgnoreCase(NuDHy_BIOT.class.getName())
                || name.equalsIgnoreCase("NuDHy_C")
                || name.equalsIgnoreCase("NuDHy_BIOT")) {
            return new NuDHy_BIOT();
        }
        throw new ClassNotFoundException(name);
    }

}
