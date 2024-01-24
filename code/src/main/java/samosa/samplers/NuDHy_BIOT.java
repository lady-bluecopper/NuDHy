package samosa.samplers;

import java.util.Arrays;
import java.util.Random;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import samosa.utils.Config;
import samosa.utils.Timer;
import samosa.helpers.ArrayHandler;
import samosa.helpers.DirectAccessArrayHandler;
import samosa.utils.Utils;


public class NuDHy_BIOT extends NuDHy {

    // probability of sampling each combo of in- and out-degree of a left node
    // when looking for edges from left to right to swap
    protected long[] Lp__map__event_id__probability = null;
    // candidate left nodes to select for each combo
    protected int[][] Lp__map__event_id__subpopulation = null;

    // probability of sampling each combo of in- and out-degree of a left node
    // when looking for edges from right to left to swap
    protected long[] Ln__map__event_id__probability = null;
    // candidate left nodes to select for each combo
    protected int[][] Ln__map__event_id__subpopulation = null;

    // probability of sampling each combo of in- and out-degree of a right node
    // when looking for edges from left to right to swap
    protected long[] Rp__map__event_id__probability = null;
    // candidate right nodes to select for each combo
    protected int[][] Rp__map__event_id__subpopulation = null;

    // probability of sampling each combo of in- and out-degree of a right node
    // when looking for edges from right to left to swap
    protected long[] Rn__map__event_id__probability = null;
    // candidate right nodes to select for each combo
    protected int[][] Rn__map__event_id__subpopulation = null;

    public NuDHy_BIOT() {
        super();
    }

    public NuDHy_BIOT(final NuDHy_BIOT obj) {
        super(obj);
        //
        this.Lp__map__event_id__probability = Arrays.copyOf(obj.Lp__map__event_id__probability,
                obj.Lp__map__event_id__probability.length);
        this.Ln__map__event_id__probability = Arrays.copyOf(obj.Ln__map__event_id__probability,
                obj.Ln__map__event_id__probability.length);
        this.Rp__map__event_id__probability = Arrays.copyOf(obj.Rp__map__event_id__probability,
                obj.Rp__map__event_id__probability.length);
        this.Rn__map__event_id__probability = Arrays.copyOf(obj.Rn__map__event_id__probability,
                obj.Rn__map__event_id__probability.length);
        //
        this.Lp__map__event_id__subpopulation = new int[obj.Lp__map__event_id__subpopulation.length][];
        this.Ln__map__event_id__subpopulation = new int[obj.Ln__map__event_id__subpopulation.length][];
        this.Rp__map__event_id__subpopulation = new int[obj.Rp__map__event_id__subpopulation.length][];
        this.Rn__map__event_id__subpopulation = new int[obj.Rn__map__event_id__subpopulation.length][];
        Utils.initializeArrayOfArrays(obj.Lp__map__event_id__subpopulation, this.Lp__map__event_id__subpopulation);
        Utils.initializeArrayOfArrays(obj.Ln__map__event_id__subpopulation, this.Ln__map__event_id__subpopulation);
        Utils.initializeArrayOfArrays(obj.Rp__map__event_id__subpopulation, this.Rp__map__event_id__subpopulation);
        Utils.initializeArrayOfArrays(obj.Rn__map__event_id__subpopulation, this.Rn__map__event_id__subpopulation);
    }

    /**
     * Computes vartheta(i, j), nu(i,j), eta(i,j) and phi(i,j) for each combo
     * of in- and out-degrees, and stores the candidate nodes to select for each combo
     */
    protected void createMappingsEventsProbabilities() {
        boolean compute_on_LEFT;
        // Compute vartheta(i, j) and eta(i,j)
        compute_on_LEFT = true;
        this.computeLorRpn(this.left_P, this.left_M, compute_on_LEFT);
        // Compute nu(i,j) and phi(i,j)
        compute_on_LEFT = false;
        this.computeLorRpn(this.right_P, this.right_M, compute_on_LEFT);
        System.out.println();
        System.out.println();
        System.out.println("this.Lp__map__event_id__probability.length " + this.Lp__map__event_id__probability.length);
        System.out.println(
                "this.Lp__map__event_id__subpopulation.length " + this.Lp__map__event_id__subpopulation.length);
        System.out.println("this.Ln__map__event_id__probability.length " + this.Ln__map__event_id__probability.length);
        System.out.println(
                "this.Ln__map__event_id__subpopulation.length " + this.Ln__map__event_id__subpopulation.length);
        System.out.println("this.Rp__map__event_id__probability.length " + this.Rp__map__event_id__probability.length);
        System.out.println(
                "this.Rp__map__event_id__subpopulation.length " + this.Rp__map__event_id__subpopulation.length);
        System.out.println("this.Rn__map__event_id__probability.length " + this.Rn__map__event_id__probability.length);
        System.out.println(
                "this.Rn__map__event_id__subpopulation.length " + this.Rn__map__event_id__subpopulation.length);
        System.out.println();
        System.out.println();

    }

    /**
     * Hash a pair of in- and out-degrees into a long.
     * @param id_l in-degree
     * @param id_r out-degree
     * @param shift shift
     * @return 
     */
    protected static long computeHash(final long id_l, final long id_r, final int shift) {
        return (id_l << shift) + id_r;
    }

    /**
     * Given a hash code, retrieves the corresponding in- and out-degrees.
     * @param key a hash code
     * @param shift shift
     * @param right_mask_shift mask
     * @param result__id_l__id_r array to store the in- and out-degree
     */
    protected static void computeReverseHash(final long key, final int shift, final long right_mask_shift,
            final long[] result__id_l__id_r) {
        result__id_l__id_r[1] = key & right_mask_shift;
        result__id_l__id_r[0] = key >> shift;
    }

    /**
     * Computes vartheta(i, j)/eta(i,j) and nu(i,j)/phi(i,j) for each combo
     * of in- and out-degrees of a left/right node, and stores the candidate 
     * nodes to select for each combo.
     * @param P out-neighbors of each node (either left or right)
     * @param N in-neighbors of each node (either left or right)
     * @param compute_on_LEFT
     */
    protected void computeLorRpn(int[][] P, int[][] N, boolean compute_on_LEFT) {
        // max in-degree
        int max_in_deg = -1;
        // max out-degree
        int max_out_deg = -1;
        for (int node_id = 0; node_id < P.length; node_id++) {
            max_in_deg = (N[node_id].length > max_in_deg ? N[node_id].length : max_in_deg);
            max_out_deg = (P[node_id].length > max_out_deg ? P[node_id].length : max_out_deg);
        }
        // initialize variables used to hash the pairs of in- and out-degrees
        int shift = 1 + ((int) Math.floor(Math.log(max_out_deg) / Math.log(2)));
        long right_mask_shift = ((long) Math.pow(2, shift)) - 1L;
        System.out.println("shift             = " + shift);
        System.out.println("right_mask_shift  = " + right_mask_shift);
        System.out.println();

        final int max = (max_in_deg >= max_out_deg ? max_in_deg : max_out_deg);
        // for each combo of in- and out-degree, list of nodes with that combo
        Long2ObjectOpenHashMap<IntArrayList> LorRpn = new Long2ObjectOpenHashMap(max);
        LorRpn.defaultReturnValue(null);
        int c_node_in_deg, c_node_out_deg;
        long key;
        IntArrayList c_set_of_nodes = null;
        for (int node_id = 0; node_id < P.length; node_id++) {
            c_node_out_deg = P[node_id].length;
            c_node_in_deg = N[node_id].length;
            // hash the combo of in- and out-degree of this node
            key = computeHash(c_node_in_deg, c_node_out_deg, shift);
            // store in the map
            if (!LorRpn.containsKey(key)) {
                c_set_of_nodes = new IntArrayList();
                LorRpn.putIfAbsent(key, c_set_of_nodes);
            } else {
                c_set_of_nodes = LorRpn.get(key);
            }
            c_set_of_nodes.add(node_id);
        }
        //
        // compute population and sub-populations.
        long p_denominator = 0;
        long n_denominator = 0;
        LongArrayList p_map__event_id__probability = new LongArrayList();
        LongArrayList n_map__event_id__probability = new LongArrayList();
        ObjectArrayList<int[]> temp_p_map__event_id__subpopulation = new ObjectArrayList();
        ObjectArrayList<int[]> temp_n_map__event_id__subpopulation = new ObjectArrayList();
        long c_pn_size = 0;
        long[] inDeg_outDeg = {-1, -1};
        IntArrayList c_subpopulation;
        
        for (long c_key : LorRpn.keySet()) {
            c_subpopulation = LorRpn.get(c_key);
            // number of nodes with that combo of in- and out-degree
            c_pn_size = c_subpopulation.size();
            // we need at least two nodes with the same combo to swap edges
            // between them
            if (c_pn_size >= 2) {
                computeReverseHash(c_key, shift, right_mask_shift, inDeg_outDeg);
                // if the out-degree is > 0, we can consider this key when
                // looking for candidate left to right edges to swap
                if (inDeg_outDeg[1] >= 1) {
                    // update denominator of probability of sampling the combo c_key
                    // when looking for edges from left to right to swap
                    // this is the denominator of vartheta(i, j) and nu(i,j)                   
                    p_denominator += (c_pn_size * (c_pn_size - 1L)) / 2L;
                    p_map__event_id__probability.add(p_denominator);
                    // store the list of nodes from which the pair is drawn 
                    // when the combo c_key is sampled
                    temp_p_map__event_id__subpopulation.add(c_subpopulation.toIntArray());
                }
                // if the in-degree is > 0, we can consider this key when
                // looking for candidate right to left edges to swap
                if (inDeg_outDeg[0] >= 1) {
                    // update denominator of probability of sampling the combo c_key
                    // when looking for edges from right to left to swap
                    // this is the denominator of eta(i,j) and phi(i,j)
                    n_denominator += (c_pn_size * (c_pn_size - 1L)) / 2L;
                    n_map__event_id__probability.add(n_denominator);
                    // store the list of nodes from which the pair is drawn 
                    // when the combo c_key is sampled
                    temp_n_map__event_id__subpopulation.add(c_subpopulation.toIntArray());
                }
            }
        }
        //
        if (compute_on_LEFT) {
            // denominator of vartheta(i, j)
            this.Lp__map__event_id__probability = p_map__event_id__probability.toLongArray();
            // convert ArrayList to array of arrays
            this.Lp__map__event_id__subpopulation = new int[temp_p_map__event_id__subpopulation.size()][];
            for (int i = 0; i < temp_p_map__event_id__subpopulation.size(); i++) {
                this.Lp__map__event_id__subpopulation[i] = temp_p_map__event_id__subpopulation.get(i);
            }
            // denominator of eta(i,j)
            this.Ln__map__event_id__probability = n_map__event_id__probability.toLongArray();
            // convert ArrayList to array of arrays
            this.Ln__map__event_id__subpopulation = new int[temp_n_map__event_id__subpopulation.size()][];
            for (int i = 0; i < temp_n_map__event_id__subpopulation.size(); i++) {
                this.Ln__map__event_id__subpopulation[i] = temp_n_map__event_id__subpopulation.get(i);
            }
        } else {
            // denominator of nu(i,j)
            this.Rp__map__event_id__probability = p_map__event_id__probability.toLongArray();
            // convert ArrayList to array of arrays
            this.Rp__map__event_id__subpopulation = new int[temp_p_map__event_id__subpopulation.size()][];
            for (int i = 0; i < temp_p_map__event_id__subpopulation.size(); i++) {
                this.Rp__map__event_id__subpopulation[i] = temp_p_map__event_id__subpopulation.get(i);
            }
            // denominator of phi(i,j)
            this.Rn__map__event_id__probability = n_map__event_id__probability.toLongArray();
            // convert ArrayList to array of arrays
            this.Rn__map__event_id__subpopulation = new int[temp_n_map__event_id__subpopulation.size()][];
            for (int i = 0; i < temp_n_map__event_id__subpopulation.size(); i++) {
                this.Rn__map__event_id__subpopulation[i] = temp_n_map__event_id__subpopulation.get(i);
            }
        }
        // DEBUG ///////////////////////////////////////////////////////////////
        int c_node_id, p_node_in_deg, p_node_out_deg;
        int[][][] all_subpopulations = new int[2][][];
        long[][] all_probabilities = new long[2][];
        if (compute_on_LEFT) {
            all_subpopulations[0] = this.Lp__map__event_id__subpopulation;
            all_subpopulations[1] = this.Ln__map__event_id__subpopulation;
            all_probabilities[0] = this.Lp__map__event_id__probability;
            all_probabilities[1] = this.Ln__map__event_id__probability;
        } else {
            all_subpopulations[0] = this.Rp__map__event_id__subpopulation;
            all_subpopulations[1] = this.Rn__map__event_id__subpopulation;
            all_probabilities[0] = this.Rp__map__event_id__probability;
            all_probabilities[1] = this.Rn__map__event_id__probability;
        }
        for (int index = 0; index < 2; index++) {
            int[][] c__map__event_id__subpopulation = all_subpopulations[index];
            long[] c__probabilities = all_probabilities[index];
            if (c__probabilities == null) {
                continue;
            }
            System.out.println("Check");
            System.out.println("index " + index);
            for (int i = 0; i < c__map__event_id__subpopulation.length; i++) {
                p_node_in_deg = -1;
                p_node_out_deg = -1;
                for (int j = 0; j < c__map__event_id__subpopulation[i].length; j++) {
                    c_node_id = c__map__event_id__subpopulation[i][j];
                    c_node_in_deg = N[c_node_id].length;
                    c_node_out_deg = P[c_node_id].length;
                    if (j > 0) {
                        // assert
                        if (!(p_node_in_deg == c_node_in_deg && p_node_out_deg == c_node_out_deg)) {
                            System.out.println("ERR");
                            System.exit(-1);
                        }
                    }
                    p_node_in_deg = c_node_in_deg;
                    p_node_out_deg = c_node_out_deg;
                }
            }
        }
        // END DEBUG ///////////////////////////////////////////////////////////
    }

    /**
     * Moves in the Markov graph for max_number_of_iterations steps.
     *
     * @param max_number_of_iterations number of iterations
     * @param seed seed for reproducibility
     * @param timer Timer to store running times
     * @return number of self-loops performed
     */
    @Override
    public int runMC(final int max_number_of_iterations, final int seed, Timer timer) {
        
        Random rnd = new Random(seed);
        final double bias = this.rectifyBias(0.5, size_D0, size_D1);
        int number_of_rejections = 0;
        //
        ArrayHandler ah = new ArrayHandler();
        int u = 0, z = 0, v = 0, w = 0;
        int[][] P = null;
        int[][] M = null;
        //
        long[] map__event_id__probability;
        int[][] map__event_id__subpopulation;
        //
        long dice;
        int event_id_from_dice;
        int[] sampled_subpopulation;
        boolean biased_coin, fair_coin;
        for (int iteration = 1; iteration <= max_number_of_iterations; iteration++) {
            timer.start();
            do {
                // LET'S TOSS A BIASED COIN
                biased_coin = (rnd.nextDouble() <= bias);
                // LET'S TOSS A FAIR COIN
                fair_coin = rnd.nextBoolean();
                // case <heads, heads>
                if (biased_coin == true && fair_coin == true) {
                    map__event_id__probability = this.Lp__map__event_id__probability;
                    map__event_id__subpopulation = this.Lp__map__event_id__subpopulation;
                    P = this.left_P;
                    M = this.right_M;
                // case <tails, heads>
                } else if (biased_coin == false && fair_coin == true) {
                    map__event_id__probability = this.Ln__map__event_id__probability;
                    map__event_id__subpopulation = this.Ln__map__event_id__subpopulation;
                    P = this.left_M;
                    M = this.right_P;
                // case <heads, tails>
                } else if (biased_coin == true && fair_coin == false) {
                    map__event_id__probability = this.Rn__map__event_id__probability;
                    map__event_id__subpopulation = this.Rn__map__event_id__subpopulation;
                    P = this.right_M;
                    M = this.left_P;
                // case <tails, tails>
                } else {
                    map__event_id__probability = this.Rp__map__event_id__probability;
                    map__event_id__subpopulation = this.Rp__map__event_id__subpopulation;
                    P = this.right_P;
                    M = this.left_M;
                }
            } while (map__event_id__probability.length == 0);
            // sample a combo of (in, out) degree
            dice = 1 + rnd.nextLong(map__event_id__probability[map__event_id__probability.length - 1]);
            event_id_from_dice = ah.getIndexToAddValueKeepingTheArraySorted(map__event_id__probability, dice);
            sampled_subpopulation = map__event_id__subpopulation[event_id_from_dice];
            // sample a pair of vertices with that combo of in- and out-degrees
            do {
                u = sampled_subpopulation[rnd.nextInt(sampled_subpopulation.length)];
                z = sampled_subpopulation[rnd.nextInt(sampled_subpopulation.length)];
            } while (u == z);
            // sample a neighbor of u not in P[z] and a neighbor of z not in P[u]
            // (if any)
            ah.pickTwoRandomElementsInTheXor(P[u], P[z], rnd);
            if (ah.result[0] == -1) {
                // no neighbor to swap; we self-loop
                number_of_rejections++;
                timer.stop();
                continue;
            }
            // perform the swap!
            v = ah.result[0];
            w = ah.result[1];
            ah.fastReplaceInSortedArray(P[u], v, w);
            ah.fastReplaceInSortedArray(M[v], u, z);
            ah.fastReplaceInSortedArray(P[z], w, v);
            ah.fastReplaceInSortedArray(M[w], z, u);
            timer.stop();
        }
        return number_of_rejections;
    }

    /**
     * Starting from the original hypergraph, move in the Markov graph and store
     * on disk the current state after k steps up to max_number_of_iterations steps.
     * @param max_number_of_iterations maximum number of steps in the Markov graph
     * @param seed seed for reproducibility
     * @param file_name path to original hypergraph
     * @return number of self-loops performed
     * @throws Exception 
     */
    @Override
    public int runMC_multiFlush(final int max_number_of_iterations,
            int seed, String file_name) throws Exception {

        parseInputFile(file_name);
        createMappingsEventsProbabilities();

        Random rnd = new Random(seed);
        final double bias = this.rectifyBias(0.5, size_D0, size_D1);
        int number_of_rejections = 0;

        DirectAccessArrayHandler ah = new DirectAccessArrayHandler();
        int u = 0, z = 0, v = 0, w = 0;
        int index_of_v_in_Pu = 0, index_of_w_in_Pz = 0;
        int[][] P = null;
        int[][] M = null;
        //
        long[] map__event_id__probability;
        int[][] map__event_id__subpopulation;
        //
        long dice;
        int event_id_from_dice;
        int[] sampled_subpopulation;
        boolean biased_coin, fair_coin;
        final int flush = max_number_of_iterations / ((int) Config.maxNumSwapsFactor * 2);

        long excTimemsec;
        final long t0 = System.currentTimeMillis();
        for (int iteration = 0; iteration <= max_number_of_iterations; iteration++) {
            if (iteration % flush == 0) {
                excTimemsec = System.currentTimeMillis() - t0;
                super.printOnFile(file_name, iteration, seed, excTimemsec,
                        number_of_rejections, Config.resultsDir);
            }
            do {
                // LET'S TOSS A BIASED COIN
                biased_coin = (rnd.nextDouble() <= bias);
                // LET'S TOSS A FAIR COIN
                fair_coin = rnd.nextBoolean();
                // case <heads, heads>
                if (biased_coin == true && fair_coin == true) {
                    map__event_id__probability = this.Lp__map__event_id__probability;
                    map__event_id__subpopulation = this.Lp__map__event_id__subpopulation;
                    P = this.left_P;
                    M = this.right_M;
                // case <tails, heads>
                } else if (biased_coin == false && fair_coin == true) {
                    map__event_id__probability = this.Ln__map__event_id__probability;
                    map__event_id__subpopulation = this.Ln__map__event_id__subpopulation;
                    P = this.left_M;
                    M = this.right_P;
                // case <heads, tails>
                } else if (biased_coin == true && fair_coin == false) {
                    map__event_id__probability = this.Rn__map__event_id__probability;
                    map__event_id__subpopulation = this.Rn__map__event_id__subpopulation;
                    P = this.right_M;
                    M = this.left_P;
                // case <tails, tails>
                } else {
                    map__event_id__probability = this.Rp__map__event_id__probability;
                    map__event_id__subpopulation = this.Rp__map__event_id__subpopulation;
                    P = this.right_P;
                    M = this.left_M;
                }
            } while (map__event_id__probability.length == 0);
            // sample a combo of (in, out) degree.
            dice = 1 + rnd.nextLong(map__event_id__probability[map__event_id__probability.length - 1]);
            event_id_from_dice = ah.getIndexToAddValueKeepingTheArraySorted(map__event_id__probability, dice);
            sampled_subpopulation = map__event_id__subpopulation[event_id_from_dice];
            // sample a pair of vertices with that combo of in- and out-degrees
            do {
                u = sampled_subpopulation[rnd.nextInt(sampled_subpopulation.length)];
                z = sampled_subpopulation[rnd.nextInt(sampled_subpopulation.length)];
            } while (u == z);
            // sample a neighbor of u not in P[z] and a neighbor of z not in P[u]
            // (if any)
            ah.pickTwoRandomElementsInTheXor(P[u], P[z], rnd);
            if (ah.result[0] == -1) {
                // no neighbor to swap; we self-loop
                number_of_rejections++;
                continue;
            }
            // perform the swap!
            v = ah.result[0];
            index_of_v_in_Pu = ah.result_indexes[0];
            w = ah.result[1];
            index_of_w_in_Pz = ah.result_indexes[1];
            ah.fastReplaceInSortedArray(P[u], v, w, index_of_v_in_Pu);
            ah.fastReplaceInSortedArray(M[v], u, z);
            ah.fastReplaceInSortedArray(P[z], w, v, index_of_w_in_Pz);
            ah.fastReplaceInSortedArray(M[w], z, u);
        }
        return number_of_rejections;
    }

    /**
     * Performs numSwaps iterations in the Markov graph to generate a random sample
     * starting from the original hypergraph read from fileName
     *
     * @param fileName path to original hypergraph
     * @param numSwaps number of iterations
     * @param seed long number
     * @param timer Timer to store running times
     * @param save whether the sample should be saved on disk
     * @throws java.lang.Exception
     */
    @Override
    public void sample(String fileName, int numSwaps, long seed, Timer timer, boolean save) throws Exception {
        // parse starting hypergraph
        parseInputFile(fileName);
        // initialize probabilities
        createMappingsEventsProbabilities();
        // initial bias
        final int int_seed = (int) seed;
        int number_of_rejections;
        boolean compute_on_LEFT;
        compute_on_LEFT = true;
        computeLorRpn(left_P, left_M, compute_on_LEFT);
        compute_on_LEFT = false;
        computeLorRpn(right_P, right_M, compute_on_LEFT);
        // start Markov chain
        long t0 = System.currentTimeMillis();
        number_of_rejections = runMC(numSwaps, int_seed, timer);
        long t1 = System.currentTimeMillis();
        System.out.println("---------------------------------------------");
        System.out.println((t1 - t0) + "msec");
        System.out.println(" max_number_of_iterations = " + numSwaps);
        System.out.println(" number_of_rejections     = " + number_of_rejections);
        System.out.println(" number_of_swaps          = " + (numSwaps - number_of_rejections));
        System.out.println(" AcceptanceRate           = " + (((float) numSwaps - number_of_rejections) / numSwaps));
        System.out.println("---------------------------------------------");

        if (save) {
            super.printOnFile(fileName, numSwaps, int_seed, (t1 - t0), number_of_rejections,
                    Config.resultsDir);
        }
    }

}
