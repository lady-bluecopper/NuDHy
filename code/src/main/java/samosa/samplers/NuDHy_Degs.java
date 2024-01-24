package samosa.samplers;

import java.util.Random;
import samosa.helpers.ArrayHandler;
import samosa.helpers.DirectAccessArrayHandler;
import samosa.utils.Config;
import samosa.utils.Timer;

public class NuDHy_Degs extends NuDHy {

    public NuDHy_Degs() {
        super();
    }

    public NuDHy_Degs(final NuDHy_Degs obj) {
        super(obj);
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
        ArrayHandler ah = new ArrayHandler();
        int u = 0, z = 0, v = 0, w = 0;
        int[][] P = null;
        int[][] M = null;
        for (int iteration = 1; iteration <= max_number_of_iterations; iteration++) {
            if (timer != null) {
                timer.start();
            }
            // LET'S TOSS A BIASED COIN
            if (rnd.nextDouble() <= bias) {
                // swap edges from left to right
                P = left_P;
                M = right_M;
            } else {
                // swap edges from right to left
                P = right_P;
                M = left_M;
            }
            // sample a pair of vertices
            do {
                u = rnd.nextInt(P.length);
                z = rnd.nextInt(P.length);
            } while (u == z);
            // sample a neighbor of u not in P[z] and a neighbor of z not in P[u]
            // (if any)
            ah.pickTwoRandomElementsInTheXor(P[u], P[z], rnd);
            if (ah.result[0] == -1) {
                // no neighbors to swap; we self-loop
                number_of_rejections++;
                if (timer != null) {
                    timer.stop();
                }
                continue;
            }
            // perform the swap!
            v = ah.result[0];
            w = ah.result[1];
            ah.fastReplaceInSortedArray(P[u], v, w);
            ah.fastReplaceInSortedArray(M[v], u, z);
            ah.fastReplaceInSortedArray(P[z], w, v);
            ah.fastReplaceInSortedArray(M[w], z, u);
            if (timer != null) {
                timer.stop();
            }
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

        Random rnd = new Random(seed);
        final double bias = this.rectifyBias(0.5, size_D0, size_D1);
        int number_of_rejections = 0;

        DirectAccessArrayHandler ah = new DirectAccessArrayHandler();
        int u = 0, z = 0, v = 0, w = 0;
        int index_of_v_in_Pu = 0, index_of_w_in_Pz = 0;
        int[][] P = null;
        int[][] M = null;
        final int flush = max_number_of_iterations / ((int) Config.maxNumSwapsFactor * 2);
        long excTimemsec;
        final long t0 = System.currentTimeMillis();
        for (int iteration = 0; iteration <= max_number_of_iterations; iteration++) {
            if (iteration % flush == 0) {
                excTimemsec = System.currentTimeMillis() - t0;
                super.printOnFile(file_name, iteration, seed, excTimemsec,
                        number_of_rejections, Config.resultsDir);
            }
            // LET'S TOSS A BIASED COIN
            if (rnd.nextDouble() <= bias) {
                // swap edges from left to right
                P = this.left_P;
                M = this.right_M;
            } else {
                // swap edges from right to left
                P = this.right_P;
                M = this.left_M;
            }
            // sample a pair of vertices
            do {
                u = rnd.nextInt(P.length);
                z = rnd.nextInt(P.length);
            } while (u == z);
            // sample a neighbor of u not in P[z] and a neighbor of z not in P[u]
            // (if any)
            ah.pickTwoRandomElementsInTheXor(P[u], P[z], rnd);
            if (ah.result[0] == -1) {
                // no neighbors to swap; we self-loop
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
     * Performs numSwaps swaps in the Markov graph to generate a random sample
     * starting from dirHyperGraph
     *
     * @param fileName path to input file
     * @param numSwaps number of swaps
     * @param seed long number
     * @param timer Timer
     * @param save sample should be saved on disk
     * @throws java.lang.Exception
     */
    @Override
    public void sample(String fileName, int numSwaps, long seed, Timer timer, boolean save) throws Exception {

        parseInputFile(fileName);

        final int int_seed = (int) seed;
        int number_of_rejections;
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
