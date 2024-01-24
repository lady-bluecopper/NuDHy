package samosa.samplers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Random;
import samosa.helpers.SimpleIntSetCountingToInfinity;
import samosa.structures.InnerOuterIDTransformer;
import samosa.utils.Config;
import samosa.utils.Timer;

public class UnpretentiousNullModel extends InnerOuterIDTransformer implements Sampler {

    // head of each hyperedge
    protected int[][] map__he_id__head;
    // tail of each hyperedge
    protected int[][] map__he_id__tail;
    // max inner id of a node in the hypergraph
    protected int max_inner_node_id;
    // structure used to initialize the heads/tails of the hyperedges
    protected SimpleIntSetCountingToInfinity set_node_id;
    // whether the head and tail of each hyperedge must be disjoint
    protected boolean force_head_tail_disjointness;

    public UnpretentiousNullModel(boolean force_head_tail_disjointness) {
        super();
        this.map__he_id__head = null;
        this.map__he_id__tail = null;
        this.max_inner_node_id = -1;
        this.set_node_id = null;
        this.force_head_tail_disjointness = force_head_tail_disjointness;
    }

    public UnpretentiousNullModel(final UnpretentiousNullModel unm) {
        super(unm);
        this.map__he_id__head = new int[unm.map__he_id__head.length][];
        for (int i = 0; i < this.map__he_id__head.length; i++) {
            this.map__he_id__head[i] = Arrays.copyOf(unm.map__he_id__head[i], unm.map__he_id__head[i].length);
        }
        this.map__he_id__tail = new int[unm.map__he_id__tail.length][];
        for (int i = 0; i < this.map__he_id__tail.length; i++) {
            if (unm.map__he_id__tail[i] != null) {
                this.map__he_id__tail[i] = Arrays.copyOf(unm.map__he_id__tail[i], unm.map__he_id__tail[i].length);
            } else {
                this.map__he_id__tail[i] = new int[0];
            }
        }
        this.max_inner_node_id = unm.max_inner_node_id;
        this.set_node_id = new SimpleIntSetCountingToInfinity(unm.set_node_id);
        this.force_head_tail_disjointness = unm.force_head_tail_disjointness;
    }
    
    /**
     * Starting from the original hypergraph, generate a random hypergraph with
     * the same head size and tail size distributions.
     * @param fileName path to original hypergraph
     * @param numSwaps parameter not used by this sampler
     * @param seed seed for reproducibility
     * @param timer object to store the running times
     * @param save whether the sample should be saved on disk
     * @throws java.lang.Exception
     */
    public void sample(String fileName, int numSwaps, 
            long seed, Timer timer, boolean save) throws Exception {

        // parse starting hypergraph
        parseInputFile(fileName);
        
        Random rnd = new Random(seed);
        int upper_bound_for_sampler = this.max_inner_node_id + 1;
        long excTimemsec;
        final long t0 = System.currentTimeMillis();
        int c_he_id = -1;
        int[] c_head = null;
        int[] c_tail = null;
        int c_node_id = -1;
        int index = -1;
        // iterate over the hyperedges to fill
        timer.start();
        for (c_he_id = 0; c_he_id < this.map__he_id__head.length; c_he_id++) {
            c_head = this.map__he_id__head[c_he_id];
            c_tail = this.map__he_id__tail[c_he_id];
            // initialize node group
            this.set_node_id.clear();
            for (index = 0; index < c_head.length; index++) {
                do {
                    // search for a new node to populate the head
                    c_node_id = rnd.nextInt(upper_bound_for_sampler);
                } while (!this.set_node_id.informativeAdd(c_node_id));
                // store node c_node_id in the head
                c_head[index] = c_node_id;
            }
            // if head and tail must be disjoint
            if (!force_head_tail_disjointness) {
                this.set_node_id.clear();
            }
            for (index = 0; index < c_tail.length; index++) {
                do {
                    // search for a new node to populate the tail
                    c_node_id = rnd.nextInt(upper_bound_for_sampler);
                } while (!this.set_node_id.informativeAdd(c_node_id));
                // store node c_node_id in the tail
                c_tail[index] = c_node_id;
            }
        }
        timer.stop();
        excTimemsec = System.currentTimeMillis() - t0;
        if (save) {
            printOnFile(fileName, (int) seed, excTimemsec, Config.resultsDir);
        }
    }

    /**
     * Save hypergraph to disk.
     * @param input_file_complete_name path to original hypergraph
     * @param random_seed seed used to generate the sample
     * @param execTimeINmsec time spent to generate the sample
     * @param output_directory directory to store the random hypergraph
     * @throws Exception 
     */
    public void printOnFile(String input_file_complete_name, 
            int random_seed,
            long execTimeINmsec, 
            String output_directory)
            throws Exception {
        // create output file name
        final int index = this.getClass().toString().lastIndexOf(".");
        String output_file_name = input_file_complete_name.replace(".tsv", "")
                .substring(input_file_complete_name.lastIndexOf('/') + 1);
        String final_out_file_name = output_directory + "/" + output_file_name + "_"
                + this.getClass().toString().substring(index + 1) + random_seed + ".tsv";
        System.out.println(final_out_file_name);
        //
        StringBuilder sb_head = new StringBuilder();
        StringBuilder sb_tail = new StringBuilder();
        BufferedWriter bw = new BufferedWriter(new FileWriter(final_out_file_name));
        String out_head, out_tail;
        for (int c_node_id = 0; c_node_id < this.map__he_id__head.length; c_node_id++) {
            // string representing the head of this hyperedge
            out_head = toStringHeadOrTailOfHyperedge(this.map__he_id__head[c_node_id], sb_head);
            // string representing the tail of this hyperedge
            out_tail = toStringHeadOrTailOfHyperedge(this.map__he_id__tail[c_node_id], sb_tail);
            bw.write(out_head);
            bw.write("\t");
            bw.write(out_tail);
            bw.write("\n");
        }
        bw.close();
    }

    /**
     * Read the original hypergraph from disk and initialize the relevant data structures.
     * @param input_file_name path to original hypergraph
     * @throws Exception
     */
    protected void parseInputFile(final String input_file_name) throws Exception {

        // Populates map of inner to outer ids and of outer to inner ids
        // and return the max id of a left node and of a right node
        int[] max_id_l__max_id_r = super.parseInputFileForCreatingOuterInnerBijection(input_file_name);
        // max id of a hypergraph node
        int max_node_id = max_id_l__max_id_r[0];
        // max id of a hyperedge
        int max_hedge_id = max_id_l__max_id_r[1];

        this.max_inner_node_id = max_node_id;
        // head of each hyperedge
        this.map__he_id__head = new int[max_hedge_id + 1][];
        // tail of each hyperedge
        this.map__he_id__tail = new int[max_hedge_id + 1][];

        BufferedReader br = new BufferedReader(new FileReader(input_file_name));
        String line;
        String[] s_0__s_1 = null;
        String[] head = null;
        String[] tail = null;
        int c_he_id = -1;
        while ((line = br.readLine()) != null) {
            c_he_id++;
            s_0__s_1 = line.split("\t");
            head = s_0__s_1[0].split(",");
            // initialize size of this head
            this.map__he_id__head[c_he_id] = new int[head.length];
            if (s_0__s_1.length == 2) {
                tail = s_0__s_1[1].split(",");
                // initialize size of this tail
                this.map__he_id__tail[c_he_id] = new int[tail.length];
            }
        }
        br.close();
        this.set_node_id = new SimpleIntSetCountingToInfinity(max_node_id);
    }

}
