package samosa.structures;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;

import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;

public abstract class InnerOuterIDTransformer {

    // inner to outer node ids
    protected String[] map__inner_outer_id;
    // outer to inner node ids
    protected Object2IntOpenHashMap<String> map__outer_inner_id;

    protected InnerOuterIDTransformer() {
        this.map__inner_outer_id = null;
        this.map__outer_inner_id = null;
    }

    protected InnerOuterIDTransformer(final InnerOuterIDTransformer dhg) {
        this.map__inner_outer_id = dhg.map__inner_outer_id;
        this.map__outer_inner_id = dhg.map__outer_inner_id;
    }

    /**
     * Populates the map of inner to outer ids and of outer to inner ids, and
     * returns the max id of a left node and of a right node.
     * @param input_file_name path to input hypergraph
     * @return array with the max id of a left node and the max id of a right node
     * @throws Exception
     */
    protected int[] parseInputFileForCreatingOuterInnerBijection(final String input_file_name) throws Exception {
        //
        ObjectOpenHashSet<String> set__outher_ids = new ObjectOpenHashSet();
        BufferedReader br = new BufferedReader(new FileReader(input_file_name));
        String line;
        int max_id_r = -1;
        String[][] array_s_0__s_1 = new String[2][];
        String[][] array_s_0__s_1_v2 = new String[1][];
        while ((line = br.readLine()) != null) {
            max_id_r++;
            String[] s_0__s_1 = line.split("\t");
            if (s_0__s_1.length == 2) {
                String[] s_0 = s_0__s_1[0].split(",");
                String[] s_1 = s_0__s_1[1].split(",");
                array_s_0__s_1[0] = s_0;
                array_s_0__s_1[1] = s_1;
            }
            if (s_0__s_1.length == 1) {
                String[] s_0 = s_0__s_1[0].split(",");
                array_s_0__s_1 = array_s_0__s_1_v2;
                array_s_0__s_1[0] = s_0;
            }
            //
            for (String[] s_x : array_s_0__s_1) {
                for (String id_as_string : s_x) {
                    set__outher_ids.add(id_as_string);
                }
            }
        }
        br.close();
        //
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
        System.out.println("max_id_l = " + max_id_l);
        System.out.println("max_id_r = " + max_id_r);
        System.out.println();
        //
        int[] max_id_l__max_id_r = new int[2];
        max_id_l__max_id_r[0] = max_id_l;
        max_id_l__max_id_r[1] = max_id_r;
        return max_id_l__max_id_r;
    }

    /**
     * Generates a string representing the head/tail of a hyperedge using the
     * outer node ids.
     * @param head_or_tail array of inner node ids
     * @param sb StringBuilder
     * @return a string representing a comma-separated list of outer ids 
     */
    protected String toStringHeadOrTailOfHyperedge(int[] head_or_tail, StringBuilder sb) {
        //
        sb.setLength(0);
        //
        final int max_index = head_or_tail.length - 1;
        for (int index = 0; index <= max_index; index++) {
            sb.append(this.map__inner_outer_id[head_or_tail[index]]);
            if (index < max_index) {
                sb.append(",");
            }
        }
        Arrays.toString(head_or_tail);
        return sb.toString();
    }

}
