package samosa.utils.io;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import samosa.structures.DirectedBipartiteGraph;
import samosa.utils.Config;
import samosa.utils.Utils;

/**
 *
 * @author giulia
 */
public class Reader {

    /**
     * Method to read an undirected hypergraph from a space-separated file.
     *
     * @param fileName path to input file
     * @return list of undirected hyperedges
     * @throws IOException
     */
    public static List<int[]> readUndirectedHyperGraph(String fileName) throws IOException {

        final BufferedReader rows = new BufferedReader(new FileReader(fileName));
        List<int[]> edges = Lists.newArrayList();
        String line;
        int[] edge;
        String[] parts;

        while ((line = rows.readLine()) != null) {
            parts = line.split(" ");
            edge = new int[parts.length];
            for (int i = 0; i < parts.length; i++) {
                edge[i] = Integer.parseInt(parts[i]);
            }
            edges.add(edge);
        }
        rows.close();
        return edges;
    }

    /**
     * Method to read the samples obtained from an observed graph with max_id_l
     * left nodes and max_id_r right nodes.
     *
     * @param br BufferedReader for input file
     * @param max_id_l number of left nodes
     * @param max_id_r number of right nodes
     * @param transformer stores the outer-inner id mapping
     * @return the directed bipartite graph read from the input file
     * @throws IOException
     */
    public static DirectedBipartiteGraph readDirectedHyperGraph(
            BufferedReader br, int max_id_l, int max_id_r,
            Transformer transformer) throws IOException {
        // initialize temp arrays
        IntOpenHashSet[] tmp_left_out = new IntOpenHashSet[max_id_l + 1];
        for (int i = 0; i < tmp_left_out.length; i++) {
            tmp_left_out[i] = new IntOpenHashSet();
        }
        IntOpenHashSet[] tmp_left_in = new IntOpenHashSet[max_id_l + 1];
        for (int i = 0; i < tmp_left_in.length; i++) {
            tmp_left_in[i] = new IntOpenHashSet();
        }
        IntOpenHashSet[] tmp_right_out = new IntOpenHashSet[max_id_r + 1];
        for (int i = 0; i < tmp_right_out.length; i++) {
            tmp_right_out[i] = new IntOpenHashSet();
        }
        IntOpenHashSet[] tmp_right_in = new IntOpenHashSet[max_id_r + 1];
        for (int i = 0; i < tmp_right_in.length; i++) {
            tmp_right_in[i] = new IntOpenHashSet();
        }
        String line;
        int c_l_id;
        int c_r_id = -1;
        while ((line = br.readLine()) != null) {
            c_r_id++;
            String[] s_0__s_1 = line.split("\t");
            String[] s_0 = s_0__s_1[0].split(",");
            if (s_0__s_1.length > 1) {
                String[] s_1 = s_0__s_1[1].split(",");
                for (String s : s_1) {
                    if (!s.isEmpty()) {
                        s = s.strip();
                        c_l_id = transformer.getInnerId(s);
                        tmp_left_in[c_l_id].add(c_r_id);
                        tmp_right_out[c_r_id].add(c_l_id);
                    }
                }
            }
            for (String s : s_0) {
                if (!s.isEmpty()) {
                    s = s.strip();
                    c_l_id = transformer.getInnerId(s);
                    tmp_left_out[c_l_id].add(c_r_id);
                    tmp_right_in[c_r_id].add(c_l_id);
                }
            }
        }
        br.close();
        int size_D0, size_D1;
        int[][] left_P = new int[max_id_l + 1][];
        size_D0 = Utils.convertHashSetIntoArray(tmp_left_out, left_P);
        tmp_left_out = null;
        int[][] left_M = new int[max_id_l + 1][];
        Utils.convertHashSetIntoArray(tmp_left_in, left_M);
        tmp_left_in = null;
        int[][] right_P = new int[max_id_r + 1][];
        size_D1 = Utils.convertHashSetIntoArray(tmp_right_out, right_P);
        tmp_right_out = null;
        int[][] right_M = new int[max_id_r + 1][];
        Utils.convertHashSetIntoArray(tmp_right_in, right_M);
        tmp_right_in = null;
        System.gc();
        return new DirectedBipartiteGraph(left_P, left_M, right_P, right_M, size_D0, size_D1);
    }

    /**
     * Method to read the observed directed hypergraph from a tab-separated file.
     * 
     * @param in_file_name path to input file
     * @param transformer transformer to store the outer-inner id mapping
     * @return the directed bipartite graph read from the input file
     * @throws IOException
     */
    public static DirectedBipartiteGraph readDirectedHyperGraph(
            String in_file_name,
            Transformer transformer) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(in_file_name));
        String line;
        int max_id_l = 0;
        int max_id_r = -1;
        int c_id_l;
        String[][] array_s_0__s_1 = new String[2][];
        while ((line = br.readLine()) != null) {
            max_id_r++;
            String[] s_0__s_1 = line.split("\t");
            String[] s_0 = s_0__s_1[0].split(",");
            if (s_0__s_1.length > 1) {
                String[] s_1 = s_0__s_1[1].split(",");
                array_s_0__s_1[1] = s_1;
            }
            array_s_0__s_1[0] = s_0;
            for (String[] s_x : array_s_0__s_1) {
                if (s_x.length > 0) {
                    for (String id_as_string : s_x) {
                        if (!id_as_string.isEmpty()) {
                            id_as_string = id_as_string.strip();
                            transformer.saveMapping(id_as_string);
                            c_id_l = transformer.getInnerId(id_as_string);
                            max_id_l = (c_id_l > max_id_l ? c_id_l : max_id_l);
                        }
                    }
                }
            }
        }
        br.close();
        // read file again, store outer node ids
        br = new BufferedReader(new FileReader(in_file_name));
        return readDirectedHyperGraph(br, max_id_l, max_id_r, transformer);
    }

    /**
     * 
     * @param complete_file_name path to file
     * @return map of field->value pairs extracted from the string complete_file_name
     */
    public static Map<String, String> parseOutputFileName(String complete_file_name) {
        Map<String, String> map__field__value = Maps.newHashMap();
        if (complete_file_name.contains("/")) {
            String[] elems = complete_file_name.split("/");
            complete_file_name = elems[elems.length - 1];
        }
        String samplerName = Config.samplerType.split("_")[0];
        if (samplerName.equalsIgnoreCase("NuDHy")) {
            String[] tokens = complete_file_name.split("__");
            map__field__value.put("Dataset", tokens[0]);
            for (int i = 1; i < tokens.length - 1; i++) {
                String[] field__value = tokens[i].split("_");
                String value = field__value[1];
                if (field__value[0].equalsIgnoreCase("algorithm")) {
                    value = field__value[1] + "_" + field__value[2];
                }
                map__field__value.put(field__value[0], value);
            }
        }
        else {
            String[] tokens = complete_file_name.split("_" + Config.samplerType);
            map__field__value.put("Dataset", tokens[0]);
            String[] tokens2 = tokens[1].split("_");
            int seed = extractNumber(tokens2[0]);
            map__field__value.put("randomSeed", String.valueOf(seed));
            map__field__value.put("algorithm", Config.samplerType);
        } 
        return map__field__value;
    }
    
    /**
     * @param inputString a string
     * @return the first number appearing in the string, if any.
     */
    public static int extractNumber(String inputString) {
        // Define a regular expression pattern for extracting numbers
        Pattern pattern = Pattern.compile("\\d+");
        // Use a Matcher to find the first match in the input string
        Matcher matcher = pattern.matcher(inputString);
        // If there is a match, convert it to an integer and return
        if (matcher.find()) {
            return Integer.parseInt(matcher.group());
        } else {
            return -1;
        }
    }

    /**
     * Method used to read the results of StructuralEntropy from disk.
     * 
     * @param resultsPath path to file
     * @param type type of occurrences to read (accepted values are "head" and "tail")
     * @return a set of node groups represented as strings
     * @throws IOException 
     */
    public static ObjectOpenHashSet<String> readOccurrencesOfSample(
            String resultsPath,
            String type) throws IOException {

        final BufferedReader rows = new BufferedReader(new FileReader(resultsPath + "__" + type + ".tsv"));
        ObjectOpenHashSet<String> subsets = new ObjectOpenHashSet();
        String line;

        while ((line = rows.readLine()) != null) {
            subsets.add(line.replace(" ", "").strip());
        }
        rows.close();
        return subsets;
    }

    /**
     * Method used to read the results of ConvergenceTopFI from disk.
     * 
     * @param itemsetPath path to file
     * @return a map of itemset -> frequency pairs
     * @throws IOException 
     */
    public static Object2IntOpenHashMap<String> readFrequentItemsets(
            String itemsetPath) throws IOException {

        final BufferedReader rows = new BufferedReader(new FileReader(itemsetPath));
        Object2IntOpenHashMap<String> itemsets = new Object2IntOpenHashMap();
        String line;

        while ((line = rows.readLine()) != null) {
            String[] lst = line.split("#SUP:");
            String[] items = lst[0].strip().split(" ");
            int support = Integer.parseInt(lst[1].strip());
            Arrays.sort(items);
            itemsets.put(Arrays.toString(items), support);
        }
        rows.close();
        return itemsets;
    }

}
