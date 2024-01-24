package samosa.utils.io;

import com.google.common.collect.Lists;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.javatuples.Pair;
import samosa.utils.Config;

/**
 *
 * @author giulia
 */
public class Writer {

    /**
     * Write a multiset of transactions to disk.
     * The elements in each transaction are separated by a space.
     * 
     * @param trans transaction database
     * @param fName output filename
     * @throws IOException 
     */
    public static void writeTransactions(int[][] trans, String fName) throws IOException {

        FileWriter fwP = new FileWriter(fName);
        StringBuilder builder;
        for (int[] t : trans) {
            if (t.length > 0) {
                builder = new StringBuilder();
                for (int v : t) {
                    builder.append(v);
                    builder.append(" ");
                }
                builder.deleteCharAt(builder.length() - 1);
                builder.append("\n");
                fwP.write(builder.toString());
            }
        }
        fwP.close();
    }
    
    /**
     * Writes the frequent itemsets on disk.
     * @param freqItemsets for each test, map size -> list of frequent itemset of that size 
     * @param resultsPath path to output file
     * @throws IOException 
     */
    public static void writeFrequentItemsets(Map<Integer, List<Set<Integer>>>[] freqItemsets, 
            String resultsPath) throws IOException {
        
        FileWriter fwP = new FileWriter(resultsPath);
        String header = "Test\tSize\tRank\tFrequentItemset\n";
        List<String> row;
        try {
            // save header of file
            fwP.write(header);
            for (int c = 0; c < freqItemsets.length; c++) {
                for (int size : freqItemsets[c].keySet()) {
                    for (int p = 0; p < freqItemsets[c].get(size).size(); p++) {
                        row = Lists.newArrayList();
                        row.add(String.valueOf(c));
                        row.add(String.valueOf(size));
                        row.add(String.valueOf(p));
                        ObjectArrayList fi = new ObjectArrayList();
                        for (int i : freqItemsets[c].get(size).get(p)) {
                            fi.add(String.valueOf(i));
                        }
                        Collections.sort(fi);
                        row.add(fi.toString());
                        fwP.write(String.join("\t", row) + "\n");
                    }
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        fwP.close();
    }
    
    /**
     * Write the number of occurrences in the samples of all the node groups of
     * size less than 4 in the head and tail of the observed hyperedges. 
     * @param resultsPath path to output file
     * @param type group type (either "head" or "tail")
     * @param allSubsets set of all the node groups in the observed hypergraph
     * @param occurrences map group -> number of samples that contain the group
     * in the head/tail of some hyperedge
     * @param sampler sampler name
     * @param append whether the results should be appended to an existing output 
     * file or the file must be created from scratch
     * @throws IOException 
     */
    public static void writeStructuralProbabilities(String resultsPath, 
            String type, ObjectOpenHashSet<String> allSubsets,
            Object2DoubleOpenHashMap<String> occurrences,
            String sampler, boolean append) throws IOException {
        
        FileWriter fwP = new FileWriter(resultsPath + "__" + type + ".tsv", append);
        String header = "Sampler\tSize\tProbability\n";
        List<String> row;
        try {
            // save header of file
            fwP.write(header);
            for (String subset : allSubsets) {
                row = Lists.newArrayList();
                row.add(sampler);
                row.add(String.valueOf(subset.split(",").length));
                row.add(String.valueOf(occurrences.getDouble(subset) / Config.numSamples));
                fwP.write(String.join("\t", row) + "\n");
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        fwP.close();
    }
    
    /**
     * Writes to disk the node groups that occur in a certain sample of an observed 
     * hypergraph.
     * 
     * @param output for each size s, node groups of size s that are present in 
     * some hyperedge head/tail of the sample
     * @param resultsPath path to output file
     * @param type group type (either "head" or "tail")
     * @throws IOException 
     */
    public static void writeOccurrencesOfSample(
            ObjectOpenHashSet<String>[] output,
            String resultsPath,
            String type) throws IOException {
        
        FileWriter fwP = new FileWriter(resultsPath + "__" + type + ".tsv");
        List<String> row;
        try {
            for (int s = 0; s < output.length; s++) {
                for (String subset : output[s]) {
                    fwP.write(subset.substring(1, subset.length() - 1) + "\n");
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        fwP.close();
    }

    /**
     * 
     * @param date time stamp
     * @param runStats list of statistics
     * @param resultsPath where to store the statistics
     * @param append whether the statistics must be appended to the result file,
     * or the result file must be re-created  
     * @throws IOException 
     */
    public static void writeStatistics(LocalDateTime date, 
            List<Map<String, String>> runStats, 
            String resultsPath,
            boolean append) throws IOException {
        
        if (runStats.isEmpty()) {
            return;
        }
        FileWriter fwP = new FileWriter(resultsPath, append);
        List<String> keys = Lists.newArrayList();
        List<String> row;
        keys.add("Date");
        try {
            // save header of file
            keys.addAll(runStats.get(0).keySet());
            String header = String.join(" ", keys);
            fwP.write(header + "\n");
            for (Map<String, String> map : runStats) {
                row = Lists.newArrayList();
                row.add(date.toString());
                for (int i = 1; i < keys.size(); i++) {
                    row.add(map.get(keys.get(i)));
                }
                fwP.write(String.join(" ", row) + "\n");
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        fwP.close();
    }
    
    /**
     * Writes a directed hypergraph to disk.
     * 
     * @param edges directed hyperedges as a collection of pairs (each pair 
     * include the head and the tail of a directed hyperedge)
     * @param outPath path to output file
     * @throws IOException 
     */
    public static void writeHypergraph(Collection<Pair<Collection<Integer>, Collection<Integer>>> edges, 
            String outPath) throws IOException {
        FileWriter fwP = new FileWriter(outPath);
        try {
            for (Pair<Collection<Integer>, Collection<Integer>> edge : edges) {
                String head = String.join(",", edge.getValue0().stream().map(i -> i.toString()).toList());
                String tail = String.join(",", edge.getValue1().stream().map(i -> i.toString()).toList());
                fwP.write(head + "\t" + tail + "\n");
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        fwP.close();
    }
    
}

