package samosa.test;

import samosa.utils.CMDLineParser;
import samosa.utils.Config;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.util.Collections;
import java.util.Enumeration;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import org.apache.commons.io.FilenameUtils;
import org.javatuples.Pair;
import samosa.fim.FreqItemsetMiner;
import samosa.structures.DirectedBipartiteGraph;
import samosa.utils.io.Reader;
import samosa.utils.io.Transformer;
import samosa.utils.io.Writer;

/**
 * This class runs the convergence experiment by measuring the average relative
 * support/frequency difference for different values of the swap number
 * multiplier/factor.
 */
public class ConvergenceFI {

    public static void main(String[] args) throws IOException {

        CMDLineParser.parse(args);

        // READ DATASET
        final String filePath = Path.of(Config.datasetsDir, Config.datasetName).toString();
        final Transformer T = new Transformer();
        final DirectedBipartiteGraph dirHyperGraph = Reader.readDirectedHyperGraph(filePath, T);
        // for each edge direction and each side of the bipartite graph
        // left-out, left-in, right-out, right-in
        final int[] combos = new int[]{0, 1};
        final String[] freqs = Config.minFreqs.split(",");
        // FIND FREQUENT ITEMSETS
        long startTime = System.currentTimeMillis();
        Map<Set<Integer>, Integer>[] allFreqItemsets = new Map[combos.length];
        for (int c : combos) {
            final int[][] transactions = dirHyperGraph.getTransactions(c);
            final double freq = Double.parseDouble(freqs[c]);
            System.out.println("Mining frequent itemsets for combo " + c + ", with frequency " + freqs[c]);
            final Map<Set<Integer>, Integer> freqItemsetToSup = FreqItemsetMiner.mine(transactions, freq);
            if (Config.k > 0) {
                Map<Set<Integer>, Integer> topKFIs = Maps.newHashMap();
                List<Entry<Set<Integer>, Integer>> tmp = Lists.newArrayList();
                freqItemsetToSup.entrySet()
                        .stream()
                        .forEach(e -> {
                            if (e.getKey().size() >= Config.size) {
                                tmp.add(e);
                            }
                        });
                if (!tmp.isEmpty()) {
                    Collections.sort(tmp, (Entry<Set<Integer>, Integer> o1, Entry<Set<Integer>, Integer> o2) 
                            -> - Integer.compare(o1.getValue(), o2.getValue()));
                    for (int id = 0; id < Math.min(Config.k, tmp.size()); id++) {
                        topKFIs.put(tmp.get(id).getKey(), tmp.get(id).getValue());
                    }
                }
                allFreqItemsets[c] = topKFIs;
            } else {
                allFreqItemsets[c] = freqItemsetToSup;
            }
            System.out.println("Num Freq. Itemset: " + allFreqItemsets[c].size());
        }
        System.out.println("\nMining done in " + (System.currentTimeMillis() - startTime));
        // INIT
        final String dataName = Config.datasetName.split("\\.")[0];
        List<Map<String, String>> convergenceStats = Lists.newArrayList();
        int counter = 0;
        final String zipPath = Path.of(Config.resultsDir, "samples", "convergence", dataName + ".zip").toString();
        ZipFile zip = new ZipFile(zipPath);
        for (Enumeration e = zip.entries(); e.hasMoreElements(); ) {
            ZipEntry entry = (ZipEntry) e.nextElement();
            String fileName = entry.getName();
            if (FilenameUtils.getExtension(fileName).equals("tsv") && 
                    fileName.contains(dataName) &&
                    !fileName.startsWith(".") &&
                    !fileName.startsWith("__MACOSX")) {
                Map<String, String> fieldValueMap = Reader.parseOutputFileName(fileName);
                if (!fieldValueMap.get("sampleNumber").equals(String.valueOf(Config.seed))) {
                    continue;
                }
                counter += 1;
                System.out.println("Running convergence for file n. " + counter);
                System.out.println(fileName);
                BufferedReader reader = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
                runConvergence(combos,
                                allFreqItemsets,
                                reader,
                                dirHyperGraph.getLeftP().length,
                                dirHyperGraph.getRightP().length,
                                T,
                                fieldValueMap);
                convergenceStats.add(fieldValueMap);
            }
        }
        // SAVE STATS
        final LocalDateTime date = LocalDateTime.now();
        String resultsBaseName = String.join("__",
                "Convergence",
                "Dataset_" + Config.datasetName,
                "sampleNumber_" + String.valueOf(Config.seed));
        if (Config.k > 0) {
            resultsBaseName = String.join("__",
                "Convergence",
                "Dataset_" + Config.datasetName,
                "sampleNumber_" + String.valueOf(Config.seed),
                "k_" + String.valueOf(Config.k),
                "size_" + String.valueOf(Config.size));
        }
        Files.createDirectories(Path.of(Config.resultsDir, "convergence"));
        final String resultsPath = Path.of(Config.resultsDir, "convergence", resultsBaseName).toString();
        Writer.writeStatistics(date, convergenceStats, resultsPath, false);
        System.out.println("Result written to " + resultsPath);
    }

    /**
     * Run the convergence experiment for the random hypergraph read by the reader.
     * @param combos which kinds of datasets of itemsets to consider
     * @param freqItemsets frequent itemsets in each kind of dataset corresponding 
     * to the observed hypergraph
     * @param reader BufferedReader to read the random hypergraph from disk
     * @param max_left_id max id of a left node
     * @param max_right_id max id of a right node
     * @param T Transformer to get the inner and outer node ids
     * @param fieldValueMap fields extracted from the file name of the random hypergraph
     * @throws IOException 
     */
    private static void runConvergence(
            int[] combos,
            Map<Set<Integer>, Integer>[] freqItemsets,
            BufferedReader reader,
            final int max_left_id,
            final int max_right_id,
            Transformer T,
            Map<String, String> fieldValueMap) throws IOException {

        // initialize statistics
        Map<String, String> stats = Maps.newHashMap();
        // load hypergraph
        DirectedBipartiteGraph dirHyperGraph = Reader.readDirectedHyperGraph(reader, max_left_id, max_right_id, T);
        for (int c : combos) {
            System.out.println("\t\tGetting sample itemset to support map for combo " + c);
            final int[][] itemMap = dirHyperGraph.getTransactions(3-c);
            final Map<Set<Integer>, Integer> sampleFreqItemsetToSup = getItemsetToSupMap(itemMap, freqItemsets[c].keySet());
            final double avgRelFreqDiff = getAvgRelFreqDiff(freqItemsets[c], sampleFreqItemsetToSup);
            System.out.println("\t\tARFD: " + avgRelFreqDiff);
            stats.put("NumFreqItemsets_" + c, String.valueOf(freqItemsets[c].size()));
            stats.put("AvgRelFreqDiff_" + c, String.valueOf(avgRelFreqDiff));
        }
        fieldValueMap.putAll(stats);
    }

    /**
     * Gets the sample's itemset to support map such that the map only contains
     * itemsets that are frequent itemsets of the observed dataset.
     *
     * @param sample the sample
     * @param freqItemsets the set of frequent itemsets in the observed dataset
     * @param T transformer to get the id mapping for the starting graph
     * @param T2 transformer to get the id mapping for the sample
     * @return a map where each key is an itemset and the value is the itemset's
     * support in the sampled dataset such that the map only contains itemsets
     * that are frequent itemsets of the observed dataset
     */
    private static Map<Set<Integer>, Integer> getItemsetToSupMap(
            int[][] sample,
            Set<Set<Integer>> freqItemsets) {

        return freqItemsets
                .parallelStream()
                .map(freqItemset -> new Pair<Set<Integer>, Integer>(freqItemset,
                getItemsetInSampleCount(sample, freqItemset)))
                .collect(Collectors.toMap(e -> e.getValue0(), e -> e.getValue1()));
    }

    /**
     * Determines how many times the itemset is in the sample.
     *
     * @param sample the sample
     * @param itemset the itemset
     * @param T transformer to get the id mapping for the starting graph
     * @param T2 transformer to get the id mapping for the sample
     * @return number of occurrences of the itemset in the sample
     */
    private static int getItemsetInSampleCount(
            int[][] sample,
            Set<Integer> itemset) {

        int firstV = itemset.iterator().next();
        int[] tmp = sample[firstV];
        if (itemset.size() > 1) {
            for (int item : itemset) {
                tmp = computeIntersection(tmp, sample[item]);
                if (tmp.length == 0) {
                    return 0;
                }
            }
        }
        return tmp.length;
    }

    /**
     * Compute the intersection between two sorted arrays of integers.
     *
     * @param values_a first array
     * @param values_b second array
     * @return the common elements in values_a and values_b
     */
    public static int[] computeIntersection(int[] values_a, int[] values_b) {
        int i = 0;
        int j = 0;
        IntArrayList tmp_inter = new IntArrayList();

        while (true) {
            if (i >= values_a.length) {
                break;
            }
            if (j >= values_b.length) {
                break;
            }
            if (values_a[i] == values_b[j]) {
                tmp_inter.add(values_b[j]);
                j++;
                i++;
                continue;
            }
            if (values_a[i] > values_b[j]) {
                j++;
                continue;
            }
            i++;
        }
        return tmp_inter.toIntArray();
    }

    /**
     * Gets the average relative frequency/support difference.
     *
     * @param fItemsetToSup a map where each key is a frequent itemset in the
     * observed dataset and the value is the frequent itemset's support
     * @param sItemsetToSup a map where each key is an itemset and the value is
     * the itemset's support in the sampled dataset such that the map only
     * contains itemsets that are frequent itemsets of the observed dataset
     * @return the average relative frequency/support difference
     */
    private static double getAvgRelFreqDiff(
            Map<Set<Integer>, Integer> fItemsetToSup,
            Map<Set<Integer>, Integer> sItemsetToSup) {

        double sumRelFreqDiff = fItemsetToSup.entrySet().stream()
                .mapToDouble(e -> 1. * Math.abs(e.getValue() - sItemsetToSup.getOrDefault(e.getKey(), 0)) / e.getValue())
                .sum();
        return sumRelFreqDiff / fItemsetToSup.size();
    }
}
