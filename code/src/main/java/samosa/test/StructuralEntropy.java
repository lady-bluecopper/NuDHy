package samosa.test;

import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import org.apache.commons.io.FilenameUtils;
import samosa.samplers.NuDHy;
import samosa.samplers.NuDHy_Degs;
import samosa.utils.CMDLineParser;
import samosa.utils.Config;
import samosa.utils.io.Reader;
import samosa.utils.io.Writer;
import org.apache.commons.math3.util.Combinations;

/**
 * This class find the probability that a node group of size lower than 4 that appears
 * in some hyperedge in the observed hypergraph, also appears in the random
 * samples generated by NuDHy.
 * 
 * @author giulia
 */
public class StructuralEntropy {
    
    public static void main(String[] args) throws IOException, Exception {

        CMDLineParser.parse(args);

        // OBSERVED HYPERGRAPH
        final String filePath = Path.of(Config.datasetsDir, Config.datasetName).toString();
        // CREATE DIRECTORY FOR TMP FILES
        Files.createDirectories(Path.of(Config.resultsDir, "entropy", "tmp"));
        // INITIALIZE EXECUTOR
        ExecutorService executor = Executors.newFixedThreadPool(Config.numThreads);
        ObjectArrayList<String> fileNames = new ObjectArrayList();
        // READ SAMPLES FROM ZIP FILE
        final String dataName = Config.datasetName.split("\\.")[0];
        final String zipPath = Path.of(Config.resultsDir, "samples", Config.samplerType, dataName + ".zip").toString();
        System.out.println("Reading samples from " + zipPath);
        ZipFile zip = new ZipFile(zipPath);
        for (Enumeration e = zip.entries(); e.hasMoreElements();) {
            ZipEntry entry = (ZipEntry) e.nextElement();
            String fileName = entry.getName();
            if (FilenameUtils.getExtension(fileName).equals("tsv")
                    && fileName.contains(dataName)
                    && !fileName.startsWith(".")
                    && !fileName.startsWith("__MACOSX")) {
                // parse file name to extract relevant fields
                Map<String, String> fieldValueMap = Reader.parseOutputFileName(fileName);
                String resultsName = dataName + "__"
                        + "randomSeed_" + fieldValueMap.get("randomSeed") + "__"
                        + "algorithm_" + fieldValueMap.get("algorithm") + "__";
                if (fieldValueMap.containsKey("iterations")) {
                    resultsName += "iterations_" + fieldValueMap.get("iterations") + "__"; 
                }
                resultsName += "occurrences";
                String resultsPath = Path.of(Config.resultsDir, "entropy", "tmp", resultsName).toString();
                fileNames.add(resultsPath);

                executor.submit(() -> {
                    try {
                        // find and save group occurrences for the sample read
                        // from fileName
                        ObjectOpenHashSet<String>[][] out;
                            out = populateOccSubsetsInSample(zipPath, fileName, filePath);
                        Writer.writeOccurrencesOfSample(out[0], resultsPath, "heads");
                        Writer.writeOccurrencesOfSample(out[1], resultsPath, "tails");
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }
                });
            }
        }
        executor.shutdown();
        try {
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        System.out.println("Output processed. " + System.currentTimeMillis());
        // GET ALL POSSIBLE NODE GROUPS IN THE OBSERVED HYPERGRAPH
        ObjectOpenHashSet<String>[] all_subsets = getSubsetsInObservedHypergraph(filePath);
        System.out.println("Computed node groups in observed hypergraph. " + System.currentTimeMillis());
        // INITALIZE RESULT OBJECT: node group -> occurs in sample? 
        Object2ObjectOpenHashMap<String, Object2DoubleOpenHashMap<String>> occurrences;
        // RESULT PATH
        String result_path = Path.of(Config.resultsDir, "entropy", Config.datasetName
                + "__numSamples_" + Config.numSamples 
                + "__algorithm_" + Config.samplerType).toString();
        System.out.println(result_path);
        String[] tasks = {"heads", "tails"};
        System.out.println("Start processing output. " + System.currentTimeMillis());
        // PROCESS OUTPUT OF SAMPLES
        for (int i = 0; i < tasks.length; i++) {
            String task = tasks[i];
            // INITALIZE RESULT OBJECT: node group -> occurs in sample? 
            occurrences = new Object2ObjectOpenHashMap();
            for (String outFilePath : fileNames) {
                Map<String, String> fieldValueMap = Reader.parseOutputFileName(outFilePath);
                String sampler = fieldValueMap.get("algorithm");
                if (!occurrences.containsKey(sampler)) {
                    Object2DoubleOpenHashMap map = new Object2DoubleOpenHashMap(2 + (int) (all_subsets[i].size() / Object2DoubleOpenHashMap.DEFAULT_LOAD_FACTOR));
                    map.defaultReturnValue(0.);
                    occurrences.put(sampler, map);
                }
                ObjectOpenHashSet<String> out = Reader.readOccurrencesOfSample(outFilePath, task);
                // UPDATE OCCURRENCE MAP
                for (String subset : out) {
                    occurrences.get(sampler).addTo(subset, 1);
                }
            }
            // SAVE DATA
            boolean isFirst = true;
            for (String sampler : occurrences.keySet()) {
                Object2DoubleOpenHashMap<String> oS = occurrences.get(sampler);
                Writer.writeStructuralProbabilities(result_path, task, all_subsets[i], oS, sampler, !isFirst);
                // set to False so that next results are appended to the same output file
                isFirst = false;
            }
            System.out.println("Computed " + task + " probabilities. " + System.currentTimeMillis());
            occurrences = null;
            System.gc();
        }
    }

    /**
     *
     * @param zipPath path to the zip file including all the samples
     * @param fileName file name of the sample to read from the zip file
     * @param dataPath path to the observed hypergraph
     * @return for heads and tails, for each size up to 4, the subsets of nodes
     * of that size in the observed hypergraph that appear in the sample read
     * from the zip file.
     * @throws IOException
     * @throws Exception
     */
    private static ObjectOpenHashSet<String>[][] populateOccSubsetsInSample(
            String zipPath,
            String fileName,
            String dataPath) throws IOException, Exception {

        ObjectOpenHashSet<String>[][] out = new ObjectOpenHashSet[2][3];
        System.out.println("Finding occurrences for: " + fileName);
        ZipFile zip = new ZipFile(zipPath);
        for (Enumeration e = zip.entries(); e.hasMoreElements();) {
            ZipEntry entry = (ZipEntry) e.nextElement();
            String fileName2 = entry.getName();
            if (!fileName.equals(fileName2)) {
                continue;
            }
            // READ SAMPLE
            BufferedReader reader = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
            BufferedReader reader2 = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
            NuDHy sample = new NuDHy_Degs();
            sample.parseInputFile(reader, reader2);
            // INITALIZE RESULT OBJECT: node groups that occur in sample
            ObjectOpenHashSet<String>[] soccurrencesInHeads = new ObjectOpenHashSet[3];
            ObjectOpenHashSet<String>[] soccurrencesInTails = new ObjectOpenHashSet[3];
            // INITALIZE TMP OBJECT: node groups that do not occur in sample 
            ObjectOpenHashSet<String>[] smissesInHeads = new ObjectOpenHashSet[3];
            ObjectOpenHashSet<String>[] smissesInTails = new ObjectOpenHashSet[3];
            for (int i = 0; i < 3; i++) {
                soccurrencesInHeads[i] = new ObjectOpenHashSet();
                soccurrencesInTails[i] = new ObjectOpenHashSet();
                smissesInHeads[i] = new ObjectOpenHashSet();
                smissesInTails[i] = new ObjectOpenHashSet();
            }
            // READ OBSERVED HYPERGRAPH
            BufferedReader br = new BufferedReader(new FileReader(dataPath));
            String line;
            String[] head_and_tail;
            String[] head;
            String[] tail;
            // iterate through the hyperedges of the observed hypergraph
            while ((line = br.readLine()) != null) {
                head_and_tail = line.split("\t");
                head = head_and_tail[0].split(",");
                tail = head_and_tail[1].split(",");
                // Search for all possible subgroups of head
                findOccSubsetsOfHyperedge(sample, head, soccurrencesInHeads, smissesInHeads, true);
                // Search for all possible subgroups of tail
                findOccSubsetsOfHyperedge(sample, tail, soccurrencesInTails, smissesInTails, false);
            }
            br.close();
            out[0] = soccurrencesInHeads;
            out[1] = soccurrencesInTails;
            return out;
        }
        return out;
    }

    /**
     *
     * @param dataPath path to the observed hypergraph
     * @return for heads and tails, for each size up to 4, the subsets of nodes
     * of that size in the observed hypergraph.
     * @throws IOException
     * @throws Exception
     */
    protected static ObjectOpenHashSet<String>[] getSubsetsInObservedHypergraph(String dataPath) throws IOException, Exception {

        ObjectOpenHashSet<String>[] out = new ObjectOpenHashSet[2];
        out[0] = new ObjectOpenHashSet();
        out[1] = new ObjectOpenHashSet();
        // READ OBSERVED HYPERGRAPH
        BufferedReader br = new BufferedReader(new FileReader(dataPath));
        String line;
        String[] head_and_tail;
        while ((line = br.readLine()) != null) {
            head_and_tail = line.split("\t");
            for (int id = 0; id < head_and_tail.length; id++) {
                String[] itemset = head_and_tail[id].split(",");
                int maxSize = Math.min(itemset.length, 4);
                for (int size = 2; size <= maxSize; size++) {
                    Combinations iterable = new Combinations(itemset.length, size);
                    String[] combo = new String[size];
                    for (int[] combination : iterable) {
                        for (int i = 0; i < size; i++) {
                            combo[i] = itemset[combination[i]];
                        }
                        Arrays.sort(combo);
                        String tmp_combo_string = Arrays.toString(combo).replace(" ", "");
                        out[id].add(tmp_combo_string.substring(1, tmp_combo_string.length() - 1));
                    }
                }
            }
        }
        br.close();
        return out;
    }

    /**
     * Populates the set of node groups in the observed hypergraph that appear
     * in some hyperedge of the sampled hypergraph.
     *
     * @param sample random directed hypergraph
     * @param itemset set of nodes under examination
     * @param occurrences subsets of itemset of size up to 4 that appear in some
     * hyperedge (head/tail) of sample
     * @param misses subsets of itemset of size up to 4 that do not appear in
     * any hyperedge (head/tail) of sample
     * @param isHead whether the subset should be searched in the heads or in
     * the tails of the hyperedges in the sample
     */
    protected static void findOccSubsetsOfHyperedge(
            NuDHy sample,
            String[] itemset,
            ObjectOpenHashSet<String>[] occurrences,
            ObjectOpenHashSet<String>[] misses,
            boolean isHead) {

        if (itemset.length < 2) {
            return;
        }
        String[] S_with_outer_ids__SIZE_2 = new String[2];
        String[] S_with_outer_ids__SIZE_2_sorted = new String[2];
        String[] S_with_outer_ids__SIZE_3 = new String[3];
        String[] S_with_outer_ids__SIZE_3_sorted = new String[3];
        String[] S_with_outer_ids__SIZE_4 = new String[4];
        String[] S_with_outer_ids__SIZE_4_sorted = new String[4];
        String subset;
        int subset_size_index;// subset_size_index ;= current subset size - 2 
        for (int i = 0; i < itemset.length; i++) {
            S_with_outer_ids__SIZE_2[0] = itemset[i];
            S_with_outer_ids__SIZE_3[0] = itemset[i];
            S_with_outer_ids__SIZE_4[0] = itemset[i];
            for (int j = i + 1; j < itemset.length; j++) {
                subset_size_index = 0;
                S_with_outer_ids__SIZE_2[1] = itemset[j];
                S_with_outer_ids__SIZE_3[1] = itemset[j];
                S_with_outer_ids__SIZE_4[1] = itemset[j];

                System.arraycopy(S_with_outer_ids__SIZE_2, 0, S_with_outer_ids__SIZE_2_sorted, 0, 2);
                Arrays.sort(S_with_outer_ids__SIZE_2_sorted);
                subset = Arrays.toString(S_with_outer_ids__SIZE_2_sorted);
                if (!occurrences[subset_size_index].contains(subset) && !misses[subset_size_index].contains(subset)) {
                    if (!sample.isSinAtLeastOneHyperedge(S_with_outer_ids__SIZE_2_sorted, isHead)) {
                        // this group does not appear in the sample and thus
                        // no extension can appear in the sample
                        misses[subset_size_index].add(subset);
                        continue;
                    }
                    occurrences[subset_size_index].add(subset);
                }
                if (itemset.length < 3) {
                    // we cannot expand it anymore
                    continue;
                }

                for (int k = j + 1; k < itemset.length; k++) {
                    subset_size_index = 1;
                    S_with_outer_ids__SIZE_3[2] = itemset[k];
                    S_with_outer_ids__SIZE_4[2] = itemset[k];

                    System.arraycopy(S_with_outer_ids__SIZE_3, 0, S_with_outer_ids__SIZE_3_sorted, 0, 3);
                    Arrays.sort(S_with_outer_ids__SIZE_3_sorted);
                    subset = Arrays.toString(S_with_outer_ids__SIZE_3_sorted);
                    if (!occurrences[subset_size_index].contains(subset) && !misses[subset_size_index].contains(subset)) {
                        if (!sample.isSinAtLeastOneHyperedge(S_with_outer_ids__SIZE_3_sorted, isHead)) {
                            // this group does not appear in the sample and thus
                            // no extension can appear in the sample
                            misses[subset_size_index].add(subset);
                            continue;
                        }
                        occurrences[subset_size_index].add(subset);
                    }
                    if (itemset.length < 4) {
                        // we cannot expand it anymore
                        return;
                    }
                    for (int l = k + 1; l < itemset.length; l++) {

                        subset_size_index = 2;
                        S_with_outer_ids__SIZE_4[3] = itemset[l];

                        System.arraycopy(S_with_outer_ids__SIZE_4, 0, S_with_outer_ids__SIZE_4_sorted, 0, 4);
                        Arrays.sort(S_with_outer_ids__SIZE_4_sorted);
                        subset = Arrays.toString(S_with_outer_ids__SIZE_4_sorted);
                        if (!occurrences[subset_size_index].contains(subset) && !misses[subset_size_index].contains(subset)) {
                            if (sample.isSinAtLeastOneHyperedge(S_with_outer_ids__SIZE_4_sorted, isHead)) {
                                occurrences[subset_size_index].add(subset);
                            } else {
                                misses[subset_size_index].add(subset);
                            }
                        }
                    }
                }
            }
        }
    }

}