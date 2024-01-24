package samosa.test;

import samosa.utils.CMDLineParser;
import samosa.utils.Config;
import samosa.utils.Timer;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.util.List;
import java.util.Map;
import samosa.samplers.NuDHy;
import samosa.samplers.NuDHy_Degs;
import samosa.samplers.Sampler;
import samosa.samplers.UnpretentiousNullModel;
import samosa.utils.io.Writer;

/**
 * This class generates numSamples random hypergraphs using the sampler samplerType.
 * NuDHy samplers perform maxNumSwapsFactor * |E| swaps before returning the current
 * state in the Markov chain.
 * @author giulia
 */
public class Sampling {

    public static void main(String[] args) throws IOException, ClassNotFoundException, Exception {

        CMDLineParser.parse(args);

        String outDir = Config.resultsDir;

        String baseName = Config.datasetName.substring(0, Config.datasetName.length() - 4);
        Config.resultsDir = Path.of(outDir, baseName).toString();
        Files.createDirectories(Path.of(outDir, baseName));

        System.out.println("Generate samples for " + Config.datasetName + " using " + Config.samplerType);
        final String filePath = Path.of(Config.datasetsDir, Config.datasetName).toString();

        // initialize sampler to use
        Sampler sampler;
        if (Config.samplerType.equalsIgnoreCase("UnpretentiousNullModel")) {
            sampler = new UnpretentiousNullModel(Config.headTailDisjoint);
        } else {
            NuDHy_Degs mold_A = new NuDHy_Degs();
            mold_A.parseInputFile(filePath);
            // number of steps in the Markov graph to perform
            Config.numSwaps = (int) Config.maxNumSwapsFactor * (mold_A.size_D0 + mold_A.size_D1);
            System.out.println(baseName + ": " + Config.numSwaps);
            sampler = NuDHy.getSampler(Config.samplerType);
        }
        System.out.println("#### SAMPLER: " + sampler.getClass().getName() + " ####");
        // we run the sampler to generate the random graphs
        List<Map<String, String>> runStats = runSampler(filePath, sampler);
        // save stats
        final LocalDateTime date = LocalDateTime.now();
        String resultsBaseName = String.join("__",
                "Data_" + Config.datasetName,
                "NumSwaps_" + String.valueOf(Config.numSwaps),
                "Sampler_" + Config.samplerType,
                "NumSamples_" + String.valueOf(Config.numSamples),
                "SamplingStats.tsv");
        final String resultsPath = Path.of(outDir, "samples", resultsBaseName).toString();
        Writer.writeStatistics(date, runStats, resultsPath, false);
        System.out.println("Result written to " + resultsPath);
    }

    /**
     * Generates numSamples random hypergraphs starting from the hypergraph
     * reads from filePath.
     * @param filePath path to original hypergraph
     * @param sampler a sampler
     * @return random hypergraph sampled by sampler starting from
     * the original hypergraph
     * @throws IOException
     */
    private static List<Map<String, String>> runSampler(String filePath,
            Sampler sampler) throws IOException, Exception {

        List<Map<String, String>> allStats = Lists.newArrayList();
        Map<String, String> timeStats;
        for (long i = 0; i < Config.numSamples; i++) {
            timeStats = sample(filePath, sampler, i);
            allStats.add(timeStats);
        }
        return allStats;
    }

    /**
     * Generates a random sample and writes it to disk.
     *
     * @param dirHyperGraph original hypergraph
     * @param sampler a sampler
     * @param seed seed for reproducibility
     * @return statistics on the sampling process
     */
    private static Map<String, String> sample(String filePath,
            Sampler sampler,
            long seed) throws IOException, Exception {
        final Timer timer = new Timer();
        final long start = System.nanoTime();
        System.out.println("\t\tNum Swaps: " + Config.numSwaps);
        final Map<String, String> timeStats = Maps.newHashMap();
        System.out.println("\t\tGetting sample");
        sampler.sample(filePath, Config.numSwaps, seed, timer, Config.store);
        timer.saveWalk(System.nanoTime() - start);
        timeStats.put("MinStepTime", String.valueOf(timer.getMin()));
        timeStats.put("C10StepTime", String.valueOf(timer.getPercentile(10)));
        timeStats.put("Q1StepTime", String.valueOf(timer.getPercentile(25)));
        timeStats.put("MedianStepTime", String.valueOf(timer.getPercentile(50)));
        timeStats.put("Q3StepTime", String.valueOf(timer.getPercentile(75)));
        timeStats.put("C90StepTime", String.valueOf(timer.getPercentile(90)));
        timeStats.put("MaxStepTime", String.valueOf(timer.getMax()));
        timeStats.put("TotalTime", String.valueOf(timer.getWalkTime()));
        return timeStats;
    }

}
