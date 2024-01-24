package samosa.test;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import samosa.utils.CMDLineParser;
import samosa.utils.Config;
import samosa.utils.Timer;
import java.time.LocalDateTime;
import java.util.List;
import java.util.Map;
import samosa.samplers.NuDHy;
import samosa.samplers.NuDHy_Degs;
import samosa.samplers.NuDHy_BIOT;
import samosa.utils.io.Writer;

/**
 * This class runs the scalability experiment, which measures the step times of
 * the NuDHy samplers.
 */
public class Scalability {

    public static void main(String[] args) throws IOException, Exception {

        CMDLineParser.parse(args);

        final NuDHy[] samplers = {
            new NuDHy_Degs(),
            new NuDHy_BIOT()
        };

        System.out.println("Executing runtime experiment for " + Config.datasetName);
        
        final String filePath = Path.of(Config.datasetsDir, Config.datasetName).toString();

        long startTime;
        List<Map<String, String>> runStats = Lists.newArrayList();
        
        for (NuDHy sampler : samplers) {
            for (long seed = 0; seed < Config.numSamples; seed++) {
                final Map<String, String> timeStats = Maps.newHashMap();
                final Timer timer = new Timer();

                startTime = System.nanoTime();
                // move in the Markov graph for numSwaps steps
                sampler.sample(filePath, Config.numSwaps, seed, timer, false);

                timeStats.put("TotalTime", String.valueOf((System.nanoTime() - startTime)));
                timeStats.put("SetupTime", String.valueOf(timer.getSetupTime()));
                timeStats.put("MinStepTime", String.valueOf(timer.getMin()));
                timeStats.put("C10StepTime", String.valueOf(timer.getPercentile(10)));
                timeStats.put("Q1StepTime", String.valueOf(timer.getPercentile(25)));
                timeStats.put("MedianStepTime", String.valueOf(timer.getPercentile(50)));
                timeStats.put("Q3StepTime", String.valueOf(timer.getPercentile(75)));
                timeStats.put("C90StepTime", String.valueOf(timer.getPercentile(90)));
                timeStats.put("MaxStepTime", String.valueOf(timer.getMax()));
                timeStats.put("Sampler", sampler.getClass().getName());
                timeStats.put("NumSwaps", String.valueOf(Config.numSwaps));
                timeStats.put("Sample", String.valueOf(seed));
                runStats.add(timeStats);
            }
            System.out.println(Config.numSamples + " created for sample " + sampler.getClass().getName());
        }
        // save stats
        final LocalDateTime date = LocalDateTime.now();
        
        final String resultsBaseName = String.join("__", 
                "Data_" + Config.datasetName, 
                "NumSwaps_" + String.valueOf(Config.numSwaps),
                "NumSamples" + String.valueOf(Config.numSamples),
                "ScalabilityStats.tsv");
        Files.createDirectories(Path.of(Config.resultsDir, "scalability"));
        final String resultsPath = Path.of(Config.resultsDir, "scalability", resultsBaseName).toString();
        Writer.writeStatistics(date, runStats, resultsPath, true);
        System.out.println("Result written to " + resultsPath);
    }
}
