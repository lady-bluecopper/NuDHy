package samosa.test;

import samosa.utils.CMDLineParser;
import samosa.utils.Config;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import samosa.samplers.NuDHy;
import samosa.samplers.NuDHy_Degs;

/**
 * Starting from the observed hypergraph, move in the Markov graph using the sampler
 * samplerType and store on disk the current state after k steps up to max_num_swaps steps.
 * The procedure is executed numSamples times.
 * The output is used in he Convergence experiment.
 */
public class SamplingForConvergence {

    public static void main(String[] args) throws IOException, ClassNotFoundException, Exception {

        CMDLineParser.parse(args);

        String baseName = Config.datasetName.substring(0, Config.datasetName.length() - 4);
        Files.createDirectories(Path.of(Config.resultsDir, baseName));
        Config.resultsDir = Path.of(Config.resultsDir, baseName).toString();

        final String filePath = Path.of(Config.datasetsDir, Config.datasetName).toString();

        NuDHy_Degs mold_A = new NuDHy_Degs();
        mold_A.parseInputFile(filePath);
        final int max_num_swaps = (int) Config.maxNumSwapsFactor * (mold_A.size_D0 + mold_A.size_D1);

        // sampler
        final NuDHy sampler = NuDHy.getSampler(Config.samplerType);
        System.out.println("#### SAMPLER: " + sampler.getClass().getName() + " ####");
        // we run the sampler to generate the random graphs
        for (int i = 0; i < Config.numSamples; i++) {
            sampler.runMC_multiFlush(max_num_swaps, i, filePath);
        }
    }

}
