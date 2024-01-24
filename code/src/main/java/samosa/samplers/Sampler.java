package samosa.samplers;

import samosa.utils.Timer;

/**
 *
 * @author giulia
 */
public interface Sampler {
    
    /**
     * Starting from the original hypergraph, generate a random hypergraph
     * and store it on disk, if needed.
     * @param fileName path to original hypergraph
     * @param numSwaps number of steps in the Markov graph (if a MC sampler is used)
     * @param seed seed for reproducibility
     * @param timer object to store the running times
     * @param save whether the sample should be saved on disk
     * @throws java.lang.Exception
     */
    public abstract void sample(String fileName, int numSwaps, long seed, Timer timer, boolean save) throws Exception;
    
}
