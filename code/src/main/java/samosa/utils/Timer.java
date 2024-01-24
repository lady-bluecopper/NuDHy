package samosa.utils;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * A class to log the setup time and step times of samplers.
 */
public class Timer {

    /**
     * An object to store the logged times and compute quartiles from.
     */
    private final DescriptiveStatistics times;

    /**
     * The start time.
     */
    private long start;

    /**
     * A variable to store the setup time.
     */
    private long setupTime;
    
    /**
     * A variable to store the time required to perform all the iterations.
     */
    private long walkTime;

    
    public Timer() {
        this.times = new DescriptiveStatistics();
    }

    /**
     * Starts the timer.
     */
    public void start() {
        this.start = System.nanoTime();
    }

    /**
     * Stops the timer and saves the elapsed time.
     * @return elapsed time
     */
    public long stop() {
        final long elapsed = System.nanoTime() - this.start;
        this.times.addValue(elapsed);
        return elapsed;
    }

    /**
     * Save time in setupTime variable.
     * @param time elapsed time
     */
    public void saveSetup(long time) {
        this.setupTime = time;
    }
    
    /**
     * Save time in walkTime variable.
     * @param time elapsed time
     */
    public void saveWalk(long time) {
        this.walkTime = time;
    }

    /**
     * 
     * @return elapsed time stored in setupTime variable
     */
    public long getSetupTime() {
        return this.setupTime;
    }
    
    /**
     * 
     * @return elapsed time stored in walkTime variable
     */
    public long getWalkTime() {
        return walkTime;
    }

    /**
     * 
     * @return min elapsed time stored so far; 0 if none
     */
    public double getMin() {
        double min = this.times.getMin();
        if (Double.compare(min, Double.NaN) != 0) {
            return min;
        }
        return 0.0;
    }

    /**
     * 
     * @return max elapsed time stored so far; 0 if none 
     */
    public double getMax() {
        double max = this.times.getMax();
        if (Double.compare(max, Double.NaN) != 0) {
            return max;
        }
        return 0.0;
    }

    /**
     * 
     * @param percentile percentile
     * @return desired percentile of the elapsed times stored so far; 0 if none
     */
    public double getPercentile(double percentile) {
        double time = this.times.getPercentile(percentile);
        if (Double.compare(time, Double.NaN) != 0) {
            return time;
        }
        return 0.0;
    }
}
