package samosa.utils;

/*
 * Copyright (C) 2022 Giulia Preti
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
/**
 * Variables to store the hyper-parameters and the experimental settings.
 */
public class Config {
    // directory where the datasets are placed
    public static String datasetsDir = "../data/";
    // path to the dataset
    public static String datasetName = "eco01100.tsv";
    // directory where the output is saved
    public static String resultsDir = "../out/";
    // type of sampler
    public static String samplerType = "UnpretentiousNullModel";
    // number of iterations to perform before returning the sample
    public static int numSwaps = 1660900;
    // multiplicator to obtain the max number of swaps to perform
    public static double maxNumSwapsFactor = 50;
    // whether heads and tails must be disjoint when generating 
    // a random hypergraph using UnpretentiousNullModel
    public static boolean headTailDisjoint = false;
    // number of random sample to generate/consider
    public static int numSamples = 33;
    // minimum frequency for an itemset to be frequent
    public static String minFreqs = "1e-3,1";
    // top-k frequent itemsets
    public static int k = 20;
    // min size of a frequent itemset to be returned
    public static int size = 3;
    // seed for reproducibility
    public static long seed = 0;
    // parallel execution
    public static int numThreads = 4;
    // store samples on disk
    public static boolean store = true;
}
