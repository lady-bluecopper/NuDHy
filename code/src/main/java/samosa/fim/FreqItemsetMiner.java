package samosa.fim;

/*
 * Copyright (C) 2022 Alexander Lee and Matteo Riondato
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
import java.io.IOException;
import java.util.Map;
import java.util.Set;

/**
 * A wrapper class for mining frequent itemsets.
 */
public class FreqItemsetMiner {

    /**
     * Mines frequent itemsets and saves them to disk.
     *
     * @param transactions database
     * @param minFreq the minimum frequency threshold
     * @param freqItemsetsPath the path to save the frequent itemsets
     * @throws java.io.IOException
     */
    public static void mine(int[][] transactions,
            double minFreq,
            String freqItemsetsPath) throws IOException {
        final AlgoNegFIN algo = new AlgoNegFIN(transactions, minFreq);
        algo.runAlgorithm(transactions, minFreq, freqItemsetsPath);
    }

    /**
     * Mines frequent itemsets and returns them as a map where each key is a
     * frequent itemset and the value is the frequent itemset's support.
     *
     * @param transactions database
     * @param minFreq the minimum frequency threshold
     * @return a map where each key is a frequent itemset and the value is the
     * frequent itemset's support
     * @throws java.io.IOException
     */
    public static Map<Set<Integer>, Integer> mine(int[][] transactions,
            double minFreq) throws IOException {

        final AlgoNegFINMod algo = new AlgoNegFINMod(transactions, minFreq);
        return algo.runAlgorithm(transactions, minFreq);
    }
}
