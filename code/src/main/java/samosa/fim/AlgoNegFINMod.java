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
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import samosa.utils.MemoryLogger;

/**
 * A variant of SPMF's {@link AlgoNegFIN} that returns the collection of
 * frequent itemsets, instead of writing them to file.
 */
public class AlgoNegFINMod extends AlgoNegFIN {

    private Map<Set<Integer>, Integer> freqItemsetToSup; // frequent itemsets to their support
    
    public AlgoNegFINMod(int[][] transactions, double minSupport) {
        super(transactions, minSupport);
    }

    /**
     * This method adds an itemset to freqItemsetsToSup + all itemsets that can
     * be made using its node list. We keep the method name even if it does not
     * do what it says in order to reduce duplicate code.
     *
     * @param curNode the current node
     * @param sameCount the same count
     */
    @Override
    public void writeItemsetsToFile(SetEnumerationTreeNode curNode, int sameCount) {

        // initialize a frequent itemset
        Set<Integer> freqItemset = Sets.newHashSet();
        Set<Integer> otherFreqItemset;

        outputCount++;
        // append items from the itemset to the StringBuilder
        for (int i = 0; i < itemsetLen; i++) {
            freqItemset.add(item[itemset[i]].index);
        }
        // add the frequent itemset and its support to the map
        freqItemsetToSup.put(freqItemset, curNode.count);

        // === Add all combination that can be made using the node list of
        // this itemset
        if (sameCount > 0) {
            // generate all subsets of the node list except the empty set
            for (long i = 1, max = 1 << sameCount; i < max; i++) {
                otherFreqItemset = Sets.newHashSet();
                for (int k = 0; k < itemsetLen; k++) {
                    otherFreqItemset.add(item[itemset[k]].index);
                }

                // we create a new subset
                for (int j = 0; j < sameCount; j++) {
                    // check if the j bit is set to 1
                    int isSet = (int) i & (1 << j);
                    if (isSet > 0) {
                        // if yes, add it to the set
                        otherFreqItemset.add(item[sameItems[j]].index);
                        // newSet.add(item[sameItems[j]].index);
                    }
                }
                freqItemsetToSup.put(otherFreqItemset, curNode.count);
                outputCount++;
            }
        }
    }

    /**
     * Run the algorithm.
     *
     * @param transactions database
     * @param minsup the minsup threshold
     * @return a map where each key is a frequent itemset and its value is its
     * support
     * @throws IOException if error while reading file
     */
    public Map<Set<Integer>, Integer> runAlgorithm(int[][] transactions, 
            double minsup)
            throws IOException {
        
        freqItemsetToSup = Maps.newHashMap();

        MemoryLogger.getInstance().reset();

        // record the start time
        startTimestamp = System.currentTimeMillis();

        itemsetLen = 0;
        itemset = new int[numOfFItem];

        // Lines 2 to 6 of algorithm 3 in the paper
        // Build BMC-tree
        constructBMCTree(transactions); 
        // Lines 12 to 19 of algorithm 3 in the paper
        // Initialize tree
        initializeSetEnumerationTree();
        sameItems = new int[numOfFItem];

        // Recursively constructing_frequent_itemset_tree the tree
        SetEnumerationTreeNode curNode = nlRoot.firstChild;
        nlRoot.firstChild = null;
        SetEnumerationTreeNode next = null;
        while (curNode != null) {
            next = curNode.next;
            // call the recursive "constructing_frequent_itemset_tree" method
            constructFreqItemsetTree(curNode, 1, 0);
            curNode.next = null;
            curNode = next;
        }
        MemoryLogger.getInstance().checkMemory();
        // record the end time
        endTimestamp = System.currentTimeMillis();
        return freqItemsetToSup;
    }
    
}
