package samosa.fim;


import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import samosa.utils.MemoryLogger;

/*
 ** The implementation of the "negFIN algorithm", the algorithm presented in:
 * "Nader Aryabarzan, Behrouz Minaei-Bidgoli, and Mohammad Teshnehlab. (2018). 
 * negFIN: An efficient algorithm for fast mining frequent itemsets. 
 * Expert System with Applications, 105, 129–143"
 *
 * This file is part of the SPMF DATA MINING SOFTWARE
 * (http://www.philippe-fournier-viger.com/spmf).
 *
 * SPMF is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * SPMF is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * SPMF. If not, see <http://www.gnu.org/licenses/>.
 */
/**
 * This implementation was obtained by converting the C++ code of the negFIN
 * algorithm to Java. The C++ code of this algorithm was provided by Nader
 * Aryabarzan, available on GitHub via https://github.com/aryabarzan/negFIN/.
 *
 * <p>
 * Both the C++/Java code of the negFIN algorithms are respectively based on the
 * C++/Java code of the "FIN algorithm", the algorithm which is presented in:
 * "Z. H. Deng and S. L. Lv. (2014). Fast mining frequent itemsets using
 * Nodesets. Expert System with Applications, 41, 4505–4512"
 *
 * @author Nader Aryabarzan (Copyright 2018) @Email aryabarzan@aut.ac.ir or
 * aryabarzan@gmail.com
 */
public class AlgoNegFIN {

    // the start time and end time of the last algorithm execution
    public long startTimestamp;
    public long endTimestamp;

    // Tree stuff
    public BMCTreeNode bmcTreeRoot; // The root of BMC_tree
    public SetEnumerationTreeNode nlRoot; // The root of set enumeration tree.

    private final int numOfTrans; // // Number of transactions
    public int numOfFItem; // Number of items
    public int outputCount = 0; // number of itemsets found
    public int minSupport; // minimum count
    public Item[] item; // list of items sorted by count
    public int[] itemset; // the current itemset
    public int itemsetLen = 0; // the size of the current itemset
    public int[] sameItems;
    public Map<Integer, List<BMCTreeNode>> mapItemNodeset; // nodessets of 1-itemsets

    BufferedWriter writer = null; // object to write the output file

    /**
     * Comparator to sort items by decreasing order of frequency
     */
    static Comparator<Item> comp = (Item a, Item b) -> ((Item) b).num - ((Item) a).num;
    
    public AlgoNegFIN(int[][] transactions, 
            double minSupport) {
    
        this.numOfTrans = transactions.length;
        // (1) Scan the database and count the count of each item.
        // The count of items is stored in map where
        // key = item value = count count
        Map<Integer, Integer> mapItemCount = Maps.newHashMap();
        for (int[] transaction : transactions) {
            // for each item in the transaction
            for (int it : transaction) {
                // increase the count count of the item by 1
                mapItemCount.put(it, mapItemCount.getOrDefault(it, 0) + 1);
            }
        }
        this.numOfFItem = mapItemCount.size();
        this.minSupport = (int) Math.ceil(minSupport * numOfTrans);
        this.item = mapItemCount.entrySet()
                .stream()
                .filter(e -> e.getValue() >= this.minSupport)
                .map(e -> new Item(e.getKey(), e.getValue()))
                .toArray(Item[]::new);
        this.numOfFItem = this.item.length;
        Arrays.sort(this.item, comp);
    }
    
    /**
     * Build the tree
     *
     * @param transactions database of transactions
     * @throws IOException if an exception while reading/writing to file
     */
    public void constructBMCTree(int[][] transactions) throws IOException {

        bmcTreeRoot = new BMCTreeNode();
        int bmcTreeNodeCount = 0;
        bmcTreeRoot.label = -1;
        bmcTreeRoot.bitmapCode = new MyBitVector(numOfFItem);

        // we will use a buffer to store each transaction that is read.
        Item[] transaction = new Item[numOfFItem];
        // for each transaction
        for (int[] transaction1 : transactions) {
            // for each item in the transaction
            int tLen = 0; // tLen
            for (int itemX : transaction1) {
                // add each item from the transaction except infrequent item
                for (int j = 0; j < numOfFItem; j++) {
                    // if the item is frequent, we add it
                    if (itemX == item[j].index) {
                        transaction[tLen] = new Item(itemX, 0 - j);
                        tLen++;
                        break;
                    }
                }
            }
            // sort the transaction
            Arrays.sort(transaction, 0, tLen, comp);
            int curPos = 0;
            BMCTreeNode curRoot = (bmcTreeRoot);
            BMCTreeNode rightSibling = null;
            while (curPos != tLen) {
                BMCTreeNode child = curRoot.firstChild;
                while (child != null) {
                    if (child.label == 0 - transaction[curPos].num) {
                        curPos++;
                        child.count++;
                        curRoot = child;
                        break;
                    }
                    if (child.rightSibling == null) {
                        rightSibling = child;
                        child = null;
                        break;
                    }
                    child = child.rightSibling;
                }
                if (child == null) {
                    break;
                }
            }
            for (int j = curPos; j < tLen; j++) {
                BMCTreeNode bmcTreeNode = new BMCTreeNode();
                bmcTreeNode.label = 0 - transaction[j].num;
                if (rightSibling != null) {
                    rightSibling.rightSibling = bmcTreeNode;
                    rightSibling = null;
                } else {
                    curRoot.firstChild = bmcTreeNode;
                }
                bmcTreeNode.rightSibling = null;
                bmcTreeNode.firstChild = null;
                bmcTreeNode.father = curRoot;
                bmcTreeNode.count = 1;
                curRoot = bmcTreeNode;
                bmcTreeNodeCount++;
            }
        }
        BMCTreeNode root = bmcTreeRoot.firstChild;
        mapItemNodeset = Maps.newHashMap();
        while (root != null) {
            root.bitmapCode = (MyBitVector) root.father.bitmapCode.clone();
            root.bitmapCode.set(root.label); // bitIndex=numOfFItem - 1 - root.label
            List<BMCTreeNode> nodeset = mapItemNodeset.get(root.label);
            if (nodeset == null) {
                nodeset = Lists.newArrayList();
                mapItemNodeset.put(root.label, nodeset);
            }
            nodeset.add(root);

            if (root.firstChild != null) {
                root = root.firstChild;
            } else {
                if (root.rightSibling != null) {
                    root = root.rightSibling;
                } else {
                    root = root.father;
                    while (root != null) {
                        if (root.rightSibling != null) {
                            root = root.rightSibling;
                            break;
                        }
                        root = root.father;
                    }
                }
            }
        }
    }

    /**
     * Initialize the tree
     */
    public void initializeSetEnumerationTree() {
        
        nlRoot = new SetEnumerationTreeNode();
        nlRoot.label = numOfFItem;
        nlRoot.firstChild = null;
        nlRoot.next = null;
        SetEnumerationTreeNode lastChild = null;
        for (int t = numOfFItem - 1; t >= 0; t--) {
            SetEnumerationTreeNode nlNode = new SetEnumerationTreeNode();
            nlNode.label = t;
            nlNode.count = 0;
            nlNode.nodeset = mapItemNodeset.get(t);
            nlNode.firstChild = null;
            nlNode.next = null;
            nlNode.count = item[t].num;
            if (nlRoot.firstChild == null) {
                nlRoot.firstChild = nlNode;
                lastChild = nlNode;
            } else {
                lastChild.next = nlNode;
                lastChild = nlNode;
            }
        }
    }

    /**
     * Recursively constructing_frequent_itemset_tree the tree to find frequent
     * itemsets
     *
     * @param curNode
     * @param level
     * @param sameCount
     * @throws IOException if error while writing itemsets to file
     */
    public void constructFreqItemsetTree(
            SetEnumerationTreeNode curNode, 
            int level, int sameCount) throws IOException {

        MemoryLogger.getInstance().checkMemory();

        SetEnumerationTreeNode sibling = curNode.next;
        SetEnumerationTreeNode lastChild = null;
        while (sibling != null) {
            SetEnumerationTreeNode child = new SetEnumerationTreeNode();

            child.nodeset = Lists.newArrayList();
            int countNegNodeset = 0;
            if (level == 1) {
                for (int i = 0; i < curNode.nodeset.size(); i++) {
                    BMCTreeNode ni = curNode.nodeset.get(i);
                    if (!ni.bitmapCode.isSet(sibling.label)) {
                        child.nodeset.add(ni);
                        countNegNodeset += ni.count;
                    }
                }
            } else {
                for (int j = 0; j < sibling.nodeset.size(); j++) {
                    BMCTreeNode nj = sibling.nodeset.get(j);
                    if (nj.bitmapCode.isSet(curNode.label)) {
                        child.nodeset.add(nj);
                        countNegNodeset += nj.count;
                    }
                }
            }
            child.count = curNode.count - countNegNodeset;

            if (child.count >= minSupport) {
                if (curNode.count == child.count) {
                    sameItems[sameCount++] = sibling.label;
                } else {
                    child.label = sibling.label;
                    child.firstChild = null;
                    child.next = null;
                    if (curNode.firstChild == null) {
                        curNode.firstChild = lastChild = child;
                    } else {
                        lastChild.next = child;
                        lastChild = child;
                    }
                }
            } else {
                child.nodeset = null;
            }
            sibling = sibling.next;
        }
        //        resultCount += Math.pow(2.0, sameCount);
        //        nlLenSum += Math.pow(2.0, sameCount) * curNode.nodeset.size();

        itemset[itemsetLen++] = curNode.label;

        // ============= Write itemset(s) to file ===========
        writeItemsetsToFile(curNode, sameCount);
        // ======== end of write to file

        SetEnumerationTreeNode child = curNode.firstChild;
        curNode.firstChild = null;

        SetEnumerationTreeNode next;
        while (child != null) {
            next = child.next;
            constructFreqItemsetTree(child, level + 1, sameCount);
            child.next = null;
            child = next;
        }
        itemsetLen--;
    }

    /**
     * This method write an itemset to file + all itemsets that can be made
     * using its node list.
     *
     * @param curNode the current node
     * @param sameCount the same count
     * @throws IOException exception if error reading/writting to file
     */
    public void writeItemsetsToFile(SetEnumerationTreeNode curNode, int sameCount) throws IOException {

        // create a stringuffer
        StringBuilder buffer = new StringBuilder();

        outputCount++;
        // append items from the itemset to the StringBuilder
        for (int i = 0; i < itemsetLen; i++) {
            buffer.append(item[itemset[i]].index);
            buffer.append(' ');
        }
        // append the count of the itemset
        buffer.append("#SUP: ");
        buffer.append(curNode.count);
        buffer.append("\n");

        // === Write all combination that can be made using the node list of
        // this itemset
        if (sameCount > 0) {
            // generate all subsets of the node list except the empty set
            for (long i = 1, max = 1 << sameCount; i < max; i++) {
                for (int k = 0; k < itemsetLen; k++) {
                    buffer.append(item[itemset[k]].index);
                    buffer.append(' ');
                }

                // we create a new subset
                for (int j = 0; j < sameCount; j++) {
                    // check if the j bit is set to 1
                    int isSet = (int) i & (1 << j);
                    if (isSet > 0) {
                        // if yes, add it to the set
                        buffer.append(item[sameItems[j]].index);
                        buffer.append(' ');
                        // newSet.add(item[sameItems[j]].index);
                    }
                }
                buffer.append("#SUP: ");
                buffer.append(curNode.count);
                buffer.append("\n");
                outputCount++;
            }
        }
        // write the strinbuffer to file and create a new line
        // so that we are ready for writing the next itemset.
        writer.write(buffer.toString());
    }

    /**
     * Print statistics about the latest execution of the algorithm to
     * System.out.
     */
    public void printStats() {
        System.out.println("========== negFIN - STATS ============");
        System.out.println(" Minsup = " + minSupport + "\n Number of transactions: " + numOfTrans);
        System.out.println(" Number of frequent  itemsets: " + outputCount);
        System.out.println(" Total time ~: " + (endTimestamp - startTimestamp) + " ms");
        System.out.println(" Max memory:" + MemoryLogger.getInstance().getMaxMemory() + " MB");
        System.out.println("=====================================");
    }

    /**
     * Run the algorithm
     *
     * @param transactions database
     * @param minsup the minsup threshold
     * @param output the output file path
     * @throws IOException if error while reading/writting to file
     */
    public void runAlgorithm(int[][] transactions, 
            double minsup, 
            String output) throws IOException {

        MemoryLogger.getInstance().reset();
        // create object for writing the output file
        writer = new BufferedWriter(new FileWriter(output));
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
        writer.close();
        MemoryLogger.getInstance().checkMemory();
        // record the end time
        endTimestamp = System.currentTimeMillis();
    }
    
    /**
     * Gets a map where each key is a frequent itemset and the value is the
     * support for the frequent itemset.
     *
     * @param freqItemsetsPath the path for the set of frequent itemsets
     * @return a map where each key is a frequent itemset and the value is the
     * support for the frequent itemset.
     */
    public Map<Set<Integer>, Integer> getFreqItemsetToSupMap(String freqItemsetsPath) {
        
        final Map<Set<Integer>, Integer> freqItemsetToSup = Maps.newHashMap();

        try {
            final BufferedReader br = new BufferedReader(new FileReader(freqItemsetsPath));
            String line = br.readLine();
            while (line != null) {
                final String[] freqItemsetAndSup = line.split(" #SUP: ");
                final Set<Integer> freqItemset = getFreqItemset(freqItemsetAndSup);
                final int sup = getSup(freqItemsetAndSup);
                freqItemsetToSup.put(freqItemset, sup);
                line = br.readLine();
            }
            br.close();
        } catch (IOException e) {
            System.err.println("Error getting frequent itemset to support map");
            e.printStackTrace();
            System.exit(1);
        }
        return freqItemsetToSup;
    }
    
    public Set<Integer> getFreqItemset(String[] freqItemsetAndSup) {
        final String freqItemsetString = freqItemsetAndSup[0];
        final Set<Integer> freqItemset = new HashSet<>();
        for (String itemString : freqItemsetString.split(" ")) {
            final int itemInt = Integer.parseInt(itemString);
            freqItemset.add(itemInt);
        }
        return freqItemset;
    }
    
    public int getSup(String[] itemsetAndSup) {
        return Integer.parseInt(itemsetAndSup[1]);
    }

    public class Item {

        public int index;
        public int num;
        
        public Item(int index, int num) {
            this.index = index;
            this.num = num;
        }
    }

    public class SetEnumerationTreeNode {

        public int label;
        public SetEnumerationTreeNode firstChild;
        public SetEnumerationTreeNode next;
        public int count;
        List<BMCTreeNode> nodeset;
    }

    public class BMCTreeNode {

        public int label;
        public BMCTreeNode firstChild;
        public BMCTreeNode rightSibling;
        public BMCTreeNode father;
        public int count;
        MyBitVector bitmapCode;
    }
}

// This class is more efficient than the built in class BitSet
class MyBitVector {

    static long[] TWO_POWER;

    static {
        TWO_POWER = new long[64];
        TWO_POWER[0] = 1;
        for (int i = 1; i < TWO_POWER.length; i++) {
            TWO_POWER[i] = TWO_POWER[i - 1] * 2;
        }
    }

    long[] bits;

    public MyBitVector(int numOfBits) {
        bits = new long[((numOfBits - 1) / 64) + 1];
    }

    public Object clone() {
        MyBitVector result = new MyBitVector(this.bits.length * 64);
        System.arraycopy(this.bits, 0, result.bits, 0, result.bits.length);
        return result;
    }

    public void set(int bitIndex) {
        bits[bitIndex / 64] |= MyBitVector.TWO_POWER[bitIndex % 64];
    }

    public boolean isSet(int bitIndex) {
        return (bits[bitIndex / 64] & MyBitVector.TWO_POWER[bitIndex % 64]) != 0;
    }
}
