package samosa.utils;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.javatuples.Pair;

/**
 *
 * @author giulia
 */
public class Utils {

    /**
     * 
     * @param rnd Random object
     * @param set collection of integers
     * @return a random element in set
     */
    public static int getRandomSetElement(Random rnd, Collection<Integer> set) {
        return set.stream().skip(rnd.nextInt(set.size())).findFirst().orElse(-1);
    }
    
    /**
     * Converts an array of sets into an array of arrays and
     * computes the sum of the sizes of all the sets.
     * @param input_set hash set
     * @param output_array array to copy set content
     * @return sum of the sizes of the sets in input_set
     */
    public static int convertHashSetIntoArray(IntOpenHashSet[] input_set, int[][] output_array) {
        int sum_of_sizes = 0;
        for (int i = 0; i < input_set.length; i++) {
            output_array[i] = input_set[i].toIntArray();
            Arrays.sort(output_array[i]);
            sum_of_sizes += output_array[i].length;
        }
        return sum_of_sizes;
    }
    
    /**
     * Converts an array of lists into an array of arrays and
     * computes the sum of the sizes of all the sets.
     * @param input_array array list
     * @param output_array array to copy list content
     * @return 
     */
    public static int convertTempLeftRightInDefinitiveLeftRight(IntArrayList[] input_array, int[][] output_array) {
        int sum_of_sizes = 0;
        for (int i = 0; i < input_array.length; i++) {
            input_array[i].sort(null);
            output_array[i] = input_array[i].toIntArray();
            sum_of_sizes += output_array[i].length;
        }
        return sum_of_sizes;
    }
    
    /**
     * Copy the content of the first array into the second array.
     * @param from_array first array
     * @param to_array second array
     */
    public static void initializeArrayOfArrays(final int[][] from_array, final int[][] to_array) {
        for (int i = 0; i < from_array.length; i++) {
            to_array[i] = Arrays.copyOf(from_array[i], from_array[i].length);
        }
    }

}
