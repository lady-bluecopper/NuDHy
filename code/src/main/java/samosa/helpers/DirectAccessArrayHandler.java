package samosa.helpers;

import java.util.Random;

public class DirectAccessArrayHandler extends ArrayHandler {

    // indices of the values in the first array and not in the second
    protected int[] index_for_values_in_FIRST_and_not_in_SECOND;
    // indices of the values in the second array and not in the first
    protected int[] index_for_values_in_SECOND_and_not_in_FIRST;
    // array to store the values chosen from the two arrays
    public final int[] result_indexes;

    public DirectAccessArrayHandler() {
        super();
        this.index_for_values_in_FIRST_and_not_in_SECOND = new int[0];
        this.index_for_values_in_SECOND_and_not_in_FIRST = new int[0];
        this.result_indexes = new int[]{-1, -1};
    }

    /**
     * Replaces b with a in values and ensures that values remains sorted.
     * @param values list of neighbors of the node
     * @param a value currently in values
     * @param b new value to add to values
     * @param index_of_a_in_values position of a in values
     */
    public void fastReplaceInSortedArray(final int[] values, final int a, final int b, final int index_of_a_in_values) {
        int index = index_of_a_in_values;
        values[index] = b;
        if (values.length <= 1) {
            return;
        }
        int tmp;
        while (true) {
            // edge case 1.
            if (index == 0) {
                if (values[index] <= values[index + 1]) {
                    return;
                }
                // swap right.
                tmp = values[index + 1];
                values[index + 1] = values[index];
                values[index] = tmp;
                index++;
                continue;
            }
            // edge case 2.
            if (index == values.length - 1) {
                if (values[index - 1] <= values[index]) {
                    return;
                }
                // swap left.
                tmp = values[index - 1];
                values[index - 1] = values[index];
                values[index] = tmp;
                index--;
                continue;
            }
            // central case.
            if (values[index - 1] > values[index]) {
                // swap left.
                tmp = values[index - 1];
                values[index - 1] = values[index];
                values[index] = tmp;
                index--;
                continue;
            }
            if (values[index] > values[index + 1]) {
                // swap right.
                tmp = values[index + 1];
                values[index + 1] = values[index];
                values[index] = tmp;
                index++;
                continue;
            }
            return;
        }
    }

    /**
     * Compute the symmetric difference between the two arrays values_a and values_b.
     * The method assumes that the two arrays are sorted.
     * @param values_a first array of values
     * @param values_b second array of values
     */
    @Override
    public void computeSymmetricDifference(final int[] values_a, final int[] values_b) {
        //
        int i = 0;
        int j = 0;
        // reinitialize the data structures, as this object is used multiple times
        if (this.values_in_FIRST_and_not_in_SECOND.length < values_a.length) {
            this.values_in_FIRST_and_not_in_SECOND = new int[values_a.length];
            this.index_for_values_in_FIRST_and_not_in_SECOND = new int[this.values_in_FIRST_and_not_in_SECOND.length];
        }
        if (this.values_in_SECOND_and_not_in_FIRST.length < values_b.length) {
            this.values_in_SECOND_and_not_in_FIRST = new int[values_b.length];
            this.index_for_values_in_SECOND_and_not_in_FIRST = new int[this.values_in_SECOND_and_not_in_FIRST.length];
        }
        this.size_values_in_FIRST_and_not_in_SECOND = 0;
        this.size_values_in_SECOND_and_not_in_FIRST = 0;
        // iterate over the two arrays
        while (true) {
            if (i >= values_a.length) {
                break;
            }
            if (j >= values_b.length) {
                break;
            }
            if (values_a[i] == values_b[j]) {
                j++;
                i++;
                continue;
            }
            // we found a value in values_b that is not in values_a
            if (values_a[i] > values_b[j]) {
                this.values_in_SECOND_and_not_in_FIRST[this.size_values_in_SECOND_and_not_in_FIRST] = values_b[j];
                this.index_for_values_in_SECOND_and_not_in_FIRST[this.size_values_in_SECOND_and_not_in_FIRST++] = j;
                j++;
                continue;
            }
            // we found a value in values_a that is not in values_b
            this.values_in_FIRST_and_not_in_SECOND[this.size_values_in_FIRST_and_not_in_SECOND] = values_a[i];
            this.index_for_values_in_FIRST_and_not_in_SECOND[this.size_values_in_FIRST_and_not_in_SECOND++] = i;
            i++;
        }
        if (i < values_a.length) {
            // all these values are only in values_a
            for (; i < values_a.length; i++) {
                this.values_in_FIRST_and_not_in_SECOND[this.size_values_in_FIRST_and_not_in_SECOND] = values_a[i];
                this.index_for_values_in_FIRST_and_not_in_SECOND[this.size_values_in_FIRST_and_not_in_SECOND++] = i;
            }
        }
        if (j < values_b.length) {
            // all these values are only in values_b
            for (; j < values_b.length; j++) {
                this.values_in_SECOND_and_not_in_FIRST[this.size_values_in_SECOND_and_not_in_FIRST] = values_b[j];
                this.index_for_values_in_SECOND_and_not_in_FIRST[this.size_values_in_SECOND_and_not_in_FIRST++] = j;
            }
        }
    }

    /**
     * Choose a random element in values_a that is not in values_b, and a
     * random element in values_b that is not in values_a.
     * @param values_a first array of values
     * @param values_b second array of values
     * @param rnd Random object
     */
    @Override
    public void pickTwoRandomElementsInTheXor(final int[] values_a, final int[] values_b, Random rnd) {
        // computes the symmetric difference between the two arrays
        this.computeSymmetricDifference(values_a, values_b);
        this.result[0] = -1;
        this.result[1] = -1;
        if (this.size_values_in_FIRST_and_not_in_SECOND == 0 || this.size_values_in_SECOND_and_not_in_FIRST == 0) {
            // there is no value in values_a that is not in values_b
            // or no value in values_b that is not in values_a
            return;
        }
        // choose a pair of values uniformly at random
        int first_index = rnd.nextInt(this.size_values_in_FIRST_and_not_in_SECOND);
        int second_index = rnd.nextInt(this.size_values_in_SECOND_and_not_in_FIRST);
        this.result[0] = this.values_in_FIRST_and_not_in_SECOND[first_index];
        this.result_indexes[0] = this.index_for_values_in_FIRST_and_not_in_SECOND[first_index];
        this.result[1] = this.values_in_SECOND_and_not_in_FIRST[second_index];
        this.result_indexes[1] = this.index_for_values_in_SECOND_and_not_in_FIRST[second_index];
    }

}
