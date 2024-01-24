package samosa.helpers;

import java.util.Arrays;
import java.util.Random;

public class ArrayHandler {

    // values in the first array that are not in the second
    protected int[] values_in_FIRST_and_not_in_SECOND;
    // number of values in the first array and not in the second
    protected int size_values_in_FIRST_and_not_in_SECOND;
    // values in the second array that are not in the first
    protected int[] values_in_SECOND_and_not_in_FIRST;
    // number of values in the second array that are not in the first
    protected int size_values_in_SECOND_and_not_in_FIRST;
    // array to store the elements chosen from the two arrays
    public final int[] result;

    public ArrayHandler() {
        this.values_in_FIRST_and_not_in_SECOND = new int[0];
        this.size_values_in_FIRST_and_not_in_SECOND = 0;
        this.values_in_SECOND_and_not_in_FIRST = new int[0];
        this.size_values_in_SECOND_and_not_in_FIRST = 0;
        this.result = new int[]{-1, -1};
    }

    /**
     * Sorts the array of values.
     * @param values array of values
     */
    public void sort(final int[] values) {
        Arrays.sort(values);
    }

    /**
     * 
     * @param values array of values
     * @param value value to search
     * @return position of value in values
     */
    public int getIndexToAddValueKeepingTheArraySorted(final int[] values, final int value) {
        int index = this.findKeyWithBinarySearch(values, value);
        index = (index < 0 ? -1 * index - 1 : index);
        return index;
    }

    /**
     * 
     * @param values array of values
     * @param value value to search
     * @return position of value in values
     */
    public int getIndexToAddValueKeepingTheArraySorted(final long[] values, final long value) {
        int index = this.findKeyWithBinarySearch(values, value);
        index = (index < 0 ? -1 * index - 1 : index);
        return index;
    }

    /**
     * 
     * @param values array of values
     * @param o value to search
     * @return position of o in values
     */
    protected int findKeyWithBinarySearch(final int[] values, final int o) {
        int midVal;
        int from = 0;
        int to = values.length - 1;
        int mid;
        while (from <= to) {
            mid = (from + to) >>> 1;
            midVal = values[mid];
            if (midVal < o) {
                from = mid + 1;
            } else if (midVal > o) {
                to = mid - 1;
            } else {
                return mid;
            }
        }
        return -(from + 1);
    }

    /**
     * 
     * @param values array of values
     * @param o value to search
     * @return position of o in values
     */
    protected int findKeyWithBinarySearch(final long[] values, final long o) {
        long midVal;
        int from = 0;
        int to = values.length - 1;
        int mid;
        while (from <= to) {
            mid = (from + to) >>> 1;
            midVal = values[mid];
            if (midVal < o) {
                from = mid + 1;
            } else if (midVal > o) {
                to = mid - 1;
            } else {
                return mid;
            }
        }
        return -(from + 1);
    }

    /**
     * Replaces b with a in values and ensures that values remains sorted.
     * @param values array of values
     * @param a value currently in values
     * @param b value to add to values
     */
    public void fastReplaceInSortedArray(final int[] values, final int a, final int b) {
        int index = this.findKeyWithBinarySearch(values, a);
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
     * Computes the number of elements in common between two arrays.
     * The method assumes that the two arrays are sorted.
     * @param values_a first array of values
     * @param values_b second array of values
     * @return number of elements in common between the two arrays
     */
    public int computeIntersectionSize(final int[] values_a, final int[] values_b) {
        int num_common_elements = 0;
        int i = 0;
        int j = 0;
        while (true) {
            if (i >= values_a.length) {
                break;
            }
            if (j >= values_b.length) {
                break;
            }
            if (values_a[i] == values_b[j]) {
                num_common_elements++;
                j++;
                i++;
                continue;
            }
            if (values_a[i] > values_b[j]) {
                j++;
                continue;
            }
            i++;
        }
        return num_common_elements;
    }

    /**
     * 
     * @param values_a first array of values
     * @param values_b second array of values
     * @return true if all the values in values_a are contained in values_b
     */
    public boolean isAcontainedInB(final int[] values_a, final int[] values_b) {
        if (values_a.length > values_b.length) {
            return false;
        }
        int midVal;
        int from = 0;
        int to = values_b.length - 1;
        int mid = 0;
        int v_a;
        for (int i = 0; i < values_a.length; i++) {
            v_a = values_a[i];
            //
            while (from <= to) {
                mid = (from + to) >>> 1;
                midVal = values_b[mid];
                if (midVal < v_a) {
                    from = mid + 1;
                } else if (midVal > v_a) {
                    to = mid - 1;
                } else {
                    // v_a is present
                    break;
                }
            }
            if (from > to) {
                return false;
            }
            to = values_b.length - 1;
            from = mid + 1;
        }
        return true;
    }

    /**
     * Compute the symmetric difference between the two arrays values_a and values_b.
     * The method assumes that the two arrays are sorted.
     * @param values_a first array of values
     * @param values_b second array of values
     */
    public void computeSymmetricDifference(final int[] values_a, final int[] values_b) {
        //
        int i = 0;
        int j = 0;
        // reinitialize the data structures, as this object is used multiple times
        if (this.values_in_FIRST_and_not_in_SECOND.length < values_a.length) {
            this.values_in_FIRST_and_not_in_SECOND = new int[values_a.length];
        }
        if (this.values_in_SECOND_and_not_in_FIRST.length < values_b.length) {
            this.values_in_SECOND_and_not_in_FIRST = new int[values_b.length];
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
                this.values_in_SECOND_and_not_in_FIRST[this.size_values_in_SECOND_and_not_in_FIRST++] = values_b[j];
                j++;
                continue;
            }
            // we found a value in values_a that is not in values_b
            this.values_in_FIRST_and_not_in_SECOND[this.size_values_in_FIRST_and_not_in_SECOND++] = values_a[i];
            i++;
        }
        if (i < values_a.length) {
            // all these values are only in values_a
            for (; i < values_a.length; i++) {
                this.values_in_FIRST_and_not_in_SECOND[this.size_values_in_FIRST_and_not_in_SECOND++] = values_a[i];
            }
        }
        if (j < values_b.length) {
            // all these values are only in values_b
            for (; j < values_b.length; j++) {
                this.values_in_SECOND_and_not_in_FIRST[this.size_values_in_SECOND_and_not_in_FIRST++] = values_b[j];
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
        this.result[0] = this.values_in_FIRST_and_not_in_SECOND[rnd
                .nextInt(this.size_values_in_FIRST_and_not_in_SECOND)];
        this.result[1] = this.values_in_SECOND_and_not_in_FIRST[rnd
                .nextInt(this.size_values_in_SECOND_and_not_in_FIRST)];
    }

}
