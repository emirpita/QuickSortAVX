#include "partition.cpp"

template<typename NumericType>
void quicksort(NumericType *array, int left, int right) {

    int i = left;
    int j = right;

    const NumericType pivot = array[(i + j) / 2];

    scalar_partition<NumericType>(array, pivot, i, j);

    if (left < j) {
        quicksort<NumericType>(array, left, j);
    }

    if (i < right) {
        quicksort<NumericType>(array, i, right);
    }
}
