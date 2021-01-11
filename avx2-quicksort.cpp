#include "avx2-partition.cpp"

namespace qs {

    namespace avx2 {

        void quicksort(uint32_t* array, int left, int right) {

            int i = left;
            int j = right;

            const uint32_t pivot = array[(i + j)/2];
            const int AVX2_REGISTER_SIZE = 8; // in 32-bit words

            if (j - i >= 2 * AVX2_REGISTER_SIZE) {
                qs::avx2::partition_epi32(array, pivot, i, j);
            } else {
                scalar_partition_epi32(array, pivot, i, j);
            }

            if (left < j) {
                quicksort(array, left, j);
            }

            if (i < right) {
                quicksort(array, i, right);
            }
        }

        void quicksort_ps(float* array, int left, int right) {

            int i = left;
            int j = right;

            const float pivot = array[(i + j)/2];
            const int AVX2_REGISTER_SIZE = 8; // in 32-bit words

            if (j - i >= 2 * AVX2_REGISTER_SIZE) {
                qs::avx2::partition_ps(array, pivot, i, j);
            } else {
                scalar_partition_ps(array, pivot, i, j);
            }

            if (left < j) {
                quicksort_ps(array, left, j);
            }

            if (i < right) {
                quicksort_ps(array, i, right);
            }
        }

        //za 16-bitne
        void quicksort_16( uint16_t *array, int left, int right) {

            int i = left;
            int j = right;

            const uint16_t pivot = array[(i + j) / 2];
            const int AVX2_REGISTER_SIZE = 8; // in 32-bit words

            if (j - i >= 2 * AVX2_REGISTER_SIZE) {
                qs::avx2::partition_epi16(array, pivot, i, j);
            } else {
                scalar_partition_16(array, pivot, i, j);
            }

            if (left < j) {
                quicksort_16(array, left, j);
            }

            if (i < right) {
                quicksort_16(array, i, right);
            }
        }
        void quicksort_64( uint64_t *array, int left, int right) {

            int i = left;
            int j = right;

            const uint64_t pivot = array[(i + j) / 2];
            const int AVX2_REGISTER_SIZE = 8; // in 32-bit words

            if (j - i >= 2 * AVX2_REGISTER_SIZE) {
                qs::avx2::partition_epi64(array, pivot, i, j);
            } else {
                scalar_partition_64(array, pivot, i, j);
            }

            if (left < j) {
                quicksort_64(array, left, j);
            }

            if (i < right) {
                quicksort_64(array, i, right);
            }
        }
    } // namespace avx2

} // namespace qs
