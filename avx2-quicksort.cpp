#include "avx2-partition.cpp"

namespace qs {

    namespace avx2 {

        void quicksort_32(uint32_t* array, int left, int right) {

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
                quicksort_32(array, left, j);
            }

            if (i < right) {
                quicksort_32(array, i, right);
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
            //std::qsort(array, right+1, sizeof(float), compare);

        }


        void quicksort_64( uint64_t *array, int left, int right) {

            int i = left;
            int j = right;

            const uint64_t pivot = array[(i + j) / 2];
            const int AVX2_REGISTER_SIZE = 4; // in 64-bit words

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


        void quicksort_pd(double *array, int left, int right) {

            int i = left;
            int j = right;

            const double pivot = array[(i + j) / 2];
            const int AVX2_REGISTER_SIZE = 4; // in 64-bit words

            if (j - i >= 2 * AVX2_REGISTER_SIZE) {
                qs::avx2::partition_pd(array, pivot, i, j);
            } else {
                scalar_partition_pd(array, pivot, i, j);
            }

            if (left < j) {
                quicksort_pd(array, left, j);
            }

            if (i < right) {
                quicksort_pd(array, i, right);
            }

        }
    } // namespace avx2

} // namespace qs
