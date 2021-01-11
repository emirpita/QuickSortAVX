#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <immintrin.h>

#include "cmdline.cpp"
#include "input_data.cpp"
#include "quicksort-all.cpp"

template<typename NumericType>
bool is_sorted(NumericType *array, size_t n) {
    assert(n > 0);
    for (size_t i = 1; i < n; i++) {
        if (array[i - 1] > array[i]) {
            printf("mismatch at %lu\n", i);
            return false;
        }
    }

    return true;
}


const size_t AVX2_REGISTER_SIZE = 16;

template<typename NumericType>
class Test {

    bool verbose;

public:
    Test(bool v = true) : verbose(v) {}

    template<typename SORT_FN>
    bool run(SORT_FN sort) {
        const size_t start = 2 * AVX2_REGISTER_SIZE;
        const size_t end   = 256*AVX2_REGISTER_SIZE;

        if (verbose) {
            putchar('\n');
        }

        for (size_t size=start; size < end; size += 1) {

            if (verbose) {
                printf("%lu/%lu\r", size, end);
                fflush(stdout);
            }

            InputAscending<NumericType> asc(size);
            InputDescending<NumericType> dsc(size);
            InputRandom<NumericType> rnd(size);
            InputRandomFew<NumericType> rndfew(size);

            if (!test(sort, asc)) {
                printf("failed for size %lu, intput ascending\n", size);
                return false;
            }

            if (!test(sort, dsc)) {
                printf("failed for size %lu, intput descending\n", size);
                return false;
            }

            if (!test(sort, rnd)) {
                printf("failed for size %lu, intput random\n", size);
                return false;
            }

            if (!test(sort, rndfew)) {
                printf("failed for size %lu, intput random few\n", size);
                return false;
            }
        } // for

        if (verbose) {
            putchar('\n');
        }

        return true;
    }


private:
    template<typename SORT_FN>
    bool test(SORT_FN sort, InputData<NumericType> &data) {
        sort(data.pointer(), 0, data.count() - 1);

        return is_sorted<NumericType>(data.pointer(), data.count());
    }
};


class Flags {
    public:
        bool avx2;

    public:
        Flags(const CommandLine& cmd) {

            enable_all(false);

            bool any_set = false;
            if (cmd.has("-avx2")) {
                avx2 = true;
                any_set = true;
            }


            if (!any_set) {
                enable_all(true);
            }
        }

        void enable_all(bool val) {
            avx2          = val;
        }
};


int main(int argc, char* argv[]) {

    CommandLine cmd(argc, argv);

    puts("Please wait, it might take a while...");
    puts("");

    bool verbose = cmd.has("-v") || cmd.has("--verbose");
    Flags flags(cmd);

    // samo da probamo za int
    Test<uint32_t> test(verbose);
    int ret = EXIT_SUCCESS;


    printf("AVX2 base version... ");
    fflush(stdout);
    if (test.run(qs::avx2::quicksort)) {
        puts("OK");
    } else {
        puts("FAILED");
        ret = EXIT_FAILURE;
    }
    return ret;
}
