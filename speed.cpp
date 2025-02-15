#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <memory>

#include "cmdline.cpp"
#ifdef WITH_RUNTIME_STATS
#   include "runtime_stats.cpp"
#endif
#include "input_data.cpp"
#include "quicksort-all.cpp"
#define USE_RDTSC // undef to get measurments in seconds
#ifdef USE_RDTSC
#   include "rdtsc.cpp"
#else
#   include "gettime.cpp"
#endif



template<typename NumericType>
class PerformanceTest final {

    int iterations;
    InputData<NumericType> &input;
    NumericType *tmp;

public:
    PerformanceTest(int n, InputData<NumericType> &input)
            : iterations(n), input(input) {

        assert(iterations > 0);
        tmp = new NumericType[input.count()];
    }

    ~PerformanceTest() {
        delete[] tmp;
    }

public:
    template <typename SORT_FUNCTION>
    uint64_t run(SORT_FUNCTION sort) {

        uint64_t time = 0;

        int k = iterations;
        while (k--) {
            memcpy(tmp, input.pointer(), input.size());

            uint64_t t1, t2;

#ifdef USE_RDTSC
            RDTSC_START(t1);
#else
            t1 = get_time();
#endif
            sort(input.pointer(), 0, input.count() - 1);
            quicksort(input.pointer(), 0, input.count() - 1);

#ifdef USE_RDTSC
            RDTSC_START(t2);
#else
            t2 = get_time();
#endif

            const uint64_t dt = t2 - t1;

            if (time == 0) {
                time = dt;
            } else if (dt < time) {
                time = dt;
            }
        }

        return time;
    }
};


enum class InputType {
    randomfew,
    randomuniq,
    random,
    ascending,
    descending,
};


const char* as_string(InputType type) {
    switch (type) {
        case InputType::randomfew:
            return "randomfew";

        case InputType::randomuniq:
            return "randomuniq";

        case InputType::random:
            return "random";

        case InputType::ascending:
            return "ascending";

        case InputType::descending:
            return "descending";

        default:
            return "<unknown>";
    }
}

template<typename NumericType>
void std_qsort_wrapper(NumericType *array, int left, int right) {

    std::qsort(array + left, right - left + 1, sizeof(NumericType), [](const void *a, const void *b) {
        NumericType a1 = *static_cast<const NumericType *>(a);
        NumericType a2 = *static_cast<const NumericType *>(b);

        if (a1 < a2) return -1;
        if (a1 > a2) return 1;
        return 0;
    });
}

template<typename NumericType>
void std_stable_sort_wrapper(NumericType *array, int left, int right) {

    std::stable_sort(array + left, array + right + 1);
}

template<typename NumericType>
void std_sort_wrapper(NumericType *array, int left, int right) {

    std::sort(array + left, array + right + 1);
}


class Flags {
    public:
        bool std_sort;
        bool std_qsort;
        bool std_stable_sort;
        bool quicksort;
        bool avx2;

    public:
        Flags(const CommandLine& cmd) {

            enable_all(false);

            bool any_set = false;
            if (cmd.has("-std-sort")) {
                std_sort = true;
                any_set = true;
            }

            if (cmd.has("-std-qsort")) {
                std_qsort = true;
                any_set = true;
            }

            if (cmd.has("-std-stable-sort") || cmd.has("-std-stable")) {
                std_stable_sort = true;
                any_set = true;
            }

            if (cmd.has("-quicksort_32")) {
                quicksort = true;
                any_set = true;
            }

            if (cmd.has("-avx2")) {
                avx2 = true;
                any_set = true;
            }

            if (!any_set) {
                enable_all(true);
            }
        }

        void enable_all(bool val) {
            std_sort      = val;
            std_qsort     = val;
            std_stable_sort = val;
            quicksort     = val;
            avx2          = val;
        }
};

template<typename NumericType>
class Test {

    std::unique_ptr<InputData<NumericType>> data;
    InputType type;
    size_t count;
    int iterations;
    Flags flags;
    uint64_t ref;

public:
    Test(InputType type, size_t count, int iterations, Flags&& flags)
        : type(type)
        , count(count)
        , iterations(iterations)
        , flags(std::move(flags)) {

        switch (type) {
            case InputType::randomfew:
                data.reset(new InputRandomFew<NumericType>(count));
                break;

            case InputType::randomuniq:
                data.reset(new InputRandomUnique<NumericType>(count));
                break;

            case InputType::random:
                data.reset(new InputRandom<NumericType>(count));
                break;

            case InputType::ascending:
                data.reset(new InputAscending<NumericType>(count));
                break;

            case InputType::descending:
                data.reset(new InputDescending<NumericType>(count));
                break;
        }
    }

    void run() {

        printf("<------------------------------------------->");
        printf("items count: %lu (%lu bytes), input %s\n", data->count(), data->size(), as_string(type));

        ref = 0;

        if (flags.std_sort) {
            measure("std::sort", std_sort_wrapper<NumericType>);
        }

        if (flags.std_qsort) {
            measure("std::qsort", std_qsort_wrapper<NumericType>);
        }

        if (flags.std_stable_sort) {
            measure("std::stable_sort", std_stable_sort_wrapper<NumericType>);
        }

        if (flags.std_qsort) {
            measure("quick sort", quicksort<NumericType>);
        }

        if (flags.avx2) {
        }
    }

    friend void runUint8_t(Test<uint8_t> &test);

    friend void runUint16_t(Test<uint16_t> &test);

    friend void runUint32_t(Test<uint32_t> &test);

    friend void runUint64_t(Test<uint64_t> &test);

    friend void runFloat(Test<float> &test);

    friend void runDouble(Test<double> &test);


private:
    template<typename SORT_FUNCTION>
    void measure(const char *name, SORT_FUNCTION sort) {

        PerformanceTest<NumericType> test(iterations, *data);

        printf("%30s ... ", name);
        fflush(stdout);
#ifdef WITH_RUNTIME_STATS
        statistics.reset();
#endif
        uint64_t time = test.run(sort);

#ifdef USE_RDTSC
        printf("%10lu cycles", time);
        if (ref > 0) {
            printf(" (%0.2f)", ref/double(time));
        }

#   ifdef WITH_RUNTIME_STATS
        if (statistics.anything_collected()) {
            printf("\n");
            printf("\t\tpartition calls: %lu (+%lu scalar)\n", statistics.partition_calls, statistics.scalar__partition_calls);
            printf("\t\titems processed: %lu (+%lu by scalar partition)\n", statistics.items_processed, statistics.scalar__items_processed);

            const size_t total_items = statistics.items_processed + statistics.scalar__items_processed;

            if (total_items != 0) {
                const double cpi = double(time)/total_items;
                printf("\t\t               : %0.4f cycles/item\n", cpi * iterations);
            }

            if (!statistics.pvbyte_histogram.empty()) {
                puts("Histogram for pvbyte:");
                statistics.pvbyte_histogram.print();
            }
        }
#   endif // WITH_RUNTIME_STATS

#else
        printf("%0.4f s", time/1000000.0);
        if (ref > 0) {
            printf(" (%0.2f)\n", ref / double(time));
        }
#endif
        putchar('\n');

        if (ref == 0) {
            ref = time;
        }
    }
}; // Test class template


void runUint32_t(Test<uint32_t> &test) {

    printf("\n<------------------------------------------->\n");
    printf("items count: %lu (%lu bytes), input %s\n", test.data->count(), test.data->size(), as_string(test.type));

    test.ref = 0;

    if (test.flags.std_sort) {
        test.measure("std::sort", std_sort_wrapper<uint32_t>);
    }

    if (test.flags.std_qsort) {
        test.measure("std::qsort", std_qsort_wrapper<uint32_t>);
    }

    if (test.flags.std_stable_sort) {
        test.measure("std::stable_sort", std_stable_sort_wrapper<uint32_t>);
    }

    if (test.flags.std_qsort) {
        test.measure("quick sort", quicksort<uint32_t>);
    }

    // AVX2 sortiranje s obzirom na tip

    if (test.flags.avx2) {
        test.measure("AVX2 quick sort for 32-bit integer", qs::avx2::quicksort_32);
    }
}

void runUint64_t(Test<uint64_t> &test) {

    printf("\n<------------------------------------------->\n");
    printf("items count: %lu (%lu bytes), input %s\n", test.data->count(), test.data->size(), as_string(test.type));

    test.ref = 0;

    if (test.flags.std_sort) {
        test.measure("std::sort", std_sort_wrapper<uint64_t>);
    }

    if (test.flags.std_qsort) {
        test.measure("std::qsort", std_qsort_wrapper<uint64_t>);
    }

    if (test.flags.std_stable_sort) {
        test.measure("std::stable_sort", std_stable_sort_wrapper<uint64_t>);
    }

    if (test.flags.std_qsort) {
        test.measure("quick sort", quicksort<uint64_t>);
    }

    // AVX2 sortiranje s obzirom na tip

    if (test.flags.avx2) {
        test.measure("AVX2 quick sort for 64-bit integer", qs::avx2::quicksort_64);
    }
}

void runFloat(Test<float> &test) {

    printf("\n<------------------------------------------->\n");
    printf("items count: %lu (%lu bytes), input %s\n", test.data->count(), test.data->size(), as_string(test.type));

    test.ref = 0;

    if (test.flags.std_sort) {
        test.measure("std::sort", std_sort_wrapper<float>);
    }

    if (test.flags.std_qsort) {
        test.measure("std::qsort", std_qsort_wrapper<float>);
    }

    if (test.flags.std_stable_sort) {
        test.measure("std::stable_sort", std_stable_sort_wrapper<float>);
    }

    if (test.flags.std_qsort) {
        test.measure("quick sort", quicksort<float>);
    }

    // AVX2 sortiranje s obzirom na tip

    if (test.flags.avx2) {
        test.measure("AVX2 quick sort for float", qs::avx2::quicksort_ps);
    }
}

void runDouble(Test<double> &test) {

    printf("\n<--------------------------------------------------->\n");
    printf("items count: %lu (%lu bytes), input %s\n", test.data->count(), test.data->size(), as_string(test.type));

    test.ref = 0;

    if (test.flags.std_sort) {
        test.measure("std::sort", std_sort_wrapper<double>);
    }

    if (test.flags.std_qsort) {
        test.measure("std::qsort", std_qsort_wrapper<double>);
    }

    if (test.flags.std_stable_sort) {
        test.measure("std::stable_sort", std_stable_sort_wrapper<double>);
    }

    if (test.flags.std_qsort) {
        test.measure("quick sort", quicksort<double>);
    }

    // AVX2 sortiranje s obzirom na tip

    if (test.flags.avx2) {
        test.measure("AVX2 quick sort for double", qs::avx2::quicksort_pd);
    }
}



// ------------------------------------------------------------


void usage() {
    puts("usage:");
    puts("speed SIZE ITERATIONS INPUT [options]");
    puts("");
    puts("where");
    puts("* SIZE       - number of elements");
    puts("* ITERATIONS - number of iterations");
    puts("* INPUT      - one of:");
    puts("                 ascending (or asc)");
    puts("                 descending (or dsc, desc)");
    puts("                 random (or rnd, rand)");
    puts("                 randomfew");
    puts("                 randomuniq");
    puts("options      - optional name of procedure(s) to run");
}


int main(int argc, char* argv[]) {

    if (argc < 4) {
        usage();
        return EXIT_FAILURE;
    }

    int count      = atoi(argv[1]);
    int iterations = atoi(argv[2]);
    InputType type;

#define is_keyword(key) (strcmp(argv[3], key) == 0)
    if (is_keyword("descending") || is_keyword("desc") || is_keyword("dsc")) {
        type = InputType::descending;
    } else if (is_keyword("ascending") || is_keyword("asc")) {
        type = InputType::ascending;
    } else if (is_keyword("random") || is_keyword("rnd") || is_keyword("rand")) {
        type = InputType::random;
    } else if (is_keyword("randomfew")) {
        type = InputType::randomfew;
    } else if (is_keyword("randomuniq")) {
        type = InputType::randomuniq;
    } else {
        usage();
        return EXIT_FAILURE;
    }
#undef is_keyword


#ifdef USE_RDTSC
    RDTSC_SET_OVERHEAD(rdtsc_overhead_func(1), iterations);
#endif
    CommandLine cmd(argc, argv);
    Flags flags(cmd);

    Test<uint32_t> test32b(type, count, iterations, std::move(flags));
        runUint32_t(test32b);


    Test<uint64_t> test64b(type, count, iterations, std::move(flags));
        runUint64_t(test64b);


    Test<float> testFloat(type, count, iterations, std::move(flags));
        runFloat(testFloat);

    Test<double> testDouble(type, count, iterations, std::move(flags));
        runDouble(testDouble);
    return EXIT_SUCCESS;
}
