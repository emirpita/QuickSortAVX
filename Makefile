.SUFFIXES:
.PHONY: all clean

FLAGS=-std=c++17 -mbmi2 -Wall -pedantic -Wextra
FLAGS_AVX2=$(FLAGS) -mavx2 -msse3 -DHAVE_AVX2_INSTRUCTIONS

DEPS_SORT=partition.cpp \
          avx2-partition.cpp \
          avx2-quicksort.cpp \
          quicksort.cpp

SPEED_DEPS=$(DEPS_SORT) speed.cpp gettime.cpp rdtsc.cpp runtime_stats.cpp
SPEED_FLAGS=-O3 -DNDEBUG

ALL=test_avx2 speed_avx2 speed_avx2_stats

all: $(ALL)

test_avx2: test.cpp input_data.cpp $(DEPS_SORT)
	#$(CXX) $(FLAGS_AVX2) -fsanitize=address test.cpp -o $@
	$(CXX) $(FLAGS_AVX2) test.cpp -o $@

speed_avx2: $(SPEED_DEPS)
	$(CXX) $(FLAGS_AVX2) $(SPEED_FLAGS) speed.cpp -o $@


speed_avx2_stats: $(SPEED_DEPS)
	$(CXX) $(FLAGS_AVX2) $(SPEED_FLAGS) -DWITH_RUNTIME_STATS speed.cpp -o $@

run_avx2: test_avx2
	sde -cnl -- ./$^

clean:
	rm -f $(ALL)
