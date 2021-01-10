#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <immintrin.h>
#include <iostream>

#include "cmdline.cpp"
#include "input_data.cpp"
#include "quicksort-all.cpp"
#include "avx2-quicksort.cpp"
#include "avx2-partition.cpp"

int main(int argc, char* argv[]) {


    float *niz = new float[66];
    for(int i = 0; i < 66; i++){
        niz[i] = rand();
    }

    qs::avx2::quicksort_ps(niz);
    for(int i=0;i<66;i++){
        std::cout<< niz[i] << " ";
    }




    return 0;
}
