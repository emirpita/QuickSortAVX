#!/bin/bash

./speed_avx2      100 3 asc | tee -a results.txt
./speed_avx2    10000 3 asc | tee -a results.txt
./speed_avx2   100000 3 asc | tee -a results.txt
./speed_avx2  1000000 3 asc | tee -a results.txt
./speed_avx2  2000000 3 asc | tee -a results.txt
./speed_avx2  5000000 3 asc | tee -a results.txt
./speed_avx2 10000000 3 asc | tee -a results.txt
./speed_avx2 20000000 3 asc | tee -a results.txt

./speed_avx2      100 3 dsc | tee -a results.txt
./speed_avx2    10000 3 dsc | tee -a results.txt
./speed_avx2   100000 3 dsc | tee -a results.txt
./speed_avx2  1000000 3 dsc | tee -a results.txt
./speed_avx2  2000000 3 dsc | tee -a results.txt
./speed_avx2  5000000 3 dsc | tee -a results.txt
./speed_avx2 10000000 3 dsc | tee -a results.txt
./speed_avx2 20000000 3 dsc | tee -a results.txt

./speed_avx2      100 3 randomuniq | tee -a results.txt
./speed_avx2    10000 3 randomuniq | tee -a results.txt
./speed_avx2   100000 3 randomuniq | tee -a results.txt
./speed_avx2  1000000 3 randomuniq | tee -a results.txt
./speed_avx2  2000000 3 randomuniq | tee -a results.txt
./speed_avx2  5000000 3 randomuniq | tee -a results.txt
./speed_avx2 10000000 3 randomuniq | tee -a results.txt
./speed_avx2 20000000 3 randomuniq | tee -a results.txt

echo "results.txt created"
