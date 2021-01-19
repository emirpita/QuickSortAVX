#!/bin/bash

echo > results.txt

./speed_avx2      100 10 asc | tee -a results.txt
./speed_avx2    10000 10 asc | tee -a results.txt
./speed_avx2   100000 10 asc | tee -a results.txt
./speed_avx2  1000000 10 asc | tee -a results.txt
./speed_avx2  2000000 10 asc | tee -a results.txt
./speed_avx2  5000000 10 asc | tee -a results.txt
./speed_avx2 10000000 10 asc | tee -a results.txt
./speed_avx2 20000000 10 asc | tee -a results.txt

./speed_avx2      100 10 dsc | tee -a results.txt
./speed_avx2    10000 10 dsc | tee -a results.txt
./speed_avx2   100000 10 dsc | tee -a results.txt
./speed_avx2  1000000 10 dsc | tee -a results.txt
./speed_avx2  2000000 10 dsc | tee -a results.txt
./speed_avx2  5000000 10 dsc | tee -a results.txt
./speed_avx2 10000000 10 dsc | tee -a results.txt
./speed_avx2 20000000 10 dsc | tee -a results.txt

./speed_avx2      100 10 random | tee -a results.txt
./speed_avx2    10000 10 random | tee -a results.txt
./speed_avx2   100000 10 random | tee -a results.txt
./speed_avx2  1000000 10 random | tee -a results.txt
./speed_avx2  2000000 10 random | tee -a results.txt
./speed_avx2  5000000 10 random | tee -a results.txt
./speed_avx2 10000000 10 random | tee -a results.txt
./speed_avx2 20000000 10 random | tee -a results.txt

echo "results.txt created"
