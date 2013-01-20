#!/bin/bash

for i in {1..48}
do
  ./prefixsum --nproc $i -i -n 1000000000 > "data/prefixsum_1000000000_"$i
  ./prefixsum --nproc $i -i -n 500000000 > "data/prefixsum_500000000_"$i
  ./prefixsum --nproc $i -i -n 100000000 > "data/prefixsum_100000000_"$i
  ./prefixsum --nproc $i -i -n 10000000 > "data/prefixsum_10000000_"$i
  ./prefixsum --nproc $i -i -n 1000000 > "data/prefixsum_10000000_"$i
done

for i in {1..48}
do
  ./prefixsum --nproc $i -p -n 1000000000 > "data/iterative_1000000000_"$i
  ./prefixsum --nproc $i -p -n 500000000 > "data/iterative_500000000_"$i
  ./prefixsum --nproc $i -p -n 100000000 > "data/iterative_100000000_"$i
  ./prefixsum --nproc $i -p -n 10000000 > "data/iterative_10000000_"$i
  ./prefixsum --nproc $i -p -n 1000000 > "data/iterative_10000000_"$i
done