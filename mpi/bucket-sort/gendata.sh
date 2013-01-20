#!/bin/bash

for y in 21474836 214748364 2147483647
do

  for x in 10 50 100 150 200 250
  do
    mpirun -nnp 16 bucket-sort -R $x -n $y > "data/bucket/bucket-sort_"$x"_"$y"_0"

    for i in {1..35}
    do
      mpirun -node 0-$i -nnp 16 bucket-sort -R $x -n $y > "data/bucket/bucket-sort_"$x"_"$y"_"$i
    done

    mpirun -nnp 16 bucket-sort -r -R $x -n $y  > "data/radix/radix-sort_"$x"_"$y"_0"

    for i in {1..35}
    do
      mpirun -node 0-$i -nnp 16 bucket-sort -r -R $x -n $y > "data/radix/radix-sort_"$x"_"$y"_"$i
    done
  done

done