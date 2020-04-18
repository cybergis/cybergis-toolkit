#!/bin/sh


for i in $(seq 2 1 6)
do
  for j in $(seq 11 1 15)
  do
    sbatch execute.sh $i $j
  done
done

