#!/bin/bash

# Usage: ./run_simulator.sh <n> <k> <dataset_name>
# Example: ./run_simulator.sh 10000 1000 10k
# ./run_simulator.sh 200 10 a
# ./run_simulator.sh 20000 1000 20k

mkdir -p data 

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <n_taxa> <k_genes> <dataset_name>"
    exit 1
fi

n=$1          # number of taxa
k=$2          # number of genes
dataset=data/${3}    # dataset name / output folder

# Run SimPhy
./simphy_lnx64 \
  -sb f:0.000001 \
  -ld f:0 \
  -lb f:0 \
  -lt f:0 \
  -rs 1 \
  -rl f:${k} \
  -rg 1 \
  -o ${dataset} \
  -sp f:${n} \
  -su f:0.00001 \
  -sg f:1 \
  -sl f:${n} \
  -st f:1000000 \
  -om 1 \
  -v 2 \
  -od 1 \
  -op 1 \
  -oc 1 \
  -on 1 \
  -cs 22 \
  -ll ${n} \
  -ls ${n}

# Run concat_gene_trees.py on replicate 1
python concat_gene_trees.py ${dataset}/1

# Reorganize all trees in dataset
python reorganize_trees.py ${dataset}
