#!/bin/bash

# run on the head node
# snakemake -p

# run on the cluster
MAXJOBS=32

# dry
echo "dry run..."
snakemake -n -p -j $MAXJOBS --cluster-config cfg/cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory}"
echo "return to start pipeline, ctrl-c to quit"
read -n 1 -p "Continue?"

# real
echo "starting live run..."
snakemake -p -j $MAXJOBS --cluster-config cfg/cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory}"

