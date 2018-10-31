#!/bin/bash

# run on the head node
# snakemake -p

# run on the cluster
MAXJOBS=64
PREFIX="msi"

rm log/* 

mkdir -p log tmp out

# dry
echo "dry run..."
#snakemake --dryrun -j $MAXJOBS --cluster-config cfg/cluster.json --rerun-incomplete --jobname "msi-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"
echo "return to start pipeline, ctrl-c to quit"
#read -n 1 -p "Continue?"

# real
echo "starting live run..."
mkdir -p out/aggregate
cp cfg/config.yaml out/aggregate/config.yaml
snakemake --verbose -p -j $MAXJOBS --cluster-config cfg/cluster.json --rerun-incomplete --jobname "msi-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"

