#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y


if [[ ( $@ == "--help") ||  $@ == "-h" ]]; then
    echo "Usage: source submit.sh CONFIG_PATH (BIND_PATH) (RUN_NAME)"
    echo
    echo "Submit script for UCL cluster"
    echo
    echo "CONFIG_PATH - Required - path to config file for snakemake run (e.g. config/config.QAPA.yaml)"
    echo "BIND_PATH - Optional - path to externally mount to singularity containers. Default is /SAN/vyplab"
    echo "RUN_NAME - Optional - argument to name run. Config file for run will be copied to folder containing cluster log files (.submissions/<date><time>/) with run name prefixed"
    echo "-h/--help - print this help message and exit"
    exit 0
fi

# Generate string to pass config file to snakemake call (provided path enclosed in double quotes)
CONFIG_PATH=\"$1\"

# Generate string to pass bind directory argument to Singularity (provided path enclosed in double quotes and prefixed with '--bind ')
if [ "$2" != "" ]; then
    IN_BIND_PATH="/SAN/vyplab"
else
    IN_BIND_PATH=$2

# Now enclose in quotes and combine with --bind prefix
SMK_BIND_PATH=\""--bind "${IN_BIND_PATH}\"

if [ "$3" != "" ]; then
    RUN_NAME="QAPA"
else
    RUN_NAME=$1
fi

echo "Constructed config file path: "$CONFIG_PATH
echo "Constructed singularity bind command"$SMK_BIND_PATH

# Create directory for cluster job submission script outputs & copy config file
FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p ${FOLDER}
cp $1 ${FOLDER}/${RUN_NAME}_${1}_config.yaml

snakemake \
-p \
--configfile=${CONFIG_PATH} \
--use-conda \
--use-singularity \
--conda-create-envs-only \
--conda-cleanup-envs \
--singularity-args=${SMK_BIND_PATH} \
--jobscript ucl_cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
--keep-going
