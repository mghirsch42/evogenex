#!/bin/bash

#SBATCH --qos=throughput
###SBATCH --exclude=heron[00-11]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=36gb
#SBATCH --exclusive
#SBATCH --time=05:00:00
#SBATCH --output=output/%x_%j.out
#SBATCH --error=output/%x_%j.err

module load R/4.1.2/3.14

TREE_FILE=$1
SINGLE_REGIME_FILE=$2
THREE_REGIME_FILE=$3
DATA_PATH=$4
OUTPUT_FILE=$5

Rscript run_scripts/three_regime_evogenex.R $TREE_FILE $SINGLE_REGIME_FILE $THREE_REGIME_FILE $DATA_PATH $OUTPUT_FILE
