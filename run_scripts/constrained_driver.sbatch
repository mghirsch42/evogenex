#!/bin/bash

#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=36gb
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=output/%x_%j.out
#SBATCH --error=output/%x_%j.err

module load R/4.1.2/3.14

TREE_FILE=$1
REGIME_FILE=$2
DATA_PATH=$3
OUTPUT_FILE=$4

Rscript run_scripts/constrained_evogenex.R \
  $TREE_FILE \
  $REGIME_FILE \
  $DATA_PATH \
  $OUTPUT_FILE
