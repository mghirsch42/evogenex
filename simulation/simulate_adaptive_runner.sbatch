#!/bin/bash

#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=36gb
#SBATCH --exclusive
#SBATCH --time=05:00:00
#SBATCH --job-name=sim
#SBATCH --output=output/%x_%j.out
#SBATCH --error=output/%x_%j.err

module load R/4.1.2/3.14

TREE_FILE=$1
REGIME_FILE=$2
OUTPUT_PATH=$3

./simulation/simulate_adaptive_driver.sh \
    $TREE_FILE \
    $REGIME_FILE \
    $OUTPUT_PATH \
