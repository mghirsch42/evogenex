#!/bin/bash

#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=36gb
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=output/%x_%j.out
#SBATCH --error=output/%x_%j.err

#module load Python3/3.8.1
#source ../env/bin/activate

#python3 tools/differential_expression.py
./eval_scripts/get_results/simulation/get_csv.sh
