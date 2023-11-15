

# Reporoducibility

## Running EvoGeneX on 24-subline data

### Prep 24-subline data

First, you will need to download the 24-subline data from the Trisicell API and save it locally in a readable format. This is done in prep_data/get_multigene_data.py. By default, this script will use the tpm-normalized data, output files with up to 1000 genes for all 55,401 genes, include only genes at least one non-zero value for each species, and log-transform the values. You can change the maxiumum number of genes in each file using the --inc flag. If you do not want to save all the genes, you can specify the gene range you are interested in using the --start and --end flags. Genes are separated into multiple files for parallelization of running EvoGeneX. If you want a single file with all the genes, use --inc 55401.

### Run EvoGeneX over the genes.

The code to run EvoGeneX is in the run_scripts/ folder. There are three files for running either adaptive or constrained versions: _runner.sh, _driver.sbatch, and _evogenex.R. The _runner.sh file loops over all gene files in the specified folder and either calls the _driver.sbatch or the _evogenex.R scripts for each file individually, depending on if you specify using slurm or the cpu (by default slurm is used). The _driver.sbatch file submits a slurm job to run the _evogenex.R script for the specified file. If you want to submit a single slurm job rather than jobs for all the gene files, you can use this script directly. The _evogenex.R script runs EvoGeneX over all the genes in the specified file. If you want to run EvoGeneX on a single gene file on the cpu, you can call this script directly.

These scripts will save the results in specified file locations.

### Run evaluation scripts over the EvoGeneX results

Evaluation scripts are in the eval_scripts/ folder. These files organize the results and add more statistics columns. The get_gene_results.py script saves files with only positive results and all the results and organizes the columns. The add_stats_cols.R script will add columns for adjusted p-values, and log_2 fold change. The get_all_gene_results.sh will loop through all regimes and run get_gene_results.py and add_stats_cols.R for all results. The output of these files will provide files like the 24 subline results files in the Supplementary Files.

### Cluster the results

The scripts in kmeans/ provides code to cluster the results and make the heatmap in Figure 2. kmeans.py runs kmeans clustering on the results (with a defined random seed for reproducibility). heatmap.py plots the clusters as a heatmap.

## Simulation 

### Simulate data

Scripts to simulate data can be found in the simulation/ folder. simulate_neutral.R will simulate neutral data using a Brownian motion process. simuate_adaptive.R will simulate constrained data (using theta ratio parameter = 0) or adaptive data (theta ratio parameter != 0). You can control all parameter values using the parameters. You can run multiple types of simulations using the _driver.sh files, where you can specify multiple parameter options and it will run simulation scripts for the cross product of all the options. If you want to run this on a cluster, the _runner.sh files start the _driver.sh jobs using slurm.

### Run EvoGeneX on simulated data

You can run EvoGeneX on the simulated data using the same scripts used to run it on the 24 subline data in the run_scripts/ folder. See above for details.

### Run differential expression on simulated data

You can run differential expression on the simulated data using the simulation/differential_expression.py script. 

### Evaluate simulation results

Scripts to organize and evaluate the simulation results are in the eval_scripts/simulation/ folder. get_sim_results.sbatch will run eval_scripts/get_results/get_gene_results.py and eval_scripts/get_results/add_stats_cols.R and generate results files. Then, get_csv.sh will collate the results for all the simulations into a single file. The get_de_results.py will collate the results of differential expression on all simulations into a single file.



---

# Description of folders

## prep_data/

The get_multigene_data.py file will download gene expression data from the Trisicell database and format it for our purposes. It will save the log of the expression data for each gene into files, including only genes with non-zero values for all sublines. It has options to download either the TPM or FPKM normalized data and to split the data into multiple files.

## run_scripts/

This folder includes the scripts used to run EvoGeneX on the gene expression data. It is set up to work either on CPU or create multiple jobs on a cluster using SLURM. If using SlURM, it will create a separate job for each data file in the data directory.

Run the scripts over a directory of data files by calling the "_runner.sh" files. Including the -s flag will run separate jobs a cluster using the "_driver.sbatch" scripts. The "_evogenex.R" files run EvoGeneX on each file.

If you would like to run EvoGeneX on a single file rather than over multiple files, you can call the "_evogenex.R" scripts directly, or use the "_driver.sbatch" script to run it on a cluster.

## eval_scripts/

muts_in_adpt.py will calculate the number of mutated genes that are adaptive.
treatment_ttest.py will perform a t-test between responder and non-responder tumors from the treatment data.

### get_results/

get_gene_results.py will loop through all result files and combine them into a single file. It will also add clear labelling as to whether the gene has adaptive/constrained expression or not.

add_stats_cols.R takes results from get_gene_results.py and adds q-values and log fold values.

get_all_gene_results.sh will loop through the experimental situations and call get_gene_results.py and add_stats_cols.R as needed.


### enrichment/

These scripts will run over-representation enrichment analysis with KEGG or GO on the results from EvoGeneX using the R package clusterProfiler (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html). The functions to run the enrichment are in ora_functions.R. Scripts to run the functions on different data sets are in exp_ora.R and val_ora.R.

## Notes

In all R files, make sure that the library paths are set correctly. If not using the default install location, set the library path using .libPaths("path_to_library"). If using the default, delete this line.