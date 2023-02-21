#/bin/bash

result_path="results/simulation/sim_clade_t10_r8/"

for f in $result_path*/; do
    echo $f
    python3 eval_scripts/get_results/get_gene_results.py $f a -s $f
    Rscript eval_scripts/get_results/add_stats_cols.R ${f}gene_info.csv
    Rscript eval_scripts/get_results/add_stats_cols.R ${f}/all_results.csv
done
