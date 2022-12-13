#/bin/bash

result_path="results/tpm_allrep2/"

for tree in "orig" "lp_scaled"; do
# for tree in "orig"; do
    # for adpt in "agg2" "agg3" "clade2" "1_4_22" "3_10_14"; do
    for adpt in "agg" "clade"; do
        echo "Tree: ${tree}; Adpt: ${adpt}"
        python3 eval_scripts/get_results/get_gene_results.py ${result_path}${tree}/adpt_${adpt}/ a -s ${result_path}${tree}/gene_lists/adpt_${adpt}/
        Rscript eval_scripts/get_results/add_stats_cols.R ${result_path}${tree}/gene_lists/adpt_${adpt}/gene_info.csv
        Rscript eval_scripts/get_results/get_top_genes.R ${result_path}${tree}/gene_lists/adpt_${adpt}/gene_info.csv ${result_path}${tree}/gene_lists/adpt_${adpt}/top_genes.csv
    done
done