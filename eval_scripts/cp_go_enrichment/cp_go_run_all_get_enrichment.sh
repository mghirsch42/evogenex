#/bin/bash

result_path="results/tpm_allrep2/"
ref_gene_list="${result_path}ref_gene_info.csv"

# for tree in "orig" "lp_scaled"; do
for tree in "orig"; do
    # for adpt in "agg" "agg2" "agg3" "clade" "clade2"; do
    # for adpt in "agg" "clade" "1_4_22" "3_10_14"; do
    for adpt in "3_10_14"; do
        echo "Tree: ${tree}; Adpt: ${adpt}"
        Rscript eval_scripts/cp_go_enrichment/cp_go_get_enrichment.R ${result_path}${tree}/gene_lists/adpt_${adpt}/gene_info.csv ${ref_gene_list} ${result_path}${tree}/gene_lists/adpt_${adpt}/cp_go.csv
    done
done