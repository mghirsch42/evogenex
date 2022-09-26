#/bin/bash

result_path="results/tpm_allrep2/"
ref_gene_list="${result_path}ref_gene_list.csv"

# for tree in "orig" "lp_scaled"; do
for tree in "orig"; do
    # for adpt in "agg" "agg2" "agg3" "clade" "clade2"; do
    for adpt in "agg" "agg2"; do
        for ann in "cell" "bio" "mol"; do
            echo "Tree: ${tree}; Adpt: ${adpt}; Ann: ${ann}"
            python eval_scripts/pnthr_get_enrichment.py ${result_path}${tree}/gene_lists/adpt_${adpt}/genes.csv ${ref_gene_list} ${ann} ${result_path}${tree}/gene_lists/adpt_${adpt}/pnthr_${ann}.csv
        done
    done
done