#/bin/bash

result_path="results/tpm_allrep2/"

for tree in "orig" "lp_scaled"; do
# for tree in "orig"; do
    for adpt in "agg2" "agg3" "clade2"; do
    # for adpt in "agg" "agg2"; do
        echo "Tree: ${tree}; Adpt: ${adpt}"
        python eval_scripts/get_gene_results.py ${result_path}${tree}/adpt_${adpt}/ a -s ${result_path}${tree}/gene_lists/adpt_${adpt}/
    done
done