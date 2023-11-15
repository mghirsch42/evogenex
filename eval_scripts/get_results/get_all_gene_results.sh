#/bin/bash

##########
### Loop over all regime results and run get_gene_results.py and add_stats_cols.R
##########

result_path="results/tpm/"

for adpt in "har" "has" "las"; do
    echo "Adpt: ${adpt}"
    python3 eval_scripts/get_results/get_gene_results.py ${result_path}${tree}/adpt_${adpt}/ a -s ${result_path}${tree}/gene_lists/adpt_${adpt}/
    Rscript eval_scripts/get_results/add_stats_cols.R ${result_path}${tree}/gene_lists/adpt_${adpt}/gene_info.csv
    Rscript eval_scripts/get_results/add_stats_cols.R ${result_path}${tree}/gene_lists/adpt_${adpt}/all_results.csv
done
