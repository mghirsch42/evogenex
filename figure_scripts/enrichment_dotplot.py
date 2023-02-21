import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os 

method = "gsea"
ds = "go"
subset = "negative"
base_path = "results/tpm_allrep2/orig/gene_lists/"

results = pd.DataFrame()

for regime in ["agg", "clade", "1_4_22", "3_10_14"]:
    if subset == "all":
        fname = base_path + "adpt_{}/enrichment/{}_{}.csv".format(regime, method, ds)
        if not os.path.isfile(fname): continue
        enrich = pd.read_csv(fname)
    else:
        fname = base_path + "adpt_{}/enrichment/{}_{}_{}.csv".format(regime, method, ds, subset)
        if not os.path.isfile(fname): continue
        enrich = pd.read_csv(fname)
    enrich["regime"] = regime
    results = results.append(enrich)

print(len(results))
if len(results) > 50 : results = results.nsmallest(50, "p.adjust")
results["sorter"] = [0] * len(results)
sort_dict = {"agg": 0, "clade": 1, "1_4_22": 2, "3_10_14": 3}
results["sorter"] = [sort_dict[results.iloc[i]["regime"]] for i in range(len(results))]
results = results.sort_values("sorter")
print(results)
plt.figure(figsize=[8,8])
if method == "gsea":
    sns.scatterplot(data = results, x = "regime", y="Description", hue="p.adjust", size="setSize", legend="brief")
else:
    sns.scatterplot(data = results, x = "regime", y="Description", hue="p.adjust", size="Count", legend="brief")
plt.title("{} {} {} adaptive genes".format(method, ds, subset))
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.yticks(fontsize=6)
plt.tight_layout()
plt.savefig("figures/enrichment/dot_{}_{}_{}.png".format(method, ds, subset))
# plt.show()
