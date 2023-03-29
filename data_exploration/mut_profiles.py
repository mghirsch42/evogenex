import pandas as pd
import trisicell as tsc

bwes = tsc.datasets.sublines_bwes()
print(type(bwes.var["amino_acid_"]))
# for cat in bwes.var["reference"].unique():
    # print(cat, (sum(bwes.var["reference"] == cat)))
