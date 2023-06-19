import pandas as pd
import trisicell as tsc
import time
bwes = tsc.datasets.sublines_bwes()
# print(bwes)
bwes_df = bwes.var
# print(bwes_df.columns)
bwes_df = bwes_df.sort_values(["chrom", "position"])

temp = bwes_df[bwes_df["chrom"] == "chr1"][["position", "alteration"]]
temp.to_csv("mut_list.csv", index=False)
# print(bwes_df.columns)
exit()
print(bwes_df[["chrom", "position"]])
print(len(bwes_df[bwes_df["chrom"] == "chr1"]))
# bwes_df = bwes.sort_values("chrom", axis=1)
# temp = bwes.var.iloc[0:1,:]
# temp.to_csv("temp.csv")
# exit()

# print(bwes.var["position"])
# for cat in bwes.var["reference"].unique():
    # print(cat, (sum(bwes.var["reference"] == cat)))

def add_trip(freqs, ref_trip, mut_trip):
    if (ref_trip, mut_trip) in freqs.keys():
        freqs[(ref_trip, mut_trip)] += 1
    else:
        freqs[(ref_trip, mut_trip)] = 1  

freqs = {}
x = 0
i = 0
mut_idx = int(bwes_df.iloc[i]["position"])
with open("GCA_000001635.8/ncbi_dataset/data/GCA_000001635.8/GCA_000001635.8_GRCm38.p6_genomic.fna", "r") as f:
    # Read the FASTA header
    line = f.readline()
    line = line.split(" ")
    chrom = line[-2][:-1]
    # start = time.perf_counter()

    # Read the chromosome nucleotides
    ref_idx = 1
    old_ref = 1
    curr_line = past_line = f.readline()
    curr_ref_trip = curr_mut_trip = None
    start_time = time.perf_counter()
    while line[0] != ">" and len(bwes_df[bwes_df["chrom"] == "chr1"]) > i:
        # curr_line = curr_line.strip().upper()
        if ref_idx - old_ref > 1000000:
            print(ref_idx)
            old_ref = ref_idx
        # If the mutation was the last character of the previous line, 
        # we need to add the first character of this line
        # if curr_ref_trip:
        #     print("We have previous reference")
        #     q = True
        #     curr_ref_trip += curr_line[0]
        #     curr_mut_trip += curr_line[0]
        #     add_trip(freqs, curr_ref_trip, curr_mut_trip)
        #     curr_ref_trip = curr_mut_trip = None

        # If the next mutation is in this line
        if ref_idx < mut_idx and mut_idx < ref_idx + len(curr_line.strip()):
            print("found", mut_idx)
            print(ref_idx, mut_idx)
            sub_idx = mut_idx - ref_idx
            print(sub_idx)
            print(curr_line)
            if curr_line[sub_idx].upper() != bwes_df.iloc[i]["reference"]:
                print("error at mutation ", i)
                x+=1
                print("reference:", curr_line[sub_idx])
                print("bwes_df:", bwes_df.iloc[i])
            # exit()

            # If the mutation is the first character of this line,
            # we need to get the last character of the previous line
            # if sub_idx == 0:
            #     print("sub index is 0")
            #     ref_trip = past_line[len(past_line-1)] + curr_line[sub_idx:sub_idx+2]
            #     mut_trip = past_line[len(past_line-1)] + bwes_df.iloc[i]["alteration"] + curr_line[sub_idx+1]
            #     add_trip(freqs, ref_trip, mut_trip)
            # If the mutation is the last character of this line,
            # we need to store it and then get the first character of the
            # next line in the next loop
            # elif sub_idx == len(curr_line) - 1:
            #     print("sub index is last")
            #     curr_ref_trip = curr_line[sub_idx-1:sub_idx+1]
            #     curr_mut_trip = curr_line[sub_idx-1] + bwes_df.iloc[i]["alteration"]
            # else:
            #     print("sub index is middle")
            #     ref_trip = curr_line[sub_idx-1:sub_idx+2].upper()
            #     mut_trip = curr_line[sub_idx-1].upper() + bwes_df.iloc[i]["alteration"] + curr_line[sub_idx+1].upper()
            #     add_trip(freqs, ref_trip, mut_trip)
            #     i += 1

            # Get the next mutation index
            # print("Increasing mutation index")
            # mut_idx = int(bwes_df.iloc[i]["position"])
        
        # Get the next reference index and update the lines
        ref_idx = ref_idx + len(curr_line)
        past_line = curr_line
        curr_line = f.readline()
        if x == 3:
            break
print(time.perf_counter() - start_time)
print(ref_idx)
    # end = time.perf_counter()
    # print(len(chromosomes["chr"+chrom]))
    # print(end - start)
    # quit()
    # print(chrom)
    # muts = bwes_df[bwes_df["chrom"] == "chr" + chrom]
    # print(muts)
#     ref_idx = 0
#     while i < 3:
#         line = f.readline().strip()
#         if mut_idx >= ref_idx and mut_idx < ref_idx + len(line):
#             print("found", mut_idx)
#             sub_idx = mut_idx - ref_idx
#             if line[sub_idx].upper() != bwes_df.iloc[i]["reference"]:
#                 print("error at mutation ", i)
#                 print("reference:", line[sub_idx])
#                 print("bwes_df:", bwes_df.iloc[i])
#                 exit()
#             ref_trip = line[sub_idx-1:sub_idx+2].upper()
#             mut_trip = line[sub_idx-1].upper() + bwes_df.iloc[i]["alteration"] + line[sub_idx+1].upper()
#             if (ref_trip, mut_trip) in freqs.keys():
#                 freqs[(ref_trip, mut_trip)] += 1
#             else:
#                 freqs[(ref_trip, mut_trip)] = 1
#             i += 1
#             mut_idx = int(bwes_df.iloc[i]["position"])
#         ref_idx = ref_idx + len(line)
# print(freqs)
# with open("chr1_freqs.csv", "w") as f:
#     for key, val in freqs.items():
#         f.write("{},{}\n".format(key, val))
# print(chromosomes.keys())