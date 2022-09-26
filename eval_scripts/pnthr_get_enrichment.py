import requests
import argparse
import json
import pandas as pd

def main(query_file, ref_file, ann_data_set, save_file):
    base_url = "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?"

    with open(query_file, "r") as f:
        lines = f.readlines()
        lines = [l.rstrip() for l in lines]
        query_list = ",".join(lines)
    # print(query_list)
    
    with open(ref_file, "r") as f:
        lines = f.readlines()
        lines = [l.rstrip() for l in lines]
        ref_list = ",".join(lines)
    # print(ref_list)

    organism = "10090" # Panther taxonomy ID of mouse
    enrichmentTestType = "FISHER"
    correction = "FDR"

    if ann_data_set == "mol":
        ads = "GO:0003674"
    elif ann_data_set == "cell":
        ads = "GO:0005575"
    else:   # bio
        ads = "GO:0008150"

    request_url = (base_url
                  + "geneInputList=" + query_list
                  + "&organism=" + organism
                  + "&refInputList=" + ref_list 
                  + "&refOrganism=" + organism
                  + "&annotDataSet=" + ads
                  + "&enrichmentTestType=" + enrichmentTestType
                  + "&correction=" + correction
    )
    # print(request_url)

    r = requests.get(request_url)
    print(r)
    # print(json.dumps(r.json(), indent=1))
    # print(json.dumps(r.json(), indent=1)[:25])
    df = pd.json_normalize(r.json()["results"]["result"])
    df = df[df["fdr"] < 0.05]
    # df = df.drop(["number_in_list"], axis=1)

    df.to_csv(save_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("query_file", type=str, action="store")
    parser.add_argument("ref_file", type=str, action="store")
    parser.add_argument("ann_data_set", type=str, action="store")
    parser.add_argument("save_file", type=str, action="store")
    args = parser.parse_args()
    if args.ann_data_set not in ["cell", "mol", "bio"]:
        print("ann_data_set must be 'cell', 'mol', or 'bio'")
        exit()
    main(args.query_file, args.ref_file, args.ann_data_set, args.save_file)