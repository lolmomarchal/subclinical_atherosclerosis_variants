import pybedtools
import pandas as pd
import numpy as np
import argparse

"""
CPG variants: CG sites are regions of DNA where a cytosine nucleotide is followed by 
a guanine nucleotide in the linear sequence of bases along its 5' → 3' direction (Positive strand). CpG sites can creates spots of potential
methylation and are important in the context of epigenetics. 

CPG variants are characterized by those in which there is some kind of removal or insertion of a CG site when comparing a variant vs the reference genome. 

EX:

REF: ACTGTGCT
Variant 1: insertion of GT -> AC<GT>TGTGCT (creates one new CG spot)
Variant 2: deletion of CT -> ACTGTG (removes one CG spot)

Reference genome is in the context of the positive strand. 

"""
def reverse_comp(seq: str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in seq)

def find_valid_instance(region, ref):
    mid_index = len(region) // 2  
    ref_len = len(ref)
    start_range = max(0, mid_index - ref_len + 1)
    end_range = min(len(region) - ref_len + 1, mid_index + 1)
    
    for start in range(start_range, end_range):
        substring = region[start:start + ref_len]
        if substring == ref and region[mid_index] in substring:  # Check if ref matches and contains middle character (REF)
            return start 
    
    return -1 

def replace_valid_instance(region, ref, alt):
    start_idx = find_valid_instance(region, ref)
    if start_idx == -1:  
        print("region not found")
        return region
    
    ref_len = len(ref)

    new_region = region[:start_idx] + alt + region[start_idx + ref_len:]
    
    return new_region

# to get regions of 5 nucleotides before and after pos
# more in case of deletions
def get_pre_ref_post_regions(df: str, fasta: str) -> pd.DataFrame:
  if ".xlsx" in df:
    df = pd.read_excel(df)
  else:
    df = pd.read_csv(df)
    
  for i, row in df.iterrows():
    # from 1 to 0 indexing 
    pos = row["POS"] -1
    ref = row["REF"]
    alt = row["ALT"]
    strand = row["STRAND"]
    variant = row["VARIANT_CLASS"]

    # Determine the regions to extract based on the variant type  
    if variant == "deletion":
        length_ = len(ref)
        if length_ < 5:
            post = 5
            chr_pos = f"{row['CHR']} {pos-5} {pos+6}"
        else:
            post = length_
            chr_pos = f"{row['CHR']} {pos-length_} {pos+length_+1}"
    else:
        post = 5
        chr_pos = f"{row['CHR']} {pos-5} {pos+6}"
    

    # Create the BedTool object and extract the sequence
    a = pybedtools.BedTool(chr_pos, from_string=True)
    a = a.sequence(fi=fasta)
    nt = open(a.seqfn).read()
    seq = nt.split("\n")[1]

    # Extract the relevant regions
    df.loc[i, "prev_region"] = seq[0:post]
    df.loc[i, "post_region"] = seq[post+1:]
    df.loc[i, "ref_region"] = seq[post]
  return df
def get_normal_alt_regions(df: pd.DataFrame):
  for i, row in df.iterrows():
    pos = row["POS"]
    ref = row["REF"]
    alt = row["ALT"]
    variant = row["VARIANT_CLASS"]
    previous_reg = row["prev_region"]
    post_reg = row["post_region"]
    ref_region = row["ref_region"]

    # get the entire normal region 
    region = previous_reg + ref_region + post_reg
    df.loc[i, "region_normal"] = region
    
    #get region alt
    if variant != "deletion":
        df.loc[i, "region_alt"] = previous_reg+ alt + post_reg
    else:
        # for deletion need to obtain the new one
        # print(f"ref: {ref}")
        # print(f"region: {region}")
        # print(f"region ref: {ref_region}")
        # print(f"ALT {alt}")
        # print(replace_valid_instance(region, ref, alt))
        df.loc[i, "region_alt"] = replace_valid_instance(region, ref, alt)
  return df
def cpg_var(df):
  for i, row in df.iterrows():
    count_ref = row["region_normal"].count("CG")
    count_alt = row["region_alt"].count("CG")
    if count_ref > count_alt:
        df.loc[i, "CpG_variant"] = "removed"
    elif  count_ref < count_alt:
        df.loc[i, "CpG_variant"] = "added"
    else:
        df.loc[i, "CpG_variant"] = "no"
    df.loc[i, "CpG_site_diff"] = count_alt - count_ref
  return df
        

def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-input_file')           
  parser.add_argument('-output_file', default = None)      
  parser.add_argument('-reference')
  return parser.parse_args()
  
def main(input_file, reference, output_file=None):
    df = get_pre_ref_post_regions(input_file, reference)
    df = get_normal_alt_regions(df)
    df = cpg_var(df)
    if output_file is not None:
        if ".xlsx" in output_file:
            df.to_excel(output_file, index=False)
        else:
            df.to_csv(output_file, index=False)
    return df

if __name__ == "__main__":
    args = get_args()
    main(args.input_file, args.reference, args.output_file)
