import os, argparse, re

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, metavar='</path/to/raw.table>', required=True)
parser.add_argument('-o', type=str, metavar='</path/to/recoded.table>', required=True)
parser.add_argument('-p', type=str, metavar='<name>', required=True)

variant_table, recoded_table, pop = vars(parser.parse_args()).values()

with open(variant_table, 'r') as raw, open(recoded_table, 'w') as recoded:    
    for line in raw:
        if line:
            scaffold, position, reference, an, dp, *raw_gts = line.strip(' \n').split('\t') # extract site info
            raw_gts = [ re.split('\||/', gt) for gt in raw_gts ] # split phased & un-phased genotypes
            ploidy = float(max(len(gt) for gt in raw_gts)) # determine ploidy from genotype
            recoded_gts = [ '-9' if '.' in gt else sum(1 for nucleotide in gt if nucleotide != reference) for gt in raw_gts ] # recode genotypes as alt allele counts
            print(pop, ploidy, scaffold, position, an, dp, *recoded_gts, sep='\t', file=recoded) # write recoded table

