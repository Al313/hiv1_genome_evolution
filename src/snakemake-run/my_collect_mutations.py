# loading packages

import re
import os
import sys
from tqdm import tqdm
import pickle as pkl
import pandas as pd
import gzip
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype



# define functions 

def flatten(l):
    return [item for sublist in l for item in sublist]

def read_files(filename,lines=None,passage=None):
    deletions = []
    insertions = []
    mutations = []

    with open(filename) as f:
        for line in f.readlines():
            fields = line.strip().split(',')
            #mutid = (fields[0],fields[1],fields[2])
            if lines is None:
                lines = fields[8].strip()
            if passage is None:
                passage = int(fields[9])
            try:
                from_pos = int(fields[1])
                end_pos = int(fields[2])
                ref_base = str(fields[3])
                alt_base = str(fields[4])
                
            except ValueError:
                continue
            except IndexError:
                print(line, lines,passage, 1)
            if fields[0] == 'D' or fields[0] == 'LD': #deletions: start, stop, fraction, reads, coverage at pos
                
                deletions.append([from_pos, end_pos, ref_base, alt_base, fields[0], float(fields[5]),int(fields[6]),int(fields[7]),
                                  lines,passage])
            elif fields[0] == 'I': #insertions: position, insterted bases, fraction, reads, coverage at pos
                
                insertions.append([from_pos, end_pos, ref_base, alt_base, fields[0], float(fields[5]),int(fields[6]),int(fields[7]),
                                   lines,passage])
            elif fields[0] == 'P':
                
                mutations.append([from_pos, end_pos, ref_base, alt_base, fields[0], float(fields[5]),int(fields[6]),int(fields[7]),
                                  lines,passage])
    try:
        if mutations[-1][0] < 9000:
            print(filename, mutations[-1][0], 2)
    except IndexError:
        print(lines,filename[filename.find('VP')+2:filename.find('.')],3)
    return deletions,insertions,mutations


if __name__ == "__main__":
    
    variants = []
    
    
    for sample_path in sys.argv[1:]:
        line_nu = sample_path.split('/')[-1][:2]
        passage = int(re.findall(r'VP([0-9]+)',sample_path)[-1])
        
        d,i,m = read_files(sample_path,line_nu,passage)
        
        variants += d,i,m
        
    variants = pd.DataFrame(flatten(variants),columns=['start','end','ref','alt','mut_type','fraction','reads','coverage','line','passage'])
    variants.sort_values(by='fraction',ascending=False,inplace=True)
    variants.to_pickle('/home/amovas/data/genome-evo-proj/results/tables/2-p/all_variants.pkl.gz',protocol=2, compression='gzip')

    # convert the pkl object to csv for further downstream analyses
    with gzip.open("/home/amovas/data/genome-evo-proj/results/tables/2-p/all_variants.pkl.gz", "rb") as f:
        object = pkl.load(f)
        
    df = pd.DataFrame(object)
    df.to_csv("/home/amovas/data/genome-evo-proj/results/tables/2-p/all_variants.csv.gz", index=False, compression='gzip')
