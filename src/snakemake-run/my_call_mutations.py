import os
import sys
import re
import pysam
#import fasta_tools as ft
import numpy as np
#np.seterr(divide='ignore', invalid='ignore')
np.set_printoptions(threshold=sys.maxsize,linewidth=sys.maxsize)
import time
from collections import Counter
from tqdm import tqdm

def find_point_muts(coverage,ancestor,min_coverage=1000):
    #is coverage high enough?
    too_low_coverage = sum(coverage,0)<min_coverage
    coverage[:,too_low_coverage] = 0

    #get the fractions of mutations
    total_coverage = sum(coverage,0)
    fractions = 1.0*coverage/total_coverage

    #set the original to 0
    bases = 'ACGT'
    ref_bases = [bases.find(i) for i in ancestor]

#    print(len(fractions[1]))
#    print(len(ref_bases))

    fractions[ref_bases,range(len(ancestor))] = 0
#    print(len(fractions[1]))

    #find the mutations with high enough fractions (we're not interested in low-frequency mutations)
    min_fractions = 0.001
    #return_muts = []
    muts = np.where(fractions>min_fractions)
    muts_to_return = []
    for new, loc in zip(muts[0],muts[1]):
        fraction = fractions[new,loc]
        n_changed = coverage[new,loc]
        ref_base = ancestor[loc]
        coverage_here = total_coverage[loc]
        if coverage_here>min_coverage:
            muts_to_return.append(('P',loc,loc+1,ref_base,bases[new],fraction,n_changed,coverage_here)) # +1 for 0-based indexing format
    return muts_to_return


def parse_cigar(cigar):
    operations = re.findall('([0-9]+)([A-Z]+)', cigar)
    numbers = {'M':0,'I':1,'D':2,'N':3,'S':4,'H':5,'P':6,'=':7,'X':8,'b':9}
    return [(numbers[i[1]],int(i[0])) for i in operations]

def deletion_from_chimera(cigar_begin, cigar_end,begin_pos,end_pos):
    begin = sum([i[1] for i in cigar_begin[:-1] ])
    end = sum([i[1] for i in cigar_end[:-1] ])
    if begin<end:
        if begin_pos+begin < end_pos+end:
            return begin_pos+begin,end_pos
    return None

def find_large_dels(samfile,start=0,stop=None,threshold=5):
    #find size of reference to use here
    if stop is not None:
        size = stop-start
    else:
        size = samfile.lengths[0]

    dels = []

    for read in samfile.fetch(contig=contig_name,start=start,stop=stop):
        if read.mapping_quality>30:
            if 'H' in read.cigarstring or 'S' in read.cigarstring:
                deletion = None
                #already get the real starts and ends of this read (not only where mapping starts and ends)
                start_offset = read.cigar[0][1] if read.cigar[0][0] in [4,5] else 0
                end_offset = read.cigar[-1][1] if read.cigar[-1][0] in [4,5] else 0
                start_read = read.reference_start-start_offset
                end_read = read.reference_end+end_offset
                read_length = sum([i[1] for i in read.cigar])
                #chimeric?
                tags = dict(read.tags)
                if 'SA' in tags.keys():
                    #find chimera
                    fields = tags['SA'].split(',') #the SA tag holds information on chimeric reads
                    if int(fields[4])>30: #if mapping quality of the chimera is high enough
                        chimera_cigar = parse_cigar(fields[3])

                        chimera_start_offset = chimera_cigar[0][1] if chimera_cigar[0][0] in [4,5] else 0
                        chimera_end_offset = chimera_cigar[-1][1] if chimera_cigar[-1][0] in [4,5] else 0

                        chimera_start_ref = int(fields[1]) - 1 # position in 0-based coordinate system
                        chimera_start = chimera_start_ref - chimera_start_offset
                        chimera_end = chimera_start + sum([i[1] for i in chimera_cigar])
                        chimera_end_ref = chimera_end-chimera_end_offset
                        # chimera_len = chimera_end - chimera_start

                        if end_offset!=0 and chimera_start_offset !=0: #original is clipped at end, chimera at beginning
                            if read.reference_end < chimera_start_ref: #the order of the mappings is right
                                if abs(end_offset + chimera_start_offset-read_length) < threshold:
                                    deletion = (read.reference_end,chimera_start_ref)
                                    # end_read = chimera_end #end of the read, for potentially marking later dels
                                    end_offset = chimera_end_offset

                        if start_offset!=0 and chimera_end_offset!=0: #chimera is clipped at end, original at begining
                            if chimera_end_ref < read.reference_start:#order is right
                                if abs(end_offset + chimera_start_offset-read_length) < threshold:
                                    deletion = (chimera_end_ref,read.reference_start)
                                    # start_read = chimera_start
                                    start_offset = chimera_start_offset

                    if deletion is not None:
                        dels.append(deletion)

    counts_dels = Counter(dels)
    dels = Counter()
    for deletion in counts_dels:
        if counts_dels[deletion]>10:
            dels[('LD',deletion[0],deletion[1],"","")]+=1
    return dels


def find_indels(samfile,coverage,contig_name,muts=None,n=None,start=0,stop=None,
                overall_overlap=5,min_coverage=1000):
    if muts is None:
        muts = Counter()
    
    if n is not None:
        counter = 0

    
    
    starts = [i[1] for i in muts if i[0] == 'D']
    ends = [i[2] for i in muts if i[0] == 'D']

    for read in samfile.fetch(contig=contig_name,start=start,stop=stop):
        if read.mapping_quality > 35: # filter on mapping quality was missing in eva's version! # changed from 30 to 35 to be consitent with coverage that is inputed into the function later on
            del_counter = 0
            cigar = read.cigar
            read_location = 0 # to keep everythin 0-based (checked with bam files when set to 1)
            reference_location = read.reference_start
            for operation in cigar:
                if operation[0] == 4:
                    read_location += operation[1]
                    
                if operation[0] == 0:
                    read_location += operation[1]
                    reference_location += operation[1]
                    
                elif operation[0] == 1:
                    insert = read.query_sequence[read_location:read_location+operation[1]]
                    reference_start = read.reference_start
                    muts[('I',reference_start+read_location+del_counter,reference_start+read_location+del_counter,'-',insert)]+=1
                    
                    read_location += operation[1]
                
                elif operation[0] == 2:
                    delstart = reference_location
                    delend = reference_location + operation[1]
                    deletion = ancestor[delstart:delend]
                    muts[ ('D',delstart,delend,deletion,'-')]+=1 
                    
                    reference_location += operation[1]
                    del_counter += operation[1]
                    
                elif operation[0] == 4 or operation[0] == 5: #clip (soft or hard)
                    for offset in range(-overall_overlap,overall_overlap):
                        if read.reference_start+offset in ends:
                            for dele in [i for i in muts if ((i[0] == 'D') and (i[2] == read.reference_start+offset))]:
                                muts[dele]+=1
                        if read.reference_end+offset in starts:
                            for dele in [i for i in muts if ((i[0] == 'D') and (i[1] == read.reference_end+offset))]:
                                muts[dele]+=1
                else:
                    print(operation)
                    raise
            if n is not None:
                counter+=1
                if counter>=n:
                    break

    muts_to_return=[]
    for i in muts:
        pos_start = i[1]
        pos_end = i[2]
        ref_base = i[3]
        coverage_here = max(coverage[pos_start-10:pos_start+10])
        if coverage_here>min_coverage:
            muttype = i[0]
            mutation = i[4]
            counts = muts[i]
            fraction = 1.0*counts/coverage_here
            if fraction > 1:
                fraction = 1
            muts_to_return.append((muttype,pos_start,pos_end,ref_base,mutation,fraction,counts,coverage_here))
    return muts_to_return


if __name__ == '__main__':
    alignment_file = sys.argv[1] #"/home/ali313/Desktop/ltee_raw/test-snakemake/HIV_LTE_mappings/III/13/13MT2EXPIIIVP190seq10102019_S1_L001_aligned.bam"      # sys.argv[1]
    ancestor_file = sys.argv[2] #"/home/ali313/Desktop/ltee_raw/test-snakemake/HIV_LTE_mappings/ancestor/ancestor_consensus.fasta"        #sys.argv[2]
    contig_name = sys.argv[3] #"NL43seq210314" #"NL43_ann_wk0virusPassRef_SEQregion"   #"K03455.1"    #sys.argv[3]
    #extra_info = sys.argv[4]

    fields = alignment_file.split('/')
    sample = fields[-1]
    passage = re.findall(r'VP[0-9]+',sample)[-1][2:]
    try:
        with open('mutations/backup/{}/{}/{}.csv'.format(fields[1],fields[2],passage)) as f:
            csv_here = f.readlines()
        still_to_do = False
    except:
        still_to_do = True

    if not still_to_do:
        for line in csv_here:
            print(line,end='')
    else:
        extra_info = fields[-2]+','+passage
        
        with open(ancestor_file) as f:
            ancestor = f.readlines()

        ancestor = ancestor[1]        
        
        # ancestor = ft.read_fasta(ancestor_file)['consensus']

        start = 0
        stop = len(ancestor)-1 # -1 to not include the new line character at the end!

        ancestor = ancestor[start:stop]
        samfile = pysam.AlignmentFile(alignment_file)

        all_coverage = np.array(samfile.count_coverage(contig_name, start=start,
                                                       stop=stop,quality_threshold=35))

        dels = find_large_dels(samfile)
        muts = find_indels(samfile,np.sum(all_coverage,0),contig_name,dels,min_coverage=100)
        muts+=find_point_muts(all_coverage, ancestor,min_coverage=100)
        for mut in muts:
            for i in mut:
                print(f"{i},",end='')
            print(extra_info)
