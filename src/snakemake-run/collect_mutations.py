import re
import os
import sys
from tqdm import tqdm
import pickle as pkl
import pandas as pd
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype


def read_files(filename,lines=None,passage=None):
    deletions = []
    insertions = []
    mutations = []

    with open(filename) as f:
        for line in f.readlines():
            fields = line.strip().split(',')
            #mutid = (fields[0],fields[1],fields[2])
            if lines is None:
                lines = fields[6].strip()
            if passage is None:
                passage = int(fields[7])
            try:
                from_pos = int(fields[1])
            except ValueError:
                continue
            except IndexError:
                print(line, lines,passage, 1)
            if fields[0] == 'D': #deletions: start, stop, fraction, reads, coverage at pos
                end_pos = int(fields[2])
                deletions.append([from_pos, end_pos, float(fields[3]),int(fields[4]),int(fields[5]),
                                  lines,passage])
            elif fields[0] == 'I': #insertions: position, insterted bases, fraction, reads, coverage at pos
                inserted = fields[2]
                insertions.append([from_pos, inserted, float(fields[3]),int(fields[4]),int(fields[5]),
                                   lines,passage])
            elif fields[0] == 'M':
                new_base = fields[2]
                mutations.append([from_pos,new_base,float(fields[3]),int(fields[4]),int(fields[5]),lines,
                                  passage])
    try:
        if mutations[-1][0] < 9000:
            print(filename, mutations[-1][0], 2)
    except IndexError:
        print(lines,filename[filename.find('VP')+2:filename.find('.')],3)
    return deletions,insertions,mutations

'''
line_translation = {'MT21':'MT21',
                    'MT22':'MT22',
                    'MT41':'MT41',
                    'MT42':'MT42',
                    'MT21_low':'MT21',
                    'MT22_low':'MT22',
                    'MT41_low':'MT41',
                    'MT42_low':'MT42',
                    '13MT2':'MT23',
                    '14MT2':'MT24',
                    '15MT4':'MT43',
                    '16MT4':'MT44',
                    'MT21MT2':'MT21MT2',
                    'MT22MT2':'MT22MT2',
                    'MT41MT2':'MT41MT2',
                    'MT42MT2':'MT42MT2',
                    'MT21MT4':'MT21MT4',
                    'MT22MT4':'MT22MT4',
                    'MT41MT4':'MT41MT4',
                    'MT42MT4':'MT42MT4',
                    '17MT2':'MT25',
                    '18MT2':'MT26',
                    '19MT2':'MT27',
                    '20MT2':'MT28',
                    'MT25':'MT25',
                    'MT26':'MT26',
                    'MT27':'MT27',
                    'MT28':'MT28'}

my_line_translation =  {'13':'MT21',
                        '14':'MT22',
                        '15':'MT41',
                        '16':'MT42',
                        '13R':'MT21_R',
                        '14R':'MT22_R',
                        '15R':'MT41_R',
                        '16R':'MT42_R'}

'''


dels,ins,muts = read_files('/home/amovas/data/genome-evo-proj/data/ancestor/ancestor_consensus.csv','ancestor',0)

for line in sys.argv[1:]:
#    line_name = my_line_translation[line.split('/')[2]]
    line_name = line.split('/')[9]
    time = int(re.findall(r'VP([0-9]+)',line)[-1])
    print(time,line_name)
    d,i,m = read_files(line,line_name,time)
    dels += d
    ins += i
    muts += m


print("ERROR NOT YET!")


# #ancestor
# time = 0
# for line in tqdm(['MT21', 'MT22', 'MT42', 'MT41',
#                            'MT21MT2','MT21MT4', 'MT22MT2', 'MT22MT4',
#                            'MT41MT2','MT41MT4','MT42MT4','MT42MT2',
#                            '13MT2','14MT2','15MT4','16MT4']):
#     if line[0] == '1':
#         line_name = line[2:]+(line[1] if (int(line[1])<5) else str((int(line[1])-2)))
#     else:
#         line_name = line
#     dels,ins,muts = read_files('../mutations/ancestor/ancestor.csv',line_name,time)
#
# #individual lines
# for line in tqdm(['MT21', 'MT22', 'MT42', 'MT41']):
#     for time in range(10,301,10):
#         fname = '/home/eva/Research/HIV_LTE_DATA/mutations/individual_lines/{}_low/VP{}.csv'.format(line,time)
#         if time < 60:
#             other_lines = ['','MT2','MT4']
#         else:
#             other_lines = ['']
#         for line2 in other_lines:
#             (d,i,m) = read_files(fname,line+line2,time)
#             dels += d
#             ins += i
#             muts += m
#
# #crosses
# for line in tqdm(['MT21MT2','MT21MT4', 'MT22MT2', 'MT22MT4','MT41MT2','MT41MT4','MT42MT4','MT42MT2']):
#     for time in range(60,271,10):
#         fname = '/home/eva/Research/HIV_LTE_DATA/mutations/cross/{}/VP{}.csv'.format(line,time)
#         (d,i,m) = read_files(fname,line,time)
#         dels += d
#         ins += i
#         muts += m
#
# #exp 3
# for line in tqdm(['13MT2','14MT2','15MT4','16MT4']):
#     for time in range(10,181,10):
#         fname = '/home/eva/Research/HIV_LTE_DATA/mutations/expIII/{}/VP{}.csv'.format(line,time)
#         line_name = line[2:]+(line[1] if (int(line[1])<5) else str((int(line[1])-2)))
#         (d,i,m) = read_files(fname,line_name,time)
#         dels += d
#         ins += i
#         muts += m
#
# #exp4
# for line in tqdm(['MT25','MT26','MT27','MT28']):
#     for time in range(10,181,10):
#         fname = '/home/eva/Research/HIV_LTE_DATA/mutations/expIV/{}/VP{}.csv'.format(line,time)
#         line_name =line
#         (d,i,m) = read_files(fname,line_name,time)
#         dels += d
#         ins += i
#         muts += m
#

all_dels = pd.DataFrame(dels,columns=['start','end','fraction','reads','coverage','line','passage'])
all_ins = pd.DataFrame(ins,columns=['pos','inserted','fraction','reads','coverage','line','passage'])
all_muts = pd.DataFrame(muts,columns=['pos','new','fraction','reads','coverage','line','passage'])


all_dels.sort_values(by='fraction',ascending=False,inplace=True)
all_ins.sort_values(by='fraction',ascending=False,inplace=True)
all_muts.sort_values(by='fraction',ascending=False,inplace=True)

# all_dels.at[all_dels.fraction>1,'fraction']=1
all_dels.to_pickle('/home/amovas/data/genome-evo-proj/results/tables/2_p/all_dels.pkl',protocol=2)
all_muts.to_pickle('/home/amovas/data/genome-evo-proj/results/tables/2_p/all_muts.pkl',protocol=2)
all_ins.to_pickle('/home/amovas/data/genome-evo-proj/results/tables/2_p/all_ins.pkl',protocol=2)


# convert the mutation pkl object to csv for further analyses
with open("/home/amovas/data/genome-evo-proj/results/tables/2_p/all_muts.pkl", "rb") as f:
    object = pkl.load(f)
    
df = pd.DataFrame(object)
df.to_csv(r'/home/amovas/data/genome-evo-proj/results/tables/2_p/all_muts.csv')



# convert the insertion pkl object to csv for further analyses
with open("/home/amovas/data/genome-evo-proj/results/tables/2_p/all_ins.pkl", "rb") as f:
    object = pkl.load(f)
    
df = pd.DataFrame(object)
df.to_csv(r'/home/amovas/data/genome-evo-proj/results/tables/2_p/all_ins.csv')



# convert the deletion pkl object to csv for further analyses
with open("/home/amovas/data/genome-evo-proj/results/tables/2_p/all_dels.pkl", "rb") as f:
    object = pkl.load(f)
    
df = pd.DataFrame(object)
df.to_csv(r'/home/amovas/data/genome-evo-proj/results/tables/2_p/all_dels.csv')
