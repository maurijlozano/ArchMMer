#!/usr/bin/env python3

#Imports
import os, sys, re, glob
#from io import StringIO 
from Bio import SeqIO
import pandas as pd

#functions
def generate_acc_prot_dict(table):
    unique_proteins = table.Proteins.unique()
    prot_dict = {}
    for prot in unique_proteins:
        acc = prot.split('|')[0]
        if not acc in prot_dict.keys():
            prot_dict[acc] = {}
            prot_dict[acc][prot.split('|')[1]] = prot.split('|')[2]
        else:
            prot_dict[acc][prot.split('|')[1]] = prot.split('|')[2]
    return(prot_dict)

def get_accessions_from_gb_files(folder,file_base_names):
    gb_files = glob.glob(os.path.join(folder,'*.gb'))
    accession_file_dict = {}
    for f in gb_files:
        fn = os.path.splitext(os.path.basename(f))[0]
        if fn in file_base_names:
            for r in SeqIO.parse(f,'gb'):
                accession_file_dict[r.id] = f
    return(accession_file_dict)

#main
if len(sys.argv) != 4:
    sys.exit(f'Usage: {sys.argv[0]} [table <csv output from Archmmer>] [folder <path to folder with genebank files: ext = *.gb>] [expand region < base pairs, int>]')

table_file = sys.argv[1]
folder = sys.argv[2]
expand_bp = int(sys.argv[3])

#get table
table = pd.read_csv(table_file)
#
prot_dict = generate_acc_prot_dict(table)
file_base_names = table.loc[:,'Genome name'].unique()
file_base_names = file_base_names.astype('str')

print(f'--> Getting accession IDs from the gb files in {folder}...')
accession_file_dict = get_accessions_from_gb_files(folder,file_base_names)

clinker_table = pd.DataFrame(columns=['cluster','base_file_name','r-id','start','end','strand','uid'])

for acc in prot_dict:
    file = accession_file_dict[acc]
    print(f'--> Processing {file}...')
    acc_table = pd.DataFrame(columns=['cluster','base_file_name','r-id','start','end','strand','uid'])
    i=0
    for r in SeqIO.parse(file,'gb'):
        for f in r.features:
            if f.type == 'CDS':
                if 'locus_tag' in f.qualifiers:
                    for u,p in prot_dict[acc].items():
                        if p == f.qualifiers['locus_tag'][0]:
                            start = int(f.location.start)
                            end = int(f.location.end)
                            if start > end:
                                cstart = end
                                cend = start
                            else:
                                cstart = start
                                cend = end
                            strand = f.location.strand
                            if strand == 1:
                                cstrand = '+'
                            else:
                                cstrand = '-'
                            info = ['1',os.path.splitext(os.path.basename(file))[0],r.id,cstart,cend,cstrand,u]
                            acc_table.loc[i,:] = info
                            i+=1
                elif 'gene' in f.qualifiers:
                    for u,p in prot_dict[acc].items():
                        if p == f.qualifiers['gene'][0]:
                            start = int(f.location.start)
                            end = int(f.location.end)
                            if start > end:
                                cstart = end
                                cend = start
                            else:
                                cstart = start
                                cend = end
                            strand = f.location.strand
                            if strand == 1:
                                cstrand = '+'
                            else:
                                cstrand = '-'
                            info = ['1',os.path.splitext(file)[0],r.id,cstart,cend,cstrand,u]
                            acc_table.loc[i,:] = info
                            i+=1
    #merge contiguous
    acc_table = acc_table.sort_values(['uid'], ignore_index=True)
    uid0 = acc_table.loc[:,'uid'][0]
    merged_table = pd.DataFrame(columns=acc_table.columns)
    merged_table.loc[0,:] = acc_table.loc[0,:]
    j=0
    for _,row in acc_table.iterrows():
        uid = row.uid
        if not uid == uid0:
            n0 = int(re.sub('UID','',uid0))
            n = int(re.sub('UID','',uid))
            if n in range(n0-2,n0+2):
                merged_table.loc[j,'end'] = row.end
            else:
                j+=1
                merged_table.loc[j,:] = row
            uid0 = uid
    clinker_table = pd.concat([clinker_table,merged_table])

#expand region in bp: expand_bp
clinker_table.start = [int(re.sub('[^0-9]','',str(s))) for s in clinker_table.start]
clinker_table.end = [int(re.sub('[^0-9]','',str(e))) for e in clinker_table.end]
clinker_table.start = [s-expand_bp if s-expand_bp > 0 else 0 for s in clinker_table.start]
clinker_table.end = [e+expand_bp for e in clinker_table.end]
#save table
clinker_table = clinker_table[['cluster','base_file_name','r-id','start','end','strand']]

#Hacer: si se superponen las regiones merge
clinker_table = clinker_table.sort_values(['r-id','start'])
merge_2t = pd.DataFrame(columns = clinker_table.columns)
j=0
if len(clinker_table) > 1:
    for i in range(0,len(clinker_table)-1):
        merge_2t.loc[j,:] = clinker_table.iloc[i,:]
        if (clinker_table.iloc[i,2] == clinker_table.iloc[i+1,2]) and (clinker_table.iloc[i,4] in range(int(clinker_table.iloc[i+1,3]),int(clinker_table.iloc[i+1,4])) ):
            merge_2t.loc[j,'end'] = int(clinker_table.iloc[i+1,4])
        else:
            j+=1
            if i == len(clinker_table)-2:
                merge_2t.loc[j,:] = clinker_table.iloc[i+1,:]
else:
    merge_2t = clinker_table

merge_2t.to_csv('table_for_clinker_plot.csv',header=False, index=False)