#!/usr/bin/env python3
#imports
import os, sys, glob
import pandas as pd
from Bio import SeqIO

#from table with regions to align make a clinker plot
print('This script requires Clinker installed and accessible...')
print(f'Checking arguments...')
if len(sys.argv) < 3:
    sys.exit(f'Usage: {sys.argv[0]} [table <csv file with cluster file_name record_id start end strand(+,-) >] [cluster_index <int>] [gb file folder <complete path>]')

table_file = sys.argv[1]
cluster_index = int(sys.argv[2])-1
if cluster_index < 0:
    sys.exit('Cluster index is 1 indexed. The first cluster should be numbered "1"...')
folder = sys.argv[3]
if not os.path.exists(folder):
    sys.exit(f'{folder} not found...')
else:
    gb_files = glob.glob(os.path.join(folder,'*.gb'))
    if len(gb_files) == 0:
        sys.exit('No gb files found...')

output_folder = 'regions'
if not os.path.exists(output_folder):
    print(f'Creating output folder: regions...')
    os.mkdir(output_folder)

table = pd.read_csv(table_file, header=None)

if len(table.columns) != 6:
    print('Error, incorrect table format...')
    sys.exit(f'Usage: {sys.argv[0]} [table <csv file with cluster file_name record_id start end strand(+,-)>] [cluster_index <int>] [gb file folder <complete path>]')

cluster_name = table.iloc[:,0].unique()[cluster_index]
ctable = table[~table[2].isna()]
ctable = ctable[ctable[0] == cluster_name]

#generate table file names. Extension = gb
ctable.loc[:,'files'] = [os.path.join(folder, f'{f}.gb') for f in ctable[1] ]
bigger_region = max([ abs(r[4]-r[3])  for _,r in ctable.iterrows() ]) * 2
#extract regions for clinker
print('Extracting regions from gb files...')
for _,reg in ctable.iterrows():
    file = reg.files
    name = reg[1]
    srecid = reg[2]
    start = int(reg[3])
    end = int(reg[4])
    strand = reg[5]
    expand_region_len = int((bigger_region - abs(end-start)) /2)
    if (start-expand_region_len) > 0:
        cstart = start-expand_region_len
    else:
        cstart = 0
    cend = end + expand_region_len
    if not os.path.exists(file):
        print(f'File {file} not found, skipping...')
    else:
        region_file_name = os.path.join(output_folder,f'{name}_{srecid}_{start}_{end}.gb')
        if not os.path.exists(region_file_name):
            for r in SeqIO.parse(file,'gb'):
                if r.id == srecid:
                    region_rec = r[cstart:cend]
                    if strand == '+': 
                        SeqIO.write(region_rec,region_file_name,'gb')
                    else:
                        region_rec_rc = region_rec.reverse_complement()
                        region_rec_rc.id = r.id
                        region_rec_rc.annotations = r.annotations
                        region_rec_rc.description = r.description
                        region_rec_rc.name = r.name
                        SeqIO.write(region_rec_rc,region_file_name,'gb')
                    break
        else:
            print(f'{region_file_name} found, skipping...')

#running clinker
print('> Running clinker...')
clinker_html_out = f'clinker_plot_cluster_{cluster_index}.html'
os.system(f'clinker {os.path.join(output_folder,"*.gb")} -p {clinker_html_out}')
print('Task done!')

