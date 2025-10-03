#!/usr/bin/env python3
'''This program was designed to identify and classify proteins with a defined architecture on bacterial genomes.'''
NAME = "ArcHMMer"
VERSION = "1.0"
REF = 'Not published'
GITHUB="https://github.com/maurijlozano"

#Imports
import argparse, requests, gzip, os, sys, re, glob, subprocess
#from io import StringIO 
from Bio import SeqIO
import pandas as pd

#functions
def parseArgs():
    '''
    Argument parsing is done.
    '''
    parser = argparse.ArgumentParser(description='This program was designed to identify and classify proteins with a defined architecture in bacteria.')
    # genome files
    inputArgs = parser.add_argument_group('Input files')
    inputArgs.add_argument("-i", "--InputFolder",help="Folder containing the genomes to analyze in Genbank format. The extension of the genbank files should be passed by -x argument.", dest="inputFolder", action='store', required=False)
    inputArgs.add_argument("--Overwrite",help="Overwrite the results folder.", dest="overwrite", action='store_true', required=False)
    inputArgs.add_argument("-x", "--GbFileExtension",help="Genbank file extension: gb, gbff, gbk.", dest="GbExtension", action='store', required=False)
    downloadOpt = parser.add_argument_group('Download genomes from NCBI')
    downloadOpt.add_argument("-G", "--Genus",help="Query genome genus for NCBI search.",action="store", dest="queryGenus", required=False)
    downloadOpt.add_argument("-n", "--numberOfSequences",help="Number of sequences to download.",action="store", dest="numSeqs", required=False)
    downloadOpt.add_argument("-f", "--DownloadFolder",help="Folder that will be used to download NCBI genomes and as working folder.", dest="downloadFolder", action='store', required=False)
    downloadOpt.add_argument("-a", "--accession",help="List of NCBI assembly accession numbers.", dest="accessions", nargs = '+', action='store', default = [], required=False)
    #reference genome. Optional
    refgenome = parser.add_argument_group('Reference genome')
    refgenome.add_argument("-r", "--RefGenome",help="Genome to use as reference for plot. The Gb file name must be passed as argument (e.g., 'genome.gb')", dest="refGenome", action='store', default='', required=False)
    refgenome.add_argument("-R", "--RefAccession",help="Genome to use as reference for plot. Genbank assembly accession number", dest="refacc", action='store', default='', required=False)
    #HMMER options
    hmmerArgs = parser.add_argument_group('HMMER options')
    hmmerArgs.add_argument("-c", "--ncpu",help="Number of CPUs to use in HMMER.", dest="ncpu", action='store', default='auto', required=False)
    hmmerArgs.add_argument("-e", "--evalue",help="E-value for hmmsearch.", dest="evalue", action='store', default=1e-5, required=False)
    hmmerArgs.add_argument("-C", "--hmm_coverage",help="HMM coverage for hmmsearch.", dest="hmm_coverage", action='store', default=70, required=False)
    #Architectures and HMM files
    archArgs = parser.add_argument_group('Domain architecture definition') 
    archArgs.add_argument("-p", "--hmm_ids",help="Architecture of the protein domains. A list of the PFAM/Superfamily/Panther/TIGRfam/CATH IDs from InterPro in the correct order and number, representing the target architecture. Eg. -p PFAM_A PFAM_B PFAM_B PFAM_C", dest="pfamid", nargs = '+', action='store', default = [], required=False)
    archArgs.add_argument("-P", "--hmm_files",help="Architecture of the protein domains. A list of the HMM files in the correct order and number, representing the target architecture. Eg. -p file_A.hmm file_B.hmm file_B.hmm file_C.hmm", dest="hmms", nargs = '+', action='store', default = [], required=False)
    archArgs.add_argument("-s", "--min_size",help="Minimum protein size for the selected architecture. <int>", dest="min_size", action='store', default = 0, required=False)
    archArgs.add_argument("-S", "--max_size",help="Maximum protein size for the selected architecture. <int>", dest="max_size", action='store', default = 50000, required=False)
    archArgsExl = parser.add_argument_group('Proteins with the specified HMMs will be excluded from the results') 
    archArgsExl.add_argument("--Exclude_hmm",help="List of HMM profiles. The proteins with a significant match to this profiles will be excluded. Eg. -E PFAM_A PFAM_B PFAM_B PFAM_C", dest="exclude_pfamids", nargs = '+', action='store', default = [], required=False)
    archArgsExl.add_argument("--Exclude_hmm_table",help="Table of HMM profiles. The proteins with a significant match to this profiles will be excluded. Path of a table with one PFAM ID per row. <table.txt>", dest="exclude_pfamids_table", action='store', required=False)
    archArgsExl.add_argument("--Exclude_hmm_files",help="List of HMM profile files. The proteins with a significant match to this profiles will be excluded. Eg. -E file_A.hmm file_B.hmm file_B.hmm file_C.hmm", dest="exclude_hmms", nargs = '+', action='store', default = [], required=False)
    #architecture files
    archArgsFiles = parser.add_argument_group('Specify the protein architecture in a txt file') 
    archArgsFiles.add_argument("--Architecture_file",help="A file with one architecture of the protein domains per line. Each architecture is defined by a list of HMM files in the correct order and number (See example file). It can also include the hmm files of domains to exclude from the architecture (space separated list) and hmm files of domains to include in the syntenic region (space separated). All these lists must be in one line and separated by a comma. The last two comma separated entries correspond to minimum and maximum protein size.", dest="architecture_file", action='store', default = [], required=False)
    archArgsFiles.add_argument("--Architecture_hmmid",help="A file with one architecture of the protein domains per line. Each architecture is defined by a list of PFAM/Superfamily/Panther/TIGRfam/CATH IDs from InterPro in the correct order and number (See example file). It can also include the PFAM ids of domains to exclude from the architecture (space separated list) and PFAM ids of domains to include in the syntenic region (space separated). All these lists must be in one line and separated by a comma. The last two comma separated entries correspond to minimum and maximum protein size.", dest="architecture_pfamid", action='store', default = [], required=False)
    #synteny
    syntenyArgs = parser.add_argument_group('Domain conservation in the genomic context')
    syntenyArgs.add_argument("--include_in_synteny_hmm",help="List of PFAM IDs for HMM profiles that should be found in the regions adjacent to the proteins of interest.", dest="include_in_synteny", nargs = '+', action='store', default = [], required=False)
    syntenyArgs.add_argument("--include_in_synteny_hmm_files",help="List of HMM profile files to look for in the regions adjacent to the proteins of interest.", dest="include_hmms", nargs = '+', action='store', default = [], required=False)
    syntenyArgs.add_argument("--force_synteny",help="Only include results with correct synteny.", dest="force_synteny", action='store_true',required=False)
    # Search parameters
    syntenyArgs.add_argument("-u", "--upstreamGenes",help="Number of upstream genes to analyze.", dest="up", action='store', default = 2, required=False)
    syntenyArgs.add_argument("-d", "--downstreamGenes",help="Number of downstream genes to analyze.", dest="down", action='store', default = 2, required=False)
    syntenyArgs.add_argument("-t", "--tolerance",help="The number of genes to look ahead for contiguity, if there is expected to detect tandem copies of the protein. Default =2, <int> .", dest="tolerance", action='store', default=2, required=False)
    syntenyArgs.add_argument("-m", "--more_copies",help="Filters out proteins with the desired architecture, but present in only one copy in a genomic region.", dest="more_than_one_contiguous", action='store_true', required=False)
    #
    args = parser.parse_args()
    return args

def printSoftName():
    print("\n\n\n")
    print("   ******************************************************************************")
    print("   *****   ArcHMMer: Search for proteins with defined domain architecture")
    print("   *****   Version: "+str(VERSION))
    print("   *****   Developed by Mauricio J. Lozano")
    print("   *****   github.com/maurijlozano")
    print("   ******************************************************************************")
    print("   Please cite: "+REF)
    print("   Downloaded from: "+GITHUB)
    print("\n\n\n")

#Download genomes from NCBI

def downLoadNCBIftp(ftpLinkGB,inputFolder,fileName):
    httpsLink = re.sub("ftp://","https://",ftpLinkGB)
    gzfileName = os.path.join(inputFolder,os.path.basename(httpsLink))
    try:
        with open(gzfileName, 'wb') as f:
            f.write(requests.get(httpsLink).content)
        with gzip.open(gzfileName, 'rb') as f:
            with open(fileName,'wb') as f_out:
                f_out.write(f.read())
        os.remove(gzfileName)
    except:
        if os.path.exists(gzfileName):
            os.remove(gzfileName)
        print("An error occurred while downloading the genome file from NCBI...")
    return

def getSeqsFromNCBI(queryGenus,inputFolder,numSeqs):
    #Retrieve fasta files for genomes of a selected genus...
    genomeIDs = []
    #NCBItools
    seqcount = 0
    Id = '"'+queryGenus+'"'
    r = requests.get(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term="{Id}"[Organism] AND (latest[filter] AND "complete genome"[filter] AND "refseq has annotation"[Properties])&retmode=json&retmax=1000&idtype=acc', headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        print(f'-> {queryGenus} not found...')
        return(genomeIDs)
    text = r.json()
    genomeIDList = text['esearchresult']['idlist']
    if (numSeqs > 0) and (len(genomeIDList) > numSeqs):
        genomeIDList = genomeIDList[0:numSeqs]
    nstrains = len(genomeIDList)
    print("-> A total of " + str(nstrains) + " " + Id + " strains where found.")
    if nstrains != 0:
        for id in genomeIDList:
            url =f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={id}&retmode=json'
            try:
                r = requests.get(url, headers={ "Content-Type" : "application/json"})
            except:
                return(genomeIDs)    
            if not r.ok:
                r.raise_for_status()
                print(f'-> {id} not found...')
                return(genomeIDs)
            text = r.json()
            try:
                ftpLink = text['result'][id]['ftppath_refseq']
                if ftpLink != "":
                    ftpLinkGB=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.gbff.gz"
                else:
                    ftpLink = text['result'][id]['ftppath_genbank']
                    if ftpLink != "":
                        ftpLinkGB=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.gbff.gz"
            except:
                print("-> Unable to find sequence in genbank format...")
                continue
            print(inputFolder,id+".gb")
            fileName = os.path.join(inputFolder,id+".gb")
            if not os.path.isfile(fileName):
                # download ftp..
                print("-> Downloading "+id+" sequence in genbank format...")
                downLoadNCBIftp(ftpLinkGB,inputFolder,fileName)
            else:
                print(f'-> {id} already downloaded...')
            genomeIDs.append(id)
        return(genomeIDs)
    else:
        print(f"-> Unable to download {queryGenus} sequences...")
        return(genomeIDs)
        #sys.exit()

def get_Genome_by_assembly_accession(downloadid,inputFolder):
    try:
        r = requests.get(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term="{downloadid}"&retmode=json&retmax=1000&idtype=acc', headers={ "Content-Type" : "application/json"})
        if not r.ok:
            r.raise_for_status()
            print(f'-> {downloadid} not found...')
            return
        text = r.json()
        id = text['esearchresult']['idlist'][0]
    except:
        print(f'-> {downloadid} not found...')
        return
    url =f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={id}&retmode=json'
    try:
        r = requests.get(url, headers={ "Content-Type" : "application/json"})
    except:
        return   
    if not r.ok:
        r.raise_for_status()
        print(f'-> {id} not found...')
        return
    text = r.json()
    try:
        ftppath = text['result'][id]['ftppath_refseq']
    except:
        return
    ftplink = ftppath+'/'+os.path.basename(ftppath)+"_genomic.gbff.gz"
    genomeID = text['result'][id]['assemblyaccession']
    fileName = os.path.join(inputFolder,genomeID+'.gb')
    if os.path.exists(fileName):
        return(fileName)
    downLoadNCBIftp(ftplink,inputFolder,fileName)
    return

def test_file_list(list):
    for file in list:
        if len(file) > 0:
            if not os.path.exists(file):
                sys.exit(f'\nError: {file} not found...')

#Download PFAM HMM profiles
def download_pfam_hmm(pfamid, pfam_file):
    if (os.path.exists(pfam_file)) and (os.path.getsize(pfam_file) > 1024):
        print(f'--> {pfamid} already downloaded...')
        return
    print(f'--> Downloading {pfamid}...')
    urls = {'PFAM':f'https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{pfamid}?annotation=hmm', 
                               'PANTHER':f'https://www.ebi.ac.uk/interpro/wwwapi//entry/panther/{pfamid}?annotation=hmm',
                               'CATH':f'https://www.ebi.ac.uk/interpro/wwwapi//entry/cathgene3d/{pfamid}?annotation=hmm',
                               'SUPERFAMILY': f'https://www.ebi.ac.uk/interpro/wwwapi//entry/ssf/{pfamid}?annotation=hmm',
                               'TIGRFam':f'https://www.ebi.ac.uk/interpro/wwwapi//entry/ncbifam/{pfamid}?annotation=hmm'
            }
    for url in urls:
        try:
            r = requests.get(urls[url], headers={'Accept-Encoding': 'gzip'}, stream=True)
            if not r.ok:
                #r.raise_for_status()
                #sys.exit(f'Pfam id: {pfamid} not found.')
                continue
            else:
                hmm = str(gzip.decompress(r.content), 'utf-8')
                with open(pfam_file,'w+') as hf:
                    hf.write(hmm)
                success = True
        except:
            success = False
    if not success:
        sys.exit(f'\nERROR: {pfamid} couldn\'t be downloaded.')
    return

def download_architecture_pfams(pfamids,type):
    hmms = []
    if len(pfamids) != 0:
        if type == 'arch':
            pretty_architecture = '-'.join([pfi for pfi in pfamids])
            print(f'-> Downloading HMM profiles for the architecture {pretty_architecture}')
        for pfamid in pfamids:
            pfamid = pfamid.upper()
            if len(pfamid) > 0:
                pfam_file = os.path.join(pfam_hmm_folder, pfamid+'.hmm')
                if not pfam_file in hmms:
                    download_pfam_hmm(pfamid, pfam_file)
                    hmms.append(pfam_file)
        test_file_list(hmms)
        return(hmms)
    else:
        pass

#Genbankfile to dict
def genbankToDict(genbankFile):
	genbankDict = {}
	description = {}
	for record in SeqIO.parse(genbankFile, "genbank"):
		featuresDict = {}
		count=1
		if record.features:
			for feat in record.features:
				if feat.type == "CDS":
					if not 'pseudo' in feat.qualifiers:
						UID = "UID"+str(count)
						proteinSeq = feat.qualifiers['translation'][0]
						if "locus_tag" in feat.qualifiers.keys():
							locusTag = feat.qualifiers["locus_tag"][0]
						elif "gene" in feat.qualifiers.keys():
							locusTag = feat.qualifiers["gene"][0]
						else:
							locusTag = "ND"
						if 'protein_id' in feat.qualifiers:
							proteinID = feat.qualifiers['protein_id'][0]
						else:
							proteinID = "ND"
						featuresDict[UID] = [locusTag,proteinID,proteinSeq,int(feat.location.start), int(feat.location.end), int(feat.location.strand), UID]
						count+=1
				if feat.type == 'source':
					if 'strain' in feat.qualifiers.keys():
						strain = feat.qualifiers['strain'][0]
						organism = feat.qualifiers['organism'][0] + ' (' + strain +')'
					elif 'isolate' in feat.qualifiers.keys():
						strain = feat.qualifiers['isolate'][0]
						organism = feat.qualifiers['organism'][0] + ' (' + strain +')'
					else:
						organism = feat.qualifiers['organism'][0]
		description[record.id] = organism
		genbankDict[record.id] = featuresDict
	return(genbankDict,description)

#Extract translation for all proteins (proteome) from genbank file
def extractProteome(genbankDict, genome_name, genome_folder):
    proteomeFile = os.path.join(genome_folder, genome_name + ".faa")
    with open(proteomeFile, "w+") as f:
        for key,value in genbankDict.items():
            for k,val in genbankDict[key].items():
                f.write(">"+key+"|"+val[6]+"|"+val[0]+"\n"+str(val[2])+"\n")
    return(proteomeFile)

def hmmsearch(proteomeFile,hmm, genome_folder, evalue, hmm_coverage):
    hmm_name = os.path.splitext(os.path.basename(hmm))[0]
    resFile = os.path.join(genome_folder, hmm_name + ".txt")
    subprocess.call(['hmmsearch', "--domtblout", resFile , hmm , proteomeFile], stdout=subprocess.DEVNULL)
    tableHeaders = ['target name', 'accession', 'tlen', 'query name', 'qaccession', 'qlen', 'E-value', 'score', 'bias', '#', 'of', 'c-Evalue', 'i-Evalue', 'dom score', 'dom bias', 'hmm from', 'hmm to', 'ali from', 'ali to', 'env from', 'env to', 'acc', 'description of target']
    try:
        hmmHits = pd.read_fwf(resFile, comment='#', sep='\\s+', header=None)
    except:
        if 'exclude' in hmm_name:
            print(f'---> No domains in the exclude list were found...')
        elif 'include' in hmm_name:
            print(f'',end='')
        else:
            print(f'---> No homologs of {hmm_name} were found...')
        return(pd.DataFrame(columns = tableHeaders),[])
    if len(hmmHits) > 0:
        renameDict = { old:new for old,new in zip(hmmHits.columns[0:len(tableHeaders)],tableHeaders)}
        hmmHits = hmmHits.rename(columns=renameDict)
    hmmHits = hmmHits[hmmHits['i-Evalue'] < evalue]
    #coverage
    hmmHits['hmm_coverage'] = (hmmHits['hmm to'] - hmmHits['hmm from'])/(hmmHits['qlen'])*100
    hmmHits = hmmHits[hmmHits['hmm_coverage'] > hmm_coverage]
    proteins = list(hmmHits['target name'].unique())
    return(hmmHits,proteins)

def extract_protein_subset(genbankDict, protein_subset):
    proteomeFile = os.path.join(genome_folder, "protein_subset.faa")
    with open(proteomeFile, "w+") as f:
        ps = {i.split('|')[1]:i.split('|')[0] for i in protein_subset}
        for prot in ps.keys():
            key = ps[prot]
            val = genbankDict[key][prot]
            f.write(">"+key+"|"+val[6]+"|"+val[0]+"\n"+str(val[2])+"\n")
        return(proteomeFile)

def get_domain_names(hmms):
    domains = []
    for hmm in hmms:
        with open(hmm,'r') as f:
            for line in f:
                if re.match('^NAME',line):
                    domain_name = re.search('^NAME +(.*)$',line).groups()[0]
                    domains.append(domain_name)
    return(domains)

def load_genome(genome_file):
    try:
        genbankDict, description = genbankToDict(genome_file)
        proteomeFile = extractProteome(genbankDict, genome_name, genome_folder)
        return(genbankDict, description, proteomeFile)
    except:
        #sys.exit("An error ocurred. Is the genome in genbank format?")
        print("An error ocurred. Is the genome in genbank format? Skipping {genome_file}...")
        return(None,None,None)

def search_by_architecture(genome_file, genome_name, hmms, genome_folder, evalue, hmm_coverage, exclude_file, genbankDict, proteomeFile, min_size, max_size):
    already_looked_for = []
    all_hits_table = pd.DataFrame()
    for hmm in hmms:
        print(f'--> HMM profile: {hmm} ...')
        try:
            if not hmm in already_looked_for:
                hmmHits, protein_subset = hmmsearch(proteomeFile, hmm, genome_folder, evalue, hmm_coverage)
                hmmHits = hmmHits[(min_size <= hmmHits['tlen']) & (hmmHits['tlen'] <= max_size )]
                print(f'---> A total of {len(hmmHits)} were found for {hmm}.')
                if len(hmmHits) == 0:
                    print(f'--> {hmm} not found on {genome_name}, skipping...')
                    return(pd.DataFrame(columns=['target name','tlen']), pd.DataFrame(columns=['target name', 'accession', 'tlen', 'query name', 'qaccession', 'qlen', 'E-value', 'score', 'bias', '#', 'of', 'c-Evalue', 'i-Evalue', 'dom score', 'dom bias', 'hmm from', 'hmm to', 'ali from', 'ali to', 'env from', 'env to', 'acc', 'description of target']))
                all_hits_table = pd.concat([all_hits_table,hmmHits])
                proteomeFile = extract_protein_subset(genbankDict, protein_subset)
                already_looked_for.append(hmm)
        except:
            already_looked_for.append(hmm)
            #sys.exit(f'An error occurred while searching for {hmm}')
            return(pd.DataFrame(columns=['target name','tlen']), pd.DataFrame(columns=['target name', 'accession', 'tlen', 'query name', 'qaccession', 'qlen', 'E-value', 'score', 'bias', '#', 'of', 'c-Evalue', 'i-Evalue', 'dom score', 'dom bias', 'hmm from', 'hmm to', 'ali from', 'ali to', 'env from', 'env to', 'acc', 'description of target']))
    #
    # #check architecture
    domains = get_domain_names(hmms)
    if len(all_hits_table) > 0:
        all_hits_table = all_hits_table[all_hits_table['target name'].isin(protein_subset)]
    else:
        print(f'--> Skipping {genome_file}, no proteins with the required domains were found...')
        return('', '')
    verified_targets = []
    for target in all_hits_table['target name'].unique():
        target_hits = all_hits_table[all_hits_table['target name'] == target].sort_values(['ali from'])
        found_domains = list(target_hits['query name'])
        if domains == found_domains:
            verified_targets.append(target)
    all_hits_table = all_hits_table[all_hits_table['target name'].isin(verified_targets)]
    #eliminate proteins with additional domains not present in flagellins
    if len(exclude_file) > 0:
        print(f'-> Looking for domains in the exclude list...')
        exclude_hmmHits, exclude_proteins = hmmsearch(proteomeFile, exclude_file, genome_folder, evalue, hmm_coverage)
        all_hits_table = all_hits_table[~all_hits_table['target name'].isin(exclude_proteins)]
    #proteins = list(all_hits_table['target name'].unique())
    proteins = all_hits_table.groupby(['target name'], as_index=False).first()[['target name','tlen']]
    return(proteins, all_hits_table)

def evaluate_contiguity(proteins, tolerance):
    protein_dict = {}
    for prot in proteins:
        rep = prot.split('|')[0]
        val = prot.split('|')[1]
        if rep in protein_dict.keys():
            protein_dict[rep].append(val)
        else:
            protein_dict[rep] = [val]
    #
    rep_dict = {}
    for k in protein_dict.keys():
        contiguity_cluster = 0         
        cluster_dict = {}
        cluster_dict[contiguity_cluster] = []
        rep_dict[k] = cluster_dict
        rep_prots = protein_dict[k]
        rep_prots = [ int(re.sub('UID','',id)) for id in rep_prots ]
        rep_prots.sort()
        nprots = len(rep_prots)
        cluster_dict[contiguity_cluster].append(rep_prots[0])
        for i,n in enumerate(rep_prots):
            j = i+1
            if j < nprots:
                if rep_prots[j] in range(n+1,n+1+tolerance):
                    if not rep_prots[j] in cluster_dict[contiguity_cluster]:
                        cluster_dict[contiguity_cluster].append(rep_prots[j])
                else:
                    contiguity_cluster += 1
                    cluster_dict[contiguity_cluster] = []
                    cluster_dict[contiguity_cluster].append(rep_prots[j])
    #
    nreps = len(rep_dict)
    if nreps == 1:
        for k in rep_dict:
            cluster_dict_k = rep_dict[k]
            nclusters = len(cluster_dict_k)
            if nclusters == 1:
                elements = list(cluster_dict_k.values())[0]
                nprots = len(elements)
                return('True', [(k, elements[0]-1, elements[-1]+1,nprots)])
            else:
                up_down = []
                for c in cluster_dict_k:
                    prots_in_c = cluster_dict_k[c]
                    nprots = len(prots_in_c)
                    up_down.append((k, prots_in_c[0]-1,prots_in_c[-1]+1,nprots))
                return('False', up_down)
    else:
        up_down = []
        for k in rep_dict:
            cluster_dict_k = rep_dict[k]
            for c in cluster_dict_k:
                prots_in_c = cluster_dict_k[c]
                nprots = len(prots_in_c)
                up_down.append((k, prots_in_c[0]-1,prots_in_c[-1]+1,nprots))
        return('False', up_down)

#Search HMMs in flanking proteins
def search_in_synteny(genbankDict, up_down, include_file, up, down, genome_folder, evalue, hmm_coverage):
    include_in_synteny_summary = []
    for cluster in up_down:
        cluster_rep = cluster[0]
        cluster_range_up = [ f'UID{i}' for i in range(cluster[1]-up,cluster[1]+1) ]
        cluster_range_down = [ f'UID{i}' for i in range(cluster[2],cluster[2]+down+1) ]
        uids_to_test = cluster_range_up + cluster_range_down
        #extract proteins
        proteomeFile = os.path.join(genome_folder, "protein_cluster.faa")
        with open(proteomeFile, "w+") as f:
            for prot in uids_to_test:
                key = cluster_rep
                if prot in genbankDict[key]:
                    val = genbankDict[key][prot]
                    f.write(">"+key+"|"+val[6]+"|"+val[0]+"\n"+str(val[2])+"\n")
        #hmm search
        hmmHits, protein_subset = hmmsearch(proteomeFile, include_file, genome_folder, evalue, hmm_coverage)
        if len(hmmHits) > 0:
            domains = get_domain_names([include_file])
            domains.sort()
            fdomains = list(hmmHits['query name'].unique())
            fdomains.sort()
            domains_string = ' ,'.join( [ f'{r["target name"]}|{r["query name"]}' for _,r in hmmHits.iterrows() ] )
            if fdomains == domains:
                include_in_synteny_summary.append(f'All syntenic genes found: {domains_string}')
            else:
                include_in_synteny_summary.append(f'Some syntenic genes found: {domains_string}')
        else:
            include_in_synteny_summary.append('')
    return(include_in_synteny_summary)

def cat(file_list, output_file):
    with open(output_file, 'w') as outfile:
        for filename in file_list:
            with open(filename, 'r') as infile:
                for line in infile:
                    outfile.write(line)

######################################################################################################
################################## Main ##############################################################
######################################################################################################

if __name__ == '__main__':
    #
    args = parseArgs()
    printSoftName()
    #
    #Get genomes
    #genome folder and files
    
    if args.inputFolder:
        inputFolder = args.inputFolder
        if not os.path.exists(inputFolder):
            sys.exit("Input folder not found.")
        #detect if there are Gb files in folder, formats .gb, gbk, ...
        if not args.GbExtension:
            GbExtension = "gb"
            genomes = glob.glob(os.path.join(inputFolder,"*"+GbExtension))
            if len(genomes) == 0:
                sys.exit("No Gb files in the input folder. Try providing the genbank file extension with -x argument.")
        else:
            GbExtension = args.GbExtension
            genomes = glob.glob(os.path.join(inputFolder,"*"+GbExtension))
        if len(genomes) == 0:
            sys.exit("No Gb files were found in the input folder.")
    elif args.queryGenus:
        #For NCBI download
        queryGenus = args.queryGenus
        if args.downloadFolder:
            inputFolder = args.downloadFolder
        else:
            inputFolder = "Results"
        if not os.path.exists(inputFolder):
            os.mkdir(inputFolder)
        else:
            from datetime import datetime
            import re
            if not args.overwrite:
                inputFolder2 = inputFolder+'_'+re.sub('[:\-_ .]','',str(datetime.now()))
                print(f'{inputFolder} folder exists. The results will be saved in {inputFolder2}...')
                inputFolder = inputFolder2
                os.mkdir(inputFolder)
        GbExtension = "gb"
        if args.numSeqs:
            numSeqs = int(args.numSeqs)
        else:
            numSeqs = -1
        print(f'Downloading genomes...')
        getSeqsFromNCBI(queryGenus,inputFolder,numSeqs)
        genomes = glob.glob(os.path.join(inputFolder,"*.gb"))
    elif (args.accessions) and (len(args.accessions) > 0):
        if args.downloadFolder:
            inputFolder = args.downloadFolder
        else:
            inputFolder = "Results"
        if not os.path.exists(inputFolder):
            os.mkdir(inputFolder)
        accessions = args.accessions
        for downloadid in accessions:
            print(f'Downloading genomes...')
            get_Genome_by_assembly_accession(downloadid,inputFolder)
            genomes = glob.glob(os.path.join(inputFolder,"*.gb"))
    else:
        sys.exit("Please input the path of the folder containing the genomes to analyze, or a Genus to download from NCBI.")
    #
    #Ref_genome file
    refGenome = args.refGenome
    refacc = args.refacc
    if refGenome:
        if os.path.exists(refGenome):
            file_name = os.path.basename(refGenome)
            new_file_name = os.path.join(inputFolder,file_name)
            if refGenome != new_file_name:
                os.rename(refGenome, new_file_name)
                refGenome = new_file_name
            print(f'\nUsing {refGenome} file as reference genome...')
            if refGenome in genomes:
                genomes.remove(refGenome)
            genomes.insert(0,refGenome)
        else:
            sys.exit(f'File {refGenome} not found...')
    elif refacc:
        if len(refacc) > 0:
            refGenome = get_Genome_by_assembly_accession(refacc,inputFolder)
        if refGenome in genomes:
            genomes.remove(refGenome)
        genomes.insert(0,refGenome)
        print(f'\nDownloading {refacc} reference genome from NCBI...')
    else:
        print(f'\nA reference genome was not specified...')
    #
    #hmms and files
    pfamids = args.pfamid #default []
    hmms = args.hmms #default []
    #
    pfam_hmm_folder = os.path.join(inputFolder,'HMM')
    if not os.path.exists(pfam_hmm_folder):
        os.mkdir(pfam_hmm_folder)
    #
    architectures = []
    excl = []
    incl = []
    #
    force_synteny = args.force_synteny
    min_size = int(args.min_size)
    max_size = int(args.max_size)
    architecture_size_lims = [(min_size,max_size)]
    #
    if args.architecture_pfamid:
        architecture_pfamid = args.architecture_pfamid
        print(f'\nArchitecture file {architecture_pfamid} will be used...')
        if not os.path.exists(architecture_pfamid):
            sys.exit(f'{architecture_pfamid} not found...')
        pfam_architectures = []
        pfam_arch_excludes = []
        pfam_arch_includes = []
        architecture_size_lims = []
        with open(architecture_pfamid,'r') as af:
            for line in af:
                if line.startswith('#'):
                    continue
                line = line.strip()
                line = re.sub('  ',' ',line)
                line = re.sub(', ',',',line)
                line = re.sub(' ,',',',line)
                pfam_lists = line.split(',')
                pfam_arch = pfam_lists[0].split(' ')
                pfam_architectures.append(pfam_arch)
                pfam_ex = pfam_lists[1].split(' ')
                pfam_arch_excludes.append(pfam_ex)
                pfam_in = pfam_lists[2].split(' ')
                pfam_arch_includes.append(pfam_in)
                mins = pfam_lists[3]
                maxs = pfam_lists[4]
                if (len(mins) >0) and (mins.isnumeric()):
                    min_size = int(pfam_lists[3])
                else:
                    min_size = 0
                if (len(maxs) >0) and (maxs.isnumeric()):
                    max_size = int(pfam_lists[4])
                else:
                    max_size = 50000
                architecture_size_lims.append((min_size,max_size))
        for pa in pfam_architectures:
            architectures.append(download_architecture_pfams(pa,'arch'))
        for pa in pfam_arch_excludes:
            excl.append(download_architecture_pfams(pa,'excl'))
        for pa in pfam_arch_includes:
            incl.append(download_architecture_pfams(pa,'inlc'))
    elif args.architecture_file:
        architecture_file = args.architecture_file
        print(f'\nArchitecture file {architecture_file} will be used...')
        if not os.path.exists(architecture_file):
            sys.exit(f'{architecture_file} not found...') 
        architecture_size_lims = []   
        with open(architecture_file,'r') as af:
            for line in af:
                if line.startswith('#'):
                    continue
                line = line.strip()
                line = re.sub('  ',' ',line)
                line = re.sub(', ',',',line)
                line = re.sub(' ,',',',line)
                pfam_lists = line.split(',')
                pfam_arch = pfam_lists[0].split(' ')
                architectures.append(pfam_arch)
                pfam_ex = pfam_lists[1].split(' ')
                excl.append(pfam_ex)
                pfam_in = pfam_lists[2].split(' ')
                incl.append(pfam_in)
                mins = pfam_lists[3]
                maxs = pfam_lists[4]
                if (len(mins) >0) and (mins.isnumeric()):
                    min_size = int(pfam_lists[3])
                else:
                    min_size = 0
                if (len(maxs) >0) and (maxs.isnumeric()):
                    max_size = int(pfam_lists[4])
                else:
                    max_size = 50000
                architecture_size_lims.append((min_size,max_size))
        for a in architectures:
            test_file_list(a)
        for ex in excl:
            test_file_list(ex)        
        for inc in incl:
            test_file_list(inc)
    elif len(pfamids) > 0:
        print('\nPFAM HMM will be used...')
        architectures = [download_architecture_pfams(pfamids,"arch")]
    elif len(hmms) > 0:
        print('\HMM files will be used...')
        test_file_list(hmms)
        architectures = [hmms]
    else:
        sys.exit('\nAn HMM file, PFAM ID or Architecture file must be provided.')
    #
    #
    print()
    if (not args.architecture_pfamid) and (not args.architecture_file):
        #exclude pfamids
        exclude_pfamids = args.exclude_pfamids
        exclude_hmms = []
        if len(args.exclude_hmms) > 0:
             ehmms = args.exclude_hmms
             for hmmf in ehmms:
                 if os.path.exists(hmmf):
                     exclude_hmms.append(hmmf)
             print(f'Proteins with a match to these HMMs files will be excluded: {", ".join(exclude_hmms)}')
        if len(exclude_pfamids) > 0:
            print(f'Proteins with a match to these PFAM HMMs will be excluded: {", ".join(exclude_pfamids)}')
            for pfamid in exclude_pfamids:
                pfamid = pfamid.upper()
                pfam_file = os.path.join(pfam_hmm_folder, pfamid+'.hmm')
                if not pfam_file in exclude_hmms:
                    download_pfam_hmm(pfamid, pfam_file)
                exclude_hmms.append(pfam_file)
        elif args.exclude_pfamids_table:
            exclude_pfamids_table = args.exclude_pfamids_table
            with open(exclude_pfamids_table,'r') as f:
                exclude_pfamids = [line.strip() for line in f]
            for pfamid in exclude_pfamids:
                pfamid = pfamid.upper()
                pfam_file = os.path.join(pfam_hmm_folder, pfamid+'.hmm')
                if not pfam_file in exclude_hmms:
                    download_pfam_hmm(pfamid, pfam_file)
                exclude_hmms.append(pfam_file)
        #
        test_file_list(exclude_hmms)
        excl = [exclude_hmms]
        #include_in_synteny
        include_in_synteny = args.include_in_synteny
        include_hmms = []
        if len(args.include_hmms) > 0:
            ihmms = args.include_hmms
            for hmmf in ihmms:
                if os.path.exists(hmmf):
                    include_hmms.append(hmmf)
            print(f'The following hmm profiles will be searched in a specified number of "up" and "down" proteins: {", ".join(include_hmms)}')
        if len(include_in_synteny) > 0:
            print(f'The following PFAM profiles will be searched in a specified number of "up" and "down" proteins: {", ".join(include_in_synteny)}')
            for pfamid in include_in_synteny:
                pfamid = pfamid.upper()
                pfam_file = os.path.join(pfam_hmm_folder, pfamid+'.hmm')
                if not pfam_file in include_hmms:
                    download_pfam_hmm(pfamid, pfam_file)
                include_hmms.append(pfam_file)
        #
        test_file_list(include_hmms)
        incl = [include_hmms]
    #
    #For each genome
    #search for profiles
    #search parameters and HMMER options
    #ncpu = int(args.ncpu)
    evalue = float(args.evalue)
    hmm_coverage = int(args.hmm_coverage)
    #
    # synteny options
    tolerance = args.tolerance
    more_than_one_contiguous = args.more_than_one_contiguous
    up = int(args.up)
    down = int(args.down)
    #
    summary_table = pd.DataFrame()
    for genome_file in genomes:
        print(f'\nProcessing {genome_file}...')
        if genome_file == refGenome:
            reference = True
        else:
            reference = False
        gsummary_table = pd.DataFrame()
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]
        genome_folder = os.path.join(inputFolder, genome_name)
        if not os.path.exists(genome_folder):
            os.mkdir(genome_folder)
        #
        genbankDict, description, proteomeFile = load_genome(genome_file)
        if not proteomeFile:
            continue
        #### 
        for arch_idx in range(len(architectures)):
            arch = architectures[arch_idx]
            arch_string = '-'.join([ os.path.splitext(os.path.basename(a))[0] for a in arch])
            exclude_hmms = excl[arch_idx]
            include_hmms = incl[arch_idx]
            if (force_synteny) and (len(''.join(include_hmms)) == 0):
                print(f'--> Architecture {arch_string} skipped, force synteny is on and there are no HMM profiles to look for its syntenic region...')
                continue
            print(f'-> Architecture: {arch_string} ...')
            min_size, max_size = architecture_size_lims[arch_idx]
            #
            if exclude_hmms != []:
                exclude_file = os.path.join(pfam_hmm_folder,f'exclude{arch_idx}.hmm')
                cat(exclude_hmms, exclude_file)
            else:
                exclude_file = ''
            #
            if include_hmms != []:
                include_file = os.path.join(pfam_hmm_folder,f'include{arch_idx}.hmm')
                cat(include_hmms, include_file)
            else:
                include_file = ''
            #
            proteins_table = pd.DataFrame()
            proteins, dom_table = search_by_architecture(genome_file, genome_name, arch, genome_folder, evalue, hmm_coverage, exclude_file, genbankDict, proteomeFile, min_size, max_size)
            if len(dom_table) > 0:
                protein_names = list(proteins['target name'])
                proteins_table['proteins'] = protein_names
                proteins_table['Replicon'] = [ p.split('|')[0]  for p in protein_names]
                proteins_table['UID'] = [ p.split('|')[1]  for p in protein_names]
                proteins_table['genome_name'] = genome_name
                proteins_table['description'] = ' ,'.join([ f'{k}|{v}' for k,v in description.items()])
                #Evaluate synteny and contiguity
                #looks for provided HMMs or PFAMs on the genomic context of n up and n down genes. Fist determine if all the hits are contiguous, then extracts adjacent genes.
                print(f'-> Looking for clusters of proteins with the same domain architecture...')
                contiguity, up_down = evaluate_contiguity(protein_names, tolerance)
                #synteny of the architecture: look for expected genes in n up and n down positions
                # HOW TO PROCESS MULTIPLE ARCHITECTURES 
                #returns a list with the genes found or '' if only some or none of the include genes was found. Order is the same of up_down.
                if len(include_file) > 0:
                    print(f'-> Looking for conserved domains in neighbor proteins...')
                    include_in_synteny = search_in_synteny(genbankDict, up_down, include_file, up, down, genome_folder, evalue, hmm_coverage)
                else:
                    include_in_synteny = ''
                proteins_table['contiguity'] = contiguity
                proteins_table['Cluster Number'] = ''
                proteins_table['Synteny'] = ''
                cnumber = 0
                for i in range(len(up_down)):
                    cluster = up_down[i]
                    cluster_rep = cluster[0]
                    cluster_range = [ f'UID{i}' for i in range(cluster[1]+1,cluster[2]) ]
                    index = proteins_table[(proteins_table['Replicon'] == cluster_rep) & (proteins_table['UID'].isin(cluster_range))]['Cluster Number'].index
                    proteins_table.loc[index, 'Cluster Number'] = cnumber
                    proteins_table.loc[index, 'Protein Copies'] = cluster[3]
                    if len(include_in_synteny) > 0:
                        proteins_table.loc[index, 'Synteny'] = include_in_synteny[i]
                    cnumber += 1
                #select the one with more copies?
                if more_than_one_contiguous:
                    proteins_table = proteins_table[proteins_table['Protein Copies'] > 1]
                if force_synteny:
                    proteins_table = proteins_table[proteins_table['Synteny'] != '']
                # actualize summary_table
                proteins_table['Architecture'] = arch_string
                proteins_table['size']  = proteins['tlen']
                gsummary_table = pd.concat([gsummary_table,proteins_table])
        if reference:
            gsummary_table['Reference'] = '*'
        else:
            gsummary_table['Reference'] = ' '
        #
        summary_table = pd.concat([summary_table,gsummary_table])
    save_file = os.path.join(inputFolder, 'Architectures.csv')
    if len(summary_table) == 0:
        print(f'\nNo homologs with the desired domain architecture were found...')
        summary_table = pd.DataFrame(columns=['proteins','Replicon','UID','genome_name','description','contiguity','Cluster Number','Synteny','Protein Copies','Architecture','size','Reference'])
    summary_table.columns = ['Proteins','Replicon','UID','Genome name','Description','Contiguity','Cluster number','Synteny','Protein Copies','Architecture','Size','Reference']
    summary_table.to_csv(save_file, index=False)
