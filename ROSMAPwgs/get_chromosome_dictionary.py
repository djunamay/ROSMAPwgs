import pandas as pd
import numpy as np
gencode = pd.read_table('./human_Release_19_GRCh37p13/gencode.v19.annotation.gff3.gz', comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])
gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)
gencode_genes['names'] = gencode_genes['attribute'].str.split(';',expand=True).loc[:,5].str.split('=', expand=True).loc[:,1]
gencode_genes['chr'] = [x.replace('chr','') for x in gencode_genes['seqname']]

chromosome_dictionary = dict(zip(gencode_genes['names'], gencode_genes['seqname']))

sele = np.array(gencode_genes[['chr', 'start','end']])
gencode_genes['location'] = [str(x[0]) + ':' + str(x[1]) + '-' + str(x[2]) for x in sele]
position_dictionary = dict(zip(gencode_genes['names'], gencode_genes['location']))

np.save('../human_Release_19_GRCh37p13/chr_dictionary.npy', chromosome_dictionary) 
np.save('../human_Release_19_GRCh37p13/pos_dictionary.npy', position_dictionary) 