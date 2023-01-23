"""
The main file, which performs qc and annotation of single-cell count matrices.
"""

import numpy as np
import synapseclient 
import synapseutils 
import argparse

from ROSMAPwgs.download import get_chromosome_annotations, extract_filenames, filter_filenames, save_files

def main(args):
    '''Given arguments from `args` performs qc and annotation.
    '''
    print('loading dictionary')
    dictionary = np.load('./human_Release_19_GRCh37p13/chr_dictionary.npy', allow_pickle=True).item()

    print('synapse logging in')
    syn = synapseclient.Synapse() 
    syn.login(args.username, args.pw) 

    print('looking up chromosomes')
    chromosomes = get_chromosome_annotations(eval(args.gene_list), dictionary)
    
    print('filtering filenames')
    walkedPath = synapseutils.walk(syn, args.synapseID)
    filenames = extract_filenames(walkedPath)
    filenames_to_download = filter_filenames(filenames, args.extension, chromosomes)
    
    print('downloading files')
    save_files(syn, args.outdir, filenames_to_download)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='perform qc and annotation of single-cell count matrix.')
    parser.add_argument('--gene_list', type=str, required=True, default=False, help='genes for which to download variant info')
    parser.add_argument('--outdir', type=str, required=True, default='./rawdata', help='where to save the data')
    parser.add_argument('--username', type=str, required=True, default=False, help='synapse username; request an account at https://www.synapse.org/#!RegisterAccount:0')
    parser.add_argument('--pw', type=str, required=True, default=True, help='synapse password; request an account at https://www.synapse.org/#!RegisterAccount:0')
    parser.add_argument('--extension', type=str, required=False, default='annotated.vcf', help='filetype to download: see README for possible extensions')
    parser.add_argument('--synapseID', type=str, required=False, default='syn11724057', help='WGS Synapse ID')
  
    args = parser.parse_args()
    ext = set(('annotated.clinical.txt', 'annotated.coding_rare.txt', 'annotated.coding.txt', 'annotated.txt', 'annotated.vcf.gz', 'annotated.vcf.gz.tbi', 'recalibrated_variants.vcf.gz', 'recalibrated_variants.vcf.gz.tbi'))
    if args.extension not in ext:
        parser.error('--extension must be one of: {annotated.clinical.txt, annotated.coding_rare.txt, annotated.coding.txt, annotated.txt, annotated.vcf.gz, annotated.vcf.gz.tbi, recalibrated_variants.vcf.gz, recalibrated_variants.vcf.gz.tbi}')
    
    main(args)
    
    