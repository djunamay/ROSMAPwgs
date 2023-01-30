"""
The main file, which performs qc and annotation of single-cell count matrices.
"""

import numpy as np
import synapseclient 
import synapseutils 
import argparse
import os
import pandas as pd
import allel
from pathlib import Path
import stat
import os as os
import distutils 

from ROSMAPwgs.download import get_chromosome_annotations, extract_filenames, filter_filenames, save_files
from ROSMAPwgs.annotate import return_all_variants_table, extract_callset_data, return_variant_indices_from_vcf, format_genotype_data, return_genotype_counts, compute_MAFs, save_paths_per_chromosome, make_executable

def main(args):
    '''Given arguments from `args` performs qc and annotation.
    '''
    print('loading dictionary')
    dictionary = np.load('./human_Release_19_GRCh37p13/chr_dictionary.npy', allow_pickle=True).item()
    if args.download is True:

        print('synapse logging in')
        syn = synapseclient.Synapse() 
        syn.login(args.username, args.pw) 

        print('looking up chromosomes')
        chromosomes = get_chromosome_annotations(eval(args.gene_list), dictionary)

        print('filtering filenames')
        walkedPath = synapseutils.walk(syn, args.synapseID)
        filenames = extract_filenames(walkedPath)
        filenames_to_download = filter_filenames(filenames, args.extension, chromosomes)

        if args.extract_HIGHandMED_annotations is True:
            filenames_to_download = np.concatenate((filenames_to_download, filter_filenames(filenames, 'annotated.coding.txt', chromosomes)))

        print('downloading files')
        save_files(syn, args.outdir, filenames_to_download)
    
    if args.extract_HIGHandMED_annotations is True:
        print('extracting medium and high annotation info')
        postion_dict = np.load('./human_Release_19_GRCh37p13/pos_dictionary.npy', allow_pickle=True).item()
        callset_names = eval(args.callset_names)
        annotation_names = eval(args.annotation_names)
        gene_list = eval(args.gene_list)
        
        files = os.listdir(args.outdir)
        files = np.array(files)[[file.startswith('DE') for file in files]]
        chromes = np.array([y.split('_')[-1] for y in [x.split('.')[0] for x in files]])
        chromes_renamed = ['chr' + x for x in chromes]
        
        gene_to_chromosome_mapping = pd.DataFrame.from_dict(dictionary, orient = 'index').loc[gene_list]
        targets = set(gene_to_chromosome_mapping[0])
        chromes_sele = np.array(chromes)[[x in targets for x in chromes_renamed]]
        
        outputs = []
        for c in np.unique(chromes_sele):
            print('--for chromosome ' + str(c))
            path_to_vcf, path_to_tbi, path_to_annotations = save_paths_per_chromosome(files[chromes==c], args.outdir)
            make_executable(path_to_tbi)
            current_genes = np.array(gene_list)[np.array(gene_to_chromosome_mapping == 'chr' + c).reshape(-1)]
            
            for gene in current_genes:
                print('---gene ' + str(gene))
                outputs.append(return_all_variants_table(path_to_vcf, path_to_annotations, path_to_tbi, postion_dict, gene, annotation_names, callset_names))
            
        if len(outputs)==1:
            result = outputs[0]
        else:
            result = pd.concat(outputs, ignore_index=True, sort=False)
        result.to_csv(args.outdir + '/HIGHandMED_coding_annotations_' + args.synapseID + '_subset.csv')
        print('done.')
        
# get main script working
# write tests for each function
# understand all the QC metrics
        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='perform qc and annotation of single-cell count matrix.')
    parser.add_argument('--gene_list', type=str, required=False, default=False, help="dtype str (list as str); genes for which to extract variant info, e.g: " + "['APOE', 'ABCA1']")
    parser.add_argument('--outdir', type=str, required=False, default=False, help='dtype str; output directory for downloaded data')
    parser.add_argument('--username', type=str, required=False, default=False, help='dtype str; synapse username; request an account at https://www.synapse.org/#!RegisterAccount:0')
    parser.add_argument('--pw', type=str, required=False, default=True, help='dtype str; synapse password; request an account at https://www.synapse.org/#!RegisterAccount:0')
    parser.add_argument('--extension', type=str, required=False, default='recalibrated_variants.vcf.gz', help='dtype str; filetype to download: --extension must be one of: {annotated.clinical.txt, annotated.coding_rare.txt, annotated.coding.txt, annotated.txt, annotated.vcf.gz, annotated.vcf.gz.tbi, recalibrated_variants.vcf.gz, recalibrated_variants.vcf.gz.tbi}')
    parser.add_argument('--synapseID', type=str, required=False, default='syn11724057', help='dtype str; WGS Synapse ID')
    parser.add_argument('--extract_HIGHandMED_annotations', type=lambda x:bool(distutils.util.strtobool(x)), required=False, default=False, help='dtype bool; whether to extract coding variants and their annotations as text file; if True, downloaded --extension must be "recalibrated_variants.vcf.gz"')
    parser.add_argument('--callset_names', type=str, required=False, default="['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT', 'variants/FILTER_PASS', 'variants/numalt', 'variants/AF']", help='dtype str (list as str); which columns of the vcf file to extract')
    parser.add_argument('--annotation_names', type=str, required=False, default="['POS', 'ID', 'REF', 'ALT', 'EFFECT', 'IMPACT', 'GENE', 'GENEID', 'HGVS_C', 'HGVC_P', 'LOF', 'NMD', '1000Gp3_AF']", help="dtype str (list as str); which columns of the annotation text file to extract; e.g. ['POS', 'ID']")
    parser.add_argument('--download', type=lambda x:bool(distutils.util.strtobool(x)), required=False, default=False, help='dtype bool; if files have already been downloaded, --outdir must indicate location of these files')

        
    args = parser.parse_args()
    ext = set(('annotated.clinical.txt', 'annotated.coding_rare.txt', 'annotated.coding.txt', 'annotated.txt', 'annotated.vcf.gz', 'annotated.vcf.gz.tbi', 'recalibrated_variants.vcf.gz', 'recalibrated_variants.vcf.gz.tbi'))
    if args.extension not in ext:
        parser.error('--extension must be one of: {annotated.clinical.txt, annotated.coding_rare.txt, annotated.coding.txt, annotated.txt, annotated.vcf.gz, annotated.vcf.gz.tbi, recalibrated_variants.vcf.gz, recalibrated_variants.vcf.gz.tbi}')
        
    main(args)
    
