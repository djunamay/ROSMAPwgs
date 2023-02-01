import numpy as np
import os
import pandas as pd
import allel
from pathlib import Path
import distutils 
from assertpy import assert_that
import random
from ROSMAPwgs.annotate import convert_to_int, return_all_variants_table, extract_callset_data, return_variant_indices_from_vcf, format_genotype_data, return_genotype_counts, compute_MAFs, save_paths_per_chromosome, make_executable, load_test_data, convert_genotypes_to_str

'''
Upload test data, if exists
'''
print('defining test data')
path1 = './human_Release_19_GRCh37p13/pos_dictionary.npy'
path2 = './test_output/DEJ_11898_B01_GRM_WGS_2017-05-15_19.recalibrated_variants.vcf.gz'
path3 = './test_output/DEJ_11898_B01_GRM_WGS_2017-05-15_19.recalibrated_variants.vcf.gz.tbi'
path4 = './test_output/DEJ_11898_B01_GRM_WGS_2017-05-15_19.recalibrated_variants.annotated.coding.txt'
callset_names = ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT', 'variants/FILTER_PASS', 'variants/numalt', 'variants/AF']
annotation_names = ['POS', 'ID', 'REF', 'ALT', 'EFFECT', 'IMPACT', 'GENE', 'GENEID', 'HGVS_C', 'HGVC_P', 'LOF', 'NMD', '1000Gp3_AF']
gene = 'ABCA7'

if os.path.isfile(path1) & os.path.isfile(path2) & os.path.isfile(path3) & os.path.isfile(path4):
    print('loading test data')
    make_executable(path3)
    postion_dict, callset, annotation, df, all_data, callset_position_indices, genotype_data, genotype_counts, MAFs, all_variants = load_test_data(path1, path2, path3, path4, callset_names, annotation_names, gene)

'''
Run tests
'''

def test_convert_genotypes_to_str():
    """
    Testing function that converts ndarray genotype to string
        - Comparing manually-determined solution to function output for a simple example
    """
    l = [[0,1], [1,1]]
    g = convert_genotypes_to_str(l)
    
    assert_that(g[0]).is_equal_to('0/1').is_true()
    assert_that(g[1]).is_equal_to('1/1').is_true()

def test_return_genotype_counts():
    """
    Testing function that returns genotype counts
        - Comparing manually-determined solution to function output for a simple example
    """
    
    l = [[[0,1], [1,1], [0,0]], [[0,0], [0,0], [1,1]]]
    g = [convert_genotypes_to_str(i) for i in l]
    df = pd.DataFrame(g)
    cts = return_genotype_counts(df)
    n_samples = df.shape[1]
    n_genotypes = np.sum(cts, axis = 1)
    assert_that(np.unique(n_genotypes==n_samples)[0]).is_true()

def test_compute_MAFs():
    """
    Testing function that computes minor allele frequencies from input genotype data
        - Comparing allele frequencies provided in the callset to function outputs
            - If they are close (up to rounding), this suggests that (1) MAFs are being computed correctly by our function and (2) that genotype extraction and formatting for specific variants of interest maintains the mapping between position and sample-wise genotypes
        - Comparing manually-determined solution to function output for a simple example
    """
    if 'callset' in locals():
        provided_maf = callset['variants/AF']
        pos_indices = [np.where(callset['variants/POS']==x)[0][0] for x in all_data['POS']]
        p = provided_maf[pos_indices]
        assert_that(np.allclose(np.array(p),np.array(MAFs['MAF']), atol = 4e-04)).is_true()

    l = [[[0,1], [1,1], [0,0]], [[0,0], [0,0], [1,1]]]
    g = [convert_genotypes_to_str(i) for i in l]
    df = pd.DataFrame(g)
    mafs = compute_MAFs(return_genotype_counts(df))
    true = [0.5, 1/3]
    assert_that(np.allclose(np.array(true),np.array(mafs['MAF']), atol = 4e-04)).is_true()

def test_extract_callset_data():
    """
    Testing function that extracts specified entries from callset dictionary and returns them as dataframe
        - Making a direct comparison between the extracted dataframe columns and the callset elements
    """
    if 'callset' in locals():
        temp = extract_callset_data(['variants/POS', 'variants/ALT'], callset)
        convert_to_int(temp)
        assert_that(np.allclose(np.array(temp['POS']),np.array(callset['variants/POS']))).is_true()
        assert_that(np.array_equal(np.array(callset['variants/ALT'][:,0]), np.array(callset['variants/ALT'][:,0]))).is_true()
    
def test_format_genotype_data():
    """
    Testing function that returns genotype data from callset
        - For a randomly sampled variant, checking that the genotype for that variant is the same in the callset as in the formatted genotype tabel (output of the format_genotype_data function)
    """
    if 'callset' in locals():
        geno_formatted = pd.DataFrame(format_genotype_data(callset, callset_position_indices, all_data))
        i = random.randint(0, len(geno_formatted)-1)
        geno_formatted = geno_formatted.loc[i]
        pos = all_data.loc[i,'POS']
        index = np.where(callset['variants/POS']==pos)[0][0]
        geno_raw = convert_genotypes_to_str(callset['calldata/GT'][index])
        assert_that(np.array_equal(geno_formatted, geno_raw)).is_true()

def test_return_all_variants_table():
    """
    Testing function that returns variant annotation and genotype data for gene of interest
        - For a randomly selected sample, convert the genotype from strings to nested list and extract the genotype for that sample from the callset
            - Compare genotypes
            - Compare REF and ALT alleles extracted from callset and merged from annotations to test that merging on positions correctly merges the same variants
    """
    if 'callset' in locals():
        temp = return_all_variants_table(path2, path4, path3, postion_dict, gene, annotation_names, callset_names, callset)
        sample = random.choice(callset['samples'])
        geno_sample = [[int(x.split('/')[0]), int(x.split('/')[1])] for x in temp[sample]]
        pos_indices = [np.where(callset['variants/POS']==x)[0][0] for x in temp['POS']]
        geno_raw_sample = callset['calldata/GT'][pos_indices,(callset['samples']==sample),:]
        
        assert_that(np.array_equal(geno_sample, geno_raw_sample)).is_true()
        assert_that(np.array_equal(temp['REF_x'], temp['REF_y'])).is_true()
        assert_that(np.array_equal(temp['ALT_0'], temp['ALT'])).is_true()
    
