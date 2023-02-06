import numpy as np
import pandas as pd
import uuid
import allel
import os as os

def format_annotations(callset, name):
    '''
    returns set of chromosome annotations for gene input list
    Args:
        gene_list list
            input genes
        dictionary dictionary
            gene-to-chromosome mapping   
    '''
    short_name = name.split('/')[1]
    temp_anno = callset[name]
    shp = temp_anno.shape
    if len(shp)==1:
        temp_anno = temp_anno.reshape(-1,1)
        columns = short_name
    else:
        columns = [short_name + '_' + str(x) for x in range(shp[1])]
    return columns, temp_anno

def flatten_list(L):
    '''
    converts 2D to 1D list
    Args:
        L 2D list to convert
    '''
    out = []
    for i in L:
        if isinstance(i, list):
            for x in i:
                out.append(x)
        else:
            out.append(i)
    return out

def return_variant_indices_from_vcf(pos):
    '''
    return dictionary mapping genomic position of each variant in the vcf file to its index in the vcf file
    more generally, maps input vector to sequence of integers from [0, N), where N is the length of the input vector
    Args:
        pos ndarray 
            int32 position indices for each variant in the vcf file, in order of appearance in the vcf file
    '''
    indices = dict(zip(pos, range(len(pos))))
    return indices

def convert_genotypes_to_str(genotype):
    '''
    converts 2D list of per-subject genotypes for a single genotype to 1D list
    e.g.: [[0,1], [1,1]] becomes ['0/1', '1/1']
    Args:
        genotype list
            2D list of genotypes for a particular locus, where each sublists consists of two elements, each encoding a single allele at that locus for a given individual
    '''
    converted = [str(i[0]) + '/' + str(i[1]) for i in genotype]
    return converted

def return_genotype_counts(genotype_data):
    '''
    returns pandas dataframe of absolute frequencies of each allele (columns) across variants of interest (rows)
    Args: 
        genotype_data pandas dataframe
            output from format_genotype_data()
    '''
    g = np.unique(genotype_data)
    counts = pd.DataFrame([[np.sum(genotype_data.iloc[i]==x) for x in g] for i in range(genotype_data.shape[0])])
    counts.columns = 'N ' + g
    return counts

def compute_MAFs(genotype_counts):
    '''
    returns pandas dataframe of minor allele frequencies across variants of interest (rows) estimated from the genotype input
    Args: 
        genotype_data pandas dataframe
            output from format_genotype_data()
    '''
    denominator = np.sum(genotype_counts*2, axis = 1)
    N_Alt_alleles_per_col = [np.sum(np.array(x.split(' ')[-1].split('/'))=='1') for x in genotype_counts.columns]
    MAF = pd.DataFrame(np.sum(genotype_counts * N_Alt_alleles_per_col, axis = 1)/denominator)
    MAF.columns = ['MAF']
    return MAF

def extract_callset_data(names, callset):
    '''
    returns pandas dataframe from dictionary for specified keys
    Args:
        callset dictionary
            vcf file loaded with allel.read_vcf()
        names list
            callset keys to extract from dictionary
    '''
    colnames = flatten_list([format_annotations(callset, i)[0] for i in names])
    vals = [format_annotations(callset, i)[1] for i in names]

    data_output = np.empty(shape=(len(vals[0]), len(colnames)), dtype=object)
    start = 0
    for i in vals:
        end = i.shape[1]+start
        data_output[:,start:end] = i
        start=end

    df = pd.DataFrame(data_output)
    df.columns = colnames
    return df 

def format_genotype_data(callset, callset_position_indices, all_data):
    '''
    returns pandas dataframe of N_variants x N_subjects, where each value represents the genotype across alleles for all variants of interest in the all_data input dataframe
    Args:
        callset dictionary
            vcf file loaded with allel.read_vcf()
            per-variant per-individual genotype information stored with key label 'calldata/GT' (array of shape N_variants x N_subjects x N_alleles)
            sample labels stored with key label 'samples' (array of shape 1 x N_subjects)
        callset_position_indices dictionary
            output from return_variant_indices_from_vcf()
        all_data pandas dataframe
            variant annotation for all variants of interest with variant position information in 'POS' column as int32
    '''
    genotype_data = callset['calldata/GT'][[callset_position_indices[x] for x in all_data['POS']]]
    genotype_data = pd.DataFrame([convert_genotypes_to_str(x) for x in genotype_data])
    genotype_data.columns = callset['samples']
    return genotype_data

def return_all_variants_table(path_to_vcf, path_to_annotations, path_to_tbi, postion_dict, gene, annotation_names, callset_names, callset=None):
    '''
    Args:
        path_to_vcf str
            path to vcf file
        path_to_annotations str
            path to variant annotation text file with variant position column labeled 'POS'
        path_to_tbi str
            path to vcf index file
        postion_dict
            dictionary mapping genes to chromosome positions (e.g. '19:45409011-45412650')
        gene str
            gene for which to return variant information
        annotation_names list
            list of column names to extract from annotation text file, must include 'POS' column with genomic variant position
        callset_names list
            list of key names to extract from vcf callset, must include 'variants/POS' and 'variants/'FILTER_PASS' with boolean indicating whether to keep (True) variants based on QC criteria
    '''
    annotation = pd.read_table(path_to_annotations, dtype=object)[annotation_names]
    
    if gene not in set(annotation['GENE']):
        raise Exception("Gene is not in annotation file. Try a different gene or select an annotation file that contains annotations for this gene.")
     
    if callset is None:
        print('loading vcf')
        callset = allel.read_vcf(path_to_vcf, fields = '*', region=postion_dict[gene], tabix=path_to_tbi) 
        
    print('vcf loaded.')
   
    print('extracting variant info')
    df = extract_callset_data(callset_names, callset)
    convert_to_int(df)
    convert_to_int(annotation)
    all_data = pd.merge(df, annotation, on = 'POS')
    all_data = all_data.loc[all_data['FILTER_PASS']]
    all_data.index = range(all_data.shape[0])
    
    callset_position_indices = return_variant_indices_from_vcf(pos = callset['variants/POS'])
    
    genotype_data = format_genotype_data(callset, callset_position_indices, all_data)
    
    genotype_counts = return_genotype_counts(genotype_data)
    MAFs = compute_MAFs(genotype_counts)
    
    all_variants = pd.concat((all_data, genotype_counts, MAFs, genotype_data), axis = 1)
    print('done')
    return all_variants

def save_paths_per_chromosome(sele_files, path):
    '''
    returns path to vcf, index, and annotation text file for input array of filenames
    Args:
        sele_paths ndarray
            file names must include a single vcf file ending in 'qz', a single vcf index file ending in 'tbi' and a single variant annotation text file ending in 'txt' for the chromosome of interest 
    '''
    path_to_vcf = path + '/' + sele_files[[s.endswith('gz') for s in sele_files]][0]
    path_to_tbi = path + '/' + sele_files[[s.endswith('tbi') for s in sele_files]][0]
    path_to_annotations = path + '/' + sele_files[[s.endswith('txt') for s in sele_files]][0]
    return path_to_vcf, path_to_tbi, path_to_annotations

def make_executable(path): # source: https://stackoverflow.com/questions/12791997/how-do-you-do-a-simple-chmod-x-from-within-python
    '''
    makes a given file chmod + x executable
    Args:
        path str
            path to file for which permissions should be changed
    '''
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)
    
def convert_to_int(df):
    '''
    Converts 'POS' column in pandas dataframe to int32
    Args:
        df pandas dataframe
            with genomic position information under column name 'POS'
    '''
    df['POS'] = df['POS'].astype('int32')
    
    
def load_test_data(path1, path2, path3, path4, callset_names, annotation_names, gene):
    '''
    Test function
    '''
    # loading the data
    annotation = pd.read_table(path4, dtype=object)[annotation_names]
    if gene not in set(annotation['GENE']):
        raise Exception("Gene is not in annotation file.")
    postion_dict = np.load(path1, allow_pickle=True).item()    
    callset = allel.read_vcf(path2, fields = '*', region=postion_dict[gene], tabix=path3) 
    # extracting callset data and merging with annotation
    df = extract_callset_data(callset_names, callset)
    convert_to_int(df)
    convert_to_int(annotation)
    all_data = pd.merge(df, annotation, on = 'POS')
    all_data = all_data.loc[all_data['FILTER_PASS']]
    all_data.index = range(all_data.shape[0])
    # extracting genotype data for variants of interest
    callset_position_indices = return_variant_indices_from_vcf(pos = callset['variants/POS'])
    genotype_data = format_genotype_data(callset, callset_position_indices, all_data)
    # computing genotype metrics
    genotype_counts = return_genotype_counts(genotype_data)
    MAFs = compute_MAFs(genotype_counts)
    # combining all datasets
    all_variants = pd.concat((all_data, genotype_counts, MAFs, genotype_data), axis = 1)
    return postion_dict, callset, annotation, df, all_data, callset_position_indices, genotype_data, genotype_counts, MAFs, all_variants
