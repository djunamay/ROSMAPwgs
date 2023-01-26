import numpy as np
import pandas as pd
import uuid
import allel

def format_annotations(callset, name):
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
    out = []
    for i in L:
        if isinstance(i, list):
            for x in i:
                out.append(x)
        else:
            out.append(i)
    return out

def return_variant_indices_from_vcf(pos):
    indices = dict(zip(pos, range(len(pos))))
    return indices

def convert_genotypes_to_str(genotype):
    converted = [str(i[0]) + '/' + str(i[1]) for i in genotype]
    return converted

def return_genotype_counts(genotype_data):
    g = np.unique(genotype_data)
    counts = pd.DataFrame([[np.sum(genotype_data.iloc[i]==x) for x in g] for i in range(genotype_data.shape[0])])
    counts.columns = 'N ' + g
    return counts

def compute_MAFs(genotype_counts):
    MAF = pd.DataFrame((genotype_counts['N 0/1']*1 + genotype_counts['N 1/1']*2)/(genotype_counts['N 0/0']*2 + genotype_counts['N 0/1']*2 + genotype_counts['N 1/1']*2))
    MAF.columns = ['MAF']
    return MAF

def extract_callset_data(names, callset):
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
    genotype_data = callset['calldata/GT'][[callset_position_indices[x] for x in all_data['POS']]]
    genotype_data = pd.DataFrame([convert_genotypes_to_str(x) for x in genotype_data])
    genotype_data.columns = callset['samples']
    return genotype_data

def return_all_variants_table(path_to_vcf, path_to_annotations, path_to_tbi, postion_dict, gene, annotation_names, callset_names):
    print('loading vcf')
    callset = allel.read_vcf(path_to_vcf, fields = '*', region=postion_dict[gene], tabix=path_to_tbi) 
    annotation = pd.read_table(path_to_annotations, dtype=object)[annotation_names]
    print('extracting variant info')
    df = extract_callset_data(callset_names, callset)
    convert_to_int(df)
    convert_to_int(annotation)
    all_data = pd.merge(df, annotation, on = 'POS')
    all_data = all_data.loc[all_data['FILTER_PASS']]
    callset_position_indices = return_variant_indices_from_vcf(pos = callset['variants/POS'])
    genotype_data = format_genotype_data(callset, callset_position_indices, all_data)
    genotype_counts = return_genotype_counts(genotype_data)
    MAFs = compute_MAFs(genotype_counts)
    all_variants = pd.concat((all_data, genotype_counts, MAFs, genotype_data), axis = 1)
    print('done')
    return all_variants

def save_paths_per_chromosome(sele_files, path):
    path_to_vcf = path + '/' + sele_files[[s.endswith('gz') for s in sele_files]][0]
    path_to_tbi = path + '/' + sele_files[[s.endswith('tbi') for s in sele_files]][0]
    path_to_annotations = path + '/' + sele_files[[s.endswith('txt') for s in sele_files]][0]
    return path_to_vcf, path_to_tbi, path_to_annotations

def make_executable(path): # source: https://stackoverflow.com/questions/12791997/how-do-you-do-a-simple-chmod-x-from-within-python
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)
    
def convert_to_int(df):
    df['POS'] = df['POS'].astype('int32')