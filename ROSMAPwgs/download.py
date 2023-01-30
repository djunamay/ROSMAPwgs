import numpy as np
import synapseclient 
import synapseutils 
import uuid

def get_chromosome_annotations(gene_list, dictionary):
    '''
    returns set of chromosome annotations for gene input list
    Args:
        gene_list list
            input genes
        dictionary dictionary
            gene-to-chromosome mapping   
    '''
    chromosomes = [dictionary[gene].replace("chr", "") for gene in gene_list]
    return set(chromosomes)

def extract_filenames(walkedPath):
    '''
    returns list of filenames from synapse generator object
    Args:
        generator object 
            output from synapseutils.walk()
    '''
    filenames = []
    for dirpath, dirname, filename in walkedPath:
        filenames.append(filename)
    return filenames

def filter_filenames(filenames, extension, chromosomes):
    '''
    returns ndarray for input filenames filtered by extension and chromosome
    Args:
        filenames list
            output from extract_filenames()
        extension str
            filenames containing this extension will be kept
        chromosomes set
            output from get_chromosome_annotations()
    '''
    annotated = []    
    for f in filenames[0]:
        curr = f[0]
        if extension in curr:
            annotated.append(f)
    index = [x[0].split('.')[0].split('_')[-1] in chromosomes for x in annotated]
    sele = np.array(annotated)[index][:,-1]
    return sele

def save_files(syn, outdir, filenames_to_download):
    '''
    downloads files from synapse to target directory
    Args:
        syn str
            synapse ID pointing to WGS data
        outdir str
            where to save the downloaded data to
        filenames_to_download ndarray
            output from filter_filenames()
    '''
    with open(outdir + '/manifest_' + str(uuid.uuid1()) + '.txt', 'w') as f:
        for i in filenames_to_download:
            files = synapseutils.syncFromSynapse(syn, i, path = outdir) 
            f.write(str(files[0]))
    print('done.')