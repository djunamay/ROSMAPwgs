import numpy as np
import synapseclient 
import synapseutils 
import uuid

def get_chromosome_annotations(gene_list, dictionary):
    '''
    returns normalized matrix of marker genes
    Args:
        filtered_counts ndarray
            2D array of counts post-cell-filtering
            N-cells x N-features
        marker_indices list
            column indices indicating marker genes in filtered_counts matrix
        marker_out numpy memmap
            N-cells x len(marker_indices)
        total_counts ndarray
            1D array of per-cell total counts, length = N-cells
    '''
    chromosomes = [dictionary[gene].replace("chr", "") for gene in gene_list]
    return set(chromosomes)

def extract_filenames(walkedPath):
    '''
    returns normalized matrix of marker genes
    Args:
        filtered_counts ndarray
            2D array of counts post-cell-filtering
            N-cells x N-features
        marker_indices list
            column indices indicating marker genes in filtered_counts matrix
        marker_out numpy memmap
            N-cells x len(marker_indices)
        total_counts ndarray
            1D array of per-cell total counts, length = N-cells
    '''
    filenames = []
    for dirpath, dirname, filename in walkedPath:
        filenames.append(filename)
    return filenames

def filter_filenames(filenames, extension, chromosomes):
    '''
    returns normalized matrix of marker genes
    Args:
        filtered_counts ndarray
            2D array of counts post-cell-filtering
            N-cells x N-features
        marker_indices list
            column indices indicating marker genes in filtered_counts matrix
        marker_out numpy memmap
            N-cells x len(marker_indices)
        total_counts ndarray
            1D array of per-cell total counts, length = N-cells
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
    returns normalized matrix of marker genes
    Args:
        filtered_counts ndarray
            2D array of counts post-cell-filtering
            N-cells x N-features
        marker_indices list
            column indices indicating marker genes in filtered_counts matrix
        marker_out numpy memmap
            N-cells x len(marker_indices)
        total_counts ndarray
            1D array of per-cell total counts, length = N-cells
    '''
    with open(outdir + '/manifest_' + str(uuid.uuid1()) + '.txt', 'w') as f:
        for i in filenames_to_download:
            files = synapseutils.syncFromSynapse(syn, i, path = outdir) 
            f.write(str(files[0]))
    print('done.')