import numpy as np
from assertpy import assert_that
from ROSMAPwgs.annotate import return_all_variants_table, extract_callset_data, return_variant_indices_from_vcf, format_genotype_data, return_genotype_counts, compute_MAFs, save_paths_per_chromosome, make_executable


def convert_genotypes_to_str():
    """
    Testing function that converts ndarray genotype to string
    """

    assert_that(marker_out.shape[1]).is_equal_to(np.count_nonzero(marker_indices)).is_true()
    assert_that(marker_out.shape[0]).is_equal_to(np.count_nonzero(keep_cells)).is_true()

def return_genotype_counts():
    """
    Testing function that returns genotype counts
    """

    assert_that(marker_out.shape[1]).is_equal_to(np.count_nonzero(marker_indices)).is_true()
    assert_that(marker_out.shape[0]).is_equal_to(np.count_nonzero(keep_cells)).is_true()

def compute_MAFs():
    """
    Testing function that computes minor allele frequencies from input genotype data
    """
    
    ## compare the minor allele frequencies to the one provided by synapse in the callset

    assert_that(marker_out.shape[1]).is_equal_to(np.count_nonzero(marker_indices)).is_true()
    assert_that(marker_out.shape[0]).is_equal_to(np.count_nonzero(keep_cells)).is_true()

def extract_callset_data():
    """
    Testing function that extracts specified entries from callset dictionary and returns them as dataframe
    """

    assert_that(marker_out.shape[1]).is_equal_to(np.count_nonzero(marker_indices)).is_true()
    assert_that(marker_out.shape[0]).is_equal_to(np.count_nonzero(keep_cells)).is_true()

def format_genotype_data():
    """
    Testing function that returns genotype data from callset
    """

    assert_that(marker_out.shape[1]).is_equal_to(np.count_nonzero(marker_indices)).is_true()
    assert_that(marker_out.shape[0]).is_equal_to(np.count_nonzero(keep_cells)).is_true()

def return_all_variants_table():
    """
    Testing function that returns variant annotation and genotype data for gene of interest
    """

    assert_that(marker_out.shape[1]).is_equal_to(np.count_nonzero(marker_indices)).is_true()
    assert_that(marker_out.shape[0]).is_equal_to(np.count_nonzero(keep_cells)).is_true()
    

# is there a way and a point to test the main function?
