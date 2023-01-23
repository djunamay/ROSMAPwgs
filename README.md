# WGS_processing

### ROSMAP WGS data source: https://www.synapse.org/#!Synapse:syn11724057

File type	Description
annotated.clinical.txt	All low frequency HIGH/MODERATE annotated variants with possible clinical impact (from ClinVar ) in text file format
annotated.coding_rare.txt	All HIGH/MODERATE annotated variants with less than 5% allele frequency in 1000genomes and ExAC in text file format
annotated.coding.txt	All annotated variants with HIGH/MODERATE impact in text file format
annotated.txt	All variants with annotations in text file format
annotated.vcf.gz	All variants with annotations in variant call format (VCF)
annotated.vcf.gz.tbi	Index file for annotated VCF file
recalibrated_variants.vcf.gz	All variants in variant call format
recalibrated_variants.vcf.gz.tbi	Index file for VCF file with all variants


### First, Create conda env
```bash
conda create --name wgs_env
conda install -c bioconda tabix
conda install -c anaconda ipykernel
python -m ipykernel install --user --name=wgs_env
pip install scikit-allel
```

### Then, install the following to map genes to chromosomes
```bash
mkdir human_Release_19_GRCh37p13
cd human_Release_19_GRCh37p13
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz
```

### To download your data of interest, run the following
```bash
python main.py --gene_list "['ABCA7']" --outdir './test_output' --username <USERNAME> --pw <PASSWORD> --extension 'annotated.coding.txt' # downloading the variant annotations of interest
python main.py --gene_list "['ABCA7']" --outdir './test_output' --username <USERNAME> --pw <PASSWORD> --extension 'annotated.vcf.gz' # downloading the variant call files

```

### You can use the following functions to further process the data

the goal is to have a table of variants for each of the ROSMAP individuals for each variant of interest + an information / annotation table for each variant
and some QC filtering?












