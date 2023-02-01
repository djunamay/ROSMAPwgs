## ROSMAPwgs:
This repository contains code to download [whole genome sequencing data from the ROSMAP longitudinal cohort study](https://www.synapse.org/#!Synapse:syn11724057) and extract coding variant annotation information for genes of interest. For more information see [here](https://www.synapse.org/#!Synapse:syn10901595).

## Getting started: 

1. Set up environment
```bash
conda create --name wgs_env
conda install -c bioconda tabix 
pip install scikit-allel
```

2. Download ensembl release GRCh37 annotations & generate dictionaries
```bash
git clone git@github.com:djunamay/ROSMAPwgs.git
mkdir ./human_Release_19_GRCh37p13
cd human_Release_19_GRCh37p13
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz
python ./ROSMAPwgs/get_dictionaries.py
```
*N.B. Be sure to download __GRCh37__ as "Paired-end 150bp reads were aligned to the GRCh37 human reference" -- see the documentation on [synapse](https://www.synapse.org/#!Synapse:syn10901595)*


3. Download vcf files for your genes of interest
```bash
python main.py --outdir './raw_data' --username <USERNAME> --pw <PASSWORD> --gene_list "['APOE', 'ABCA1']" --extension 'recalibrated_variants.vcf.gz' --extract_HIGHandMED_annotations False --download True
```
*N.B. These files are available to download from synapse:*

| File type  | Description |
| ------------- | ------------- |
| annotated.clinical.txt  | All low frequency HIGH/MODERATE annotated variants with possible clinical impact (from ClinVar) in text file format  |
| annotated.coding_rare.txt  | All HIGH/MODERATE annotated variants with less than 5% allele frequency in 1000genomes and ExAC in text file format  |
| annotated.coding.txt  | All annotated variants with HIGH/MODERATE impact in text file format  |
| annotated.txt  | All variants with annotations in text file format  |
| annotated.vcf.gz  | All variants with annotations in variant call format (VCF)  |
| annotated.vcf.gz.tbi  | Index file for annotated VCF file  |
| recalibrated_variants.vcf.gz  | All variants in variant call format  |
| recalibrated_variants.vcf.gz.tbi  | Index file for VCF file with all variants  |

	
4. Extract High/Moderate variant annotations for genes of interest
```bash
python main.py --outdir './raw_data' --username <USERNAME> --pw <PASSWORD> --gene_list "['APOE', 'ABCA1']" --extension 'annotated.coding.txt' --extract_HIGHandMED_annotations False --download True
python main.py --outdir './raw_data' --gene_list "['APOE', 'ABCA1']" --extract_HIGHandMED_annotations True --download False
```

*N.B.*
- For now, this pipeline only extracts variant info for coding regions of interest (by gene) - it just takes a few modifications to generalize this to all genomic regions of interest. We'll add that soon!
- If `--download` is `True`, `--outdir` must contain the following files for each chromosome of interest: recalibrated_variants.vcf.gz, recalibrated_variants.vcf.gz.tbi, recalibrated_variants.annotated.coding.txt. These files can be downloaded to the `--outdir` by specifying --`skip_download` as `False` and `--extension` as recalibrated_variants.vcf.gz
- Make sure the genes you specify in `--gene_list` are located on the chromosomes for which you have files in the `--outdir`
- This pipeline intentionally allows separation of the download and annotation extraction components. This is useful if you are running this code on a server and only only have internet access on the login node. In this case, run the first command on the login node and the second command on a compute node. Otherwise, just execute together with `--download True` and `--extract_HIGHandMED_annotations True`

5. To run tests:
```bash
python main.py --outdir './test_output' --username <USERNAME> --pw <PASSWORD> --gene_list "['APOE']" --extension 'recalibrated_variants.vcf.gz' --extract_HIGHandMED_annotations False --download True

pip install assertpy
pip install pytest
python -m pytest test.py
```
*N.B. If you have the specified paths path1-4 (see test.py) in the `./test_output` directory, additional tests will be performed and the testing will be slower. Remove these files or rename the parent folder if you don't want the additional tests to be run.*













