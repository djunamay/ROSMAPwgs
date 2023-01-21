# WGS_processing

### ROSMAP WGS data source: https://www.synapse.org/#!Synapse:syn11724057

### programmatic access to ROSMAP WGS data

To download the files of interest:
```python
import synapseclient 
import synapseutils 
 
# specify the following 
path = '/home/gridsan/djuna/github/ROSMAPwgs/raw_data/syn11724057'
username = 
pw = 
targets = set(['1', '2'])
search_files = 'annotated.vcf'

syn = synapseclient.Synapse() 
syn.login(usernam, pw) 

walkedPath = synapseutils.walk(syn, "syn11724057")

filenames = []
for dirpath, dirname, filename in walkedPath:
    filenames.append(filename)
    
annotated = []    
for f in filenames[0]:
    curr = f[0]
    if search_files in curr:
        annotated.append(f)
        

index = [x[0].split('.')[0].split('_')[-1] in targets for x in annotated]
sele = np.array(annotated)[index][:,-1]

with open(path + '/manifest.txt', 'w') as f:
    for i in sele:
        files = synapseutils.syncFromSynapse(syn, i, path = path) 
        f.write(str(files[0]))

```

before run the above need to extract 



### pipeline:

1. From a list of genes, get chromosome numbers or take user-input chromosome numbers
2. specify which type of files want to download for each chromosome; annotated ones or not..
3. download these; QC; annotate; process to get the variants of interest