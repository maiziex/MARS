
 :octopus: ***H***ominid ***A***ncest***r***al ***P***opulation analysis (HARP) # :gorrila:
## Install through Bioconda (The latest version 1.0.0):
```
conda install HARP
```
(Please ensure <a href="https://bioconda.github.io/user/install.html#set-up-channels">channels</a> are properly setup for bioconda before installing) 

```
HARP_step1 --help
HARP_step2 --help
# You can also check the below corresponding scripts for more details
```

## Dependencies through Github install:
HARP utilizes <a href="https://www.python.org/downloads/">Python3</a>, <a href="https://github.com/lh3/minimap2/tree/master/misc">paftools (Called haploid assemblies based variants)</a>, <a href="https://tandem.bu.edu/trf/trf.html">trf (Tandem Repeats Finder)</a>, <a href="http://samtools.sourceforge.net/">SAMtools</a>, and <a href="https://github.com/lh3/minimap2">minimap2</a>. To be able to execute the above programs by typing their name on the command line, the program executables must be in one of the directories listed in the PATH environment variable (".bashrc"). <br />
Or you could just run "./install.sh" to install them, but make sure you have installed "conda" and "wget" first. 

## Install through Github:
```
git clone https://github.com/maiziex/HARP.git
cd HARP
chmod +x install.sh
./install.sh
```

## source folder:
After running "./install.sh", a folder "source" would be download, it includes human GRCh38 reference fasta file, and reference fasta files for Gorrila, Orangutan, Chimpanzee, and Macaca which you can use for ancestral call. You could also just download them by yourself from the corresponding official websites. 

## Running The Code:
Put the "HARP/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or use the fullpath of "HARP_step1.py" and "HARP_step2.py"


### Step 1: Assemlby-based structural variants calling for population
```
HARP_step1.py --assembly_dir Aquila_results_30samples --ref_file ./source/genome.fa  --SV_len 20 --num_threads 2 --sample_list 'HG00250','HG00353','HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','HG00851','HG01971','HG02623','HG03115','HG03838','NA12878','NA18552','NA19068','NA19238','NA19239','NA19240','NA19440','NA19789','NA20587','NA24143','NA24149','NA24385','hgp','HLA1','HLA2','HLA3','HLA4','HLA5','HLA7','HLA9','HLA10'  --chr_start 21 --chr_end 21 --out_dir Results_SV_calls
```
#### *Required parameters
##### --assembly_dir: "Aquila_results_30samples" is the input folder where you store the diploid assembled contig files for each sample.  

##### --ref_file: "./source/genome.fa" is the human reference fasta (hg38) file which can be download by running "./install.sh". 

#####  --sample_list: 'HG00250','HG00353','HG00512' are the sample names corresponding to your contig files, which is the prefix of the contig files. 

#### *Optional parameters
#####  --out_dir: default = ./Results_SV_calls, it is the folder name you can define to store the final results.  

#####  --SV_len: default = 20, it is the SV size you can define.

#####  --num_threads: default = 2, it is the number of threads you can define to perform assembly-based variant calling, which corresponds to number of samples.

### Step 2: Generate the population multiple-alignments files for each SV, and SVs with ancestral state 
```
HARP_step2.py  --in_dir Results_SV_calls --assembly_dir Aquila_results_30samples --ref_file ./source/genome.fa  --sample_list 'HG00250','HG00353','HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','HG00851','HG01971','HG02623','HG03115','HG03838','NA12878','NA18552','NA19068','NA19238','NA19239','NA19240','NA19440','NA19789','NA20587','NA24143','NA24149','NA24385','hgp','HLA1','HLA2','HLA3','HLA4','HLA5','HLA7','HLA9','HLA10' --chr_start 21 --chr_end 21 --out_dir Results_MSA_HARP --Ape_ref_list "./source/Gorilla_gorilla_ref.fasta","./source/pan_troglodytes_ref.fasta","./source/pongo_abelii_ref.fasta","./source/macaca_mulatta_ref.fasta" --num_threads 15
```
#### *Required parameters
##### --in_dir: Results_SV_calls is the folder to store SV calling results from step1.
##### --assembly_dir: "Aquila_results_30samples" is the input folder where you store the diploid assembled contig files for each sample.  

##### --ref_file: "./source/genome.fa" is the human reference fasta (hg38) file which can be download by running "./install.sh". 

#####  --sample_list: 'HG00250','HG00353','HG00512' are the sample names corresponding to your contig files, which is the prefix of the contig files. 
##### --Ape_ref_list: "Gorilla_gorilla_ref.fasta", "pan_troglodytes_ref.fasta", "pongo_abelii_ref.fasta", and "macaca_mulatta_ref.fasta" are the reference fasta files for each Ape. Each reference file is seperately by comma (",") 

#### *Optional parameters
#####  --out_dir: default = ./Results_MSA_HARP, it is the folder name you can define to store the final results.  

#####  --SV_len: default = 20, it is the SV size you can define.

#####  --num_threads: default = 15, it is the number of threads you can define to perform multialignment by muscle, which corresponds to number of SV.

## Final Output:

## Cite HARP:
#### 
##### <a href="https://www.biorxiv.org/content/10.1101/660605v1">bioRxiv link</a>


## Troubleshooting:
##### Please submit issues on the github page for <a href="https://github.com/maiziex/HARP/issues">HARP</a>. 
##### Or contact with me through <a href="xzhou15@cs.stanford.edu">xzhou15@cs.stanford.edu</a>




