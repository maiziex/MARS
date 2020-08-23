***M***ultiple ***A***lignment-based ***R***efinement of ***S***vs (MARS)
<p align="center">
	<img src="https://github.com/maiziex/MARS/blob/master/source/msa3.png"  width="600" height="400">
	<p align="center">
		<em></em>
	</p>
</p>

***H***ominid ***A***ncest***r***al ***P***opulation analysis (HARP) mode (set "--HARP_flag" to 1 )
<p align="center">
	<img src="https://github.com/maiziex/MARS/blob/master/source/HARP_icon.png"  width="400" height="200">
	<p align="center">
		<em></em>
	</p>
	
</p>
<p align="center">
	<img src="https://github.com/maiziex/MARS/blob/master/source/msa2.png"  width="500" height="250">
	<p align="center">
		<em></em>
	</p>
</p>


   
## Install through Bioconda (The latest version 1.0.0):
```
conda install MARS
```
(Please ensure <a href="https://bioconda.github.io/user/install.html#set-up-channels">channels</a> are properly setup for bioconda before installing) 

```
MARS_step1 --help
MARS_step2 --help
# You can also check the below corresponding scripts for more details
```

## Dependencies through Github install:
MARS utilizes <a href="https://www.python.org/downloads/">Python3</a>, <a href="https://github.com/lh3/minimap2/tree/master/misc">paftools</a>, <a href="http://samtools.sourceforge.net/">SAMtools</a>, and <a href="https://github.com/lh3/minimap2">minimap2</a>. To be able to execute the above programs by typing their name on the command line, the program executables must be in one of the directories listed in the PATH environment variable (".bashrc"). <br />
Or you could just run "./install.sh" to install them, but make sure you have installed "conda" and "wget" first. 

## Install through Github:
```
git clone https://github.com/maiziex/MARS.git
cd MARS
chmod +x install.sh
./install.sh
```

## source folder:
After running "./install.sh", a folder "source" would be download, and reference fasta files for Gorrila, Orangutan, Chimpanzee, and Macaca which you can use for ancestral call. You could also just download them by yourself from the corresponding official websites. 

## Running The Code:
Put the "MARS/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or use the fullpath of "MARS_step1.py" and "MARS_step2.py"


### Step 1: Assembly-based structural variants calling for population
```
MARS_step1.py  --assembly_dir Aquila_results_30samples --ref_file refdata-GRCh38-2.1.0/fasta/genome.fa  --SV_len 20 --num_threads 2 <br> --sample_list 'HG00250','HG00353','HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','HG00851','HG01971','HG02623','HG03115','HG03838','NA12878','NA18552','NA19068','NA19238','NA19239','NA19240','NA19440','NA19789','NA20587','NA24143','NA24149','NA24385','HGP10X','SL10X1','SL10X2','SL10X3','SL10X4','SL10X5','SL10X7','SL10X9','SL10X10'  --chr_start 1 --chr_end 23 --out_dir Results_SV_calls
```
#### *Required parameters
##### --assembly_dir: "Aquila_results_30samples" is the input folder where you store the diploid assembled contig files for each sample.  

##### --ref_file: "refdata-GRCh38-2.1.0/fasta/genome.fa" is the human reference fasta (hg38) file. 

#####  --sample_list: 'HG00250','HG00353','HG00512' ... are the sample names corresponding to your contig files, which is the prefix of the contig files. 

#### *Optional parameters
#####  --out_dir: default = ./Results_SV_calls, it is the folder name you can define to store the final results.  

#####  --SV_len: default = 20, it is the SV size you can define.

#####  --num_threads: default = 2, it is the number of threads you can define to perform assembly-based variant calling, which corresponds to number of samples.

### Step 2: Merge SV and Generate the multiple samples alignment-based SV files 
```
MARS_step2.py  --in_dir Results_SV_calls --assembly_dir Aquila_results_30samples --ref_file refdata-GRCh38-2.1.0/fasta/genome.fa  --sample_list 'HG00250','HG00353','HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','HG00851','HG01971','HG02623','HG03115','HG03838','NA12878','NA18552','NA19068','NA19238','NA19239','NA19240','NA19440','NA19789','NA20587','NA24143','NA24149','NA24385','HGP10X','SL10X1','SL10X2','SL10X3','SL10X4','SL10X5','SL10X7','SL10X9','SL10X10' --chr_start 1 --chr_end 23 --out_dir Results_MSA_MARS --Ape_ref_list "Gorilla_gorilla_ref.fasta","pan_troglodytes_ref.fasta","pongo_abelii_ref.fasta","macaca_mulatta_ref.fasta"  -gnomad_flag 1 --gnomad_dir gnomAD_hg38_snp --HARP_flag 1 --num_threads_bychr 3
```
#### *Required parameters
##### --in_dir: Results_SV_calls is the folder to store SV calling results from step1.
##### --assembly_dir: "Aquila_results_30samples" is the input folder where you store the diploid assembled contig files for each sample.  

##### --ref_file: "refdata-GRCh38-2.1.0/fasta/genome.fa" is the human reference fasta (hg38) file. 

#####  --sample_list: 'HG00250','HG00353','HG00512'... are the sample names corresponding to your contig files, which is the prefix of the contig files. 

#### *Optional parameters
#####  --out_dir: default = ./Results_MSA_MARS, it is the folder name you can define to store the final results.  

#####  --SV_len: default = 20, it is the SV size you can define.

#####  --num_threads_bychr: default = 2, it is the number of chromosomes you can define to perform MARS parallely.

#####  --gnomad_flag_linked_snp; -gnomad_flag: default = 0. If flag set to 1, MARS will output linked dbSNP for each SV.
#####  --gnomad_dir: If flag set to 1, the users need to download gnomad VCF files and give a path to the folder which stores these gnomad VCF files. 
#####  --HARP_flag: default = 0. If flag set to 1, 
#####  --Ape_ref_list: If --HARP_flag set to 1, the users need to iniatize the Ape reference genomes they want to use. "Gorilla_gorilla_ref.fasta", "pan_troglodytes_ref.fasta", "pongo_abelii_ref.fasta", and "macaca_mulatta_ref.fasta" are the reference fasta files for each Ape. Each reference file is seperately by comma (",") 
#####  --num_threads_bychr: it is the number of threads you can define to perform assembly-based variant calling, which corresponds to number of chromosomes.

## Output files:
1. SV with 10bp left and right flanking regions around breakpoints: a txt file
<p align="center">
	<img src="https://github.com/maiziex/MARS/blob/master/source/msa3.png"  width="600" height="400">
	<p align="center">
		<em></em>
	</p>
</p>
2. SV with around 500bp left and right flanking regions: a MSA html file:
<a href="https://github.com/maiziex/MARS/blob/master/source/chr21_23807653_23807725_del.html">a html file</a>






## Troubleshooting:
##### Please submit issues on the github page for <a href="https://github.com/maiziex/MARS/issues">MARS</a>. 





