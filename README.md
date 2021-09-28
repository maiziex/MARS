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
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/mars.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/mars)  
```
conda install MARS
```
(Please ensure <a href="https://bioconda.github.io/user/install.html#set-up-channels">channels</a> are properly setup for bioconda before installing) 

```
MARS_step1 --help
MARS_step2 --help
# You can also check the below corresponding scripts for more details
```

## Install through Github:
```
git clone https://github.com/maiziex/MARS.git
cd MARS
chmod +x install.sh
./install.sh
```

## Dependencies through Github install:
MARS utilizes <a href="https://www.python.org/downloads/">Python3</a>, <a href="https://github.com/lh3/minimap2/tree/master/misc">paftools</a>, <a href="http://samtools.sourceforge.net/">SAMtools</a>, and <a href="https://github.com/lh3/minimap2">minimap2</a>. To be able to execute the above programs by typing their name on the command line, the program executables must be in one of the directories listed in the PATH environment variable (".bashrc"). <br />
Or you could just run "./install.sh" to install them, but make sure you have installed "conda" and "wget" first. 

To use MARS, please also make sure install the below python packages through pip or conda. \
pip install Bio \
Conda install pysam\
Conda install pandas\
conda install minimap2\
conda install openpyxl 

## Example data in Zenodo:
Please download the example data <a href="https://zenodo.org/record/5512900">from Zenodo</a>.
```
Ape_ref
|-Gorilla_gorilla_ref.fasta
|-pongo_abelii_ref.fasta
|-pan_troglodytes_ref.fasta
|-macaca_mulatta_ref.fasta

refdata-GRCh38-2.1.0
|-fasta
     |-genome.fa     
     
gnomAD_hg38_snp
|-gnomad.genomes.r2.1.1.sites.22.liftover_grch38.vcf.bgz
     
Aquila_results_4samples
|-HG00353
|      |-Assembly_Contigs_files
|                |-Aquila_Contig_chr22.fasta
|-HG00250
|      |-Assembly_Contigs_files
|                |-Aquila_Contig_chr22.fasta			 
|-HG00513
|      |-Assembly_Contigs_files
|                |-Aquila_Contig_chr22.fasta
|-HG00512
      |-Assembly_Contigs_files
                 |-Aquila_Contig_chr22.fasta			 		 
		 
```

## Running The Code:
Put the "MARS/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or use the fullpath of "MARS_step1.py" and "MARS_step2.py"


### Step 1: Assembly-based structural variants calling for population
```
MARS_step1.py  --assembly_dir Aquila_results_4samples --ref_file refdata-GRCh38-2.1.0/fasta/genome.fa  --SV_len 20  --sample_list 'HG00250','HG00353','HG00512','HG00513' --chr_start 22 --chr_end 22 --out_dir MARS_step1_results
```
#### *Required parameters
**--assembly_dir:** "Aquila_results_4samples" is the input folder where you store the diploid assembled contig files for each sample by Aquila/Aquila_stLFR.  

**--ref_file:** "refdata-GRCh38-2.1.0/fasta/genome.fa" is the human reference fasta (hg38) file. 

**--sample_list:** 'HG00250','HG00353','HG00512' ... are the sample names corresponding to your contig files, which is the prefix of the contig files. 

#### *Optional parameters
**--out_dir:** default = ./MARS_step1_results, it is the folder name you can define to store the final results.  

**--SV_len:** default = 20, it is the SV size you can define.

**--num_threads:** default = 2, it is the number of threads you can define to perform assembly-based variant calling, which corresponds to number of samples.

**--chr_start --chr_end:** For example: use "--chr_start 1 --chr_end 5"  will perform results for chromsomes 1,2,3,4,5. Use "--chr_start 22 --chr_end 22" will only assembly chromosome 22. 


### Step 2: Merge SV and Generate the multiple samples alignment-based SV files 
```
MARS_step2.py  --in_dir MARS_step1_results --assembly_dir Aquila_results_4samples --ref_file refdata-GRCh38-2.1.0/fasta/genome.fa  --sample_list 'HG00250','HG00353','HG00512','HG00513' --chr_start 22 --chr_end 22 --out_dir MARS_step2_results 
```

OR add "-gnomad_flag" 

```
MARS_step2.py  --in_dir MARS_step1_results --assembly_dir Aquila_results_4samples --ref_file refdata-GRCh38-2.1.0/fasta/genome.fa  --sample_list 'HG00250','HG00353','HG00512','HG00513' --chr_start 22 --chr_end 22 --out_dir MARS_step2_results -gnomad_flag 1 --gnomad_dir gnomAD_hg38_snp 
```

OR add "--HARP_flag" 

```
MARS_step2.py  --in_dir MARS_step1_results --assembly_dir Aquila_results_4samples --ref_file refdata-GRCh38-2.1.0/fasta/genome.fa  --sample_list 'HG00250','HG00353','HG00512','HG00513' --chr_start 22 --chr_end 22 --out_dir MARS_step2_results  --HARP_flag 1 --Ape_ref_list "Ape_ref/Gorilla_gorilla_ref.fasta","Ape_ref/pan_troglodytes_ref.fasta","Ape_ref/pongo_abelii_ref.fasta","Ape_ref/macaca_mulatta_ref.fasta"  
```

OR add "--HARP_flag" and "-gnomad_flag"

```
MARS_step2.py  --in_dir MARS_step1_results --assembly_dir Aquila_results_4samples --ref_file refdata-GRCh38-2.1.0/fasta/genome.fa  --sample_list 'HG00250','HG00353','HG00512','HG00513' --chr_start 22 --chr_end 22 --out_dir MARS_step2_results  --HARP_flag 1 --Ape_ref_list "Ape_ref/Gorilla_gorilla_ref.fasta","Ape_ref/pan_troglodytes_ref.fasta","Ape_ref/pongo_abelii_ref.fasta","Ape_ref/macaca_mulatta_ref.fasta"  -gnomad_flag 1 --gnomad_dir gnomAD_hg38_snp 
```

#### *Required parameters
**--in_dir:** "MARS_step1_results" is the folder to store SV calling results from step1.

**--assembly_dir:** "Aquila_results_4samples" is the input folder where you store the diploid assembled contig files for each sample.  

**--ref_file:** "refdata-GRCh38-2.1.0/fasta/genome.fa" is the human reference fasta (hg38) file. 

**--sample_list:** 'HG00250','HG00353','HG00512'... are the sample names corresponding to your contig files, which is the prefix of the contig files. 

#### *Optional parameters
**--out_dir:** default = ./MARS_step2_results, it is the folder name you can define to store the final results.  

**--SV_len:** default = 20, it is the SV size you can define.

**--num_threads_bychr:** default = 2, it is the number of chromosomes you can define to perform MARS parallely.

**--chr_start --chr_end:** For example: use "--chr_start 1 --chr_end 5"  will perform results for chromsomes 1,2,3,4,5. Use "--chr_start 22 --chr_end 22" will only assembly chromosome 22. 

**--gnomad_flag_linked_snp; -gnomad_flag:** default = 0. If flag set to 1, MARS will output linked dbSNP from gnomad for each SV.

**--gnomad_dir:** If "gnomad_flag" set to 1, the users need to download gnomad VCF files and give a path to the folder which stores these gnomAD VCF files. (For example: use "wget https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.10.liftover_grch38.vcf.bgz" to download gnomAD VCF file for chr10) For details, please check <a href="https://gnomad.broadinstitute.org/downloads/">gnomAD downloads website</a>) 

**--HARP_flag:** default = 0. If flag set to 1, it will output SV with ancestral state: derived deletions and derived insertions. 

**--Ape_ref_list:** If --HARP_flag set to 1, the users need to iniatize the Ape reference genomes they want to use. "Gorilla_gorilla_ref.fasta", "pan_troglodytes_ref.fasta", "pongo_abelii_ref.fasta", and "macaca_mulatta_ref.fasta" are the reference fasta files for each Ape. Each reference file is seperately by comma (",") 


## Output files:
#### 1. In the "MARS_step2_results/Final_tables" folder:
1. SV_msa_table_chr*.txt  
2. derived_del_msa_table_chr*.txt (set "--HARP_flag = 1")
3. derived_ins_msa_table_chr*.txt (set "--HARP_flag = 1")
4. derived_complex_msa_table_chr*.txt  (set "--HARP_flag = 1")
5. SV_linked_with_gnomad_SNP_chr*.txt (set "-gnomad_flag = 1")
6. Consensus_seq_for_SV_chr*.fasta


#### 2. In the "MARS_step2_results/MSA_SV_results_chr*/MSA_SV_files" folder, it includes txt/html files for each SV:
1. SV with 10bp left and right flanking regions around breakpoints, for example a txt file called by "chr22_29430067_29430067_ins_SV_Ref_bk_2.txt",
"chr22_29430067_29430067_ins" is the name of the SV which corresponds to the first field/column in "SV_msa_table_chr*.txt" file.
<p align="center">
	<img src="https://github.com/maiziex/MARS/blob/master/source/msa3.png"  width="600" height="400">
	<p align="center">
		<em></em>
	</p>
</p>
2. SV with around 500bp left and right flanking regions: a MSA html file:
<a href="https://github.com/maiziex/MARS/blob/master/source/chr22_23530481_23530481_ins.html">a html file</a>






## Troubleshooting:
##### Please submit issues on the github page for <a href="https://github.com/maiziex/MARS/issues">MARS</a>. 





