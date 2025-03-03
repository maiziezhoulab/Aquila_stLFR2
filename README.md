# :milky_way: Aquila_stLFR2 :eagle: 

```
#Download the reference file (hg38)
wget http://xinzhouneuroscience.org/wp-content/uploads/2019/05/source.tar.gz

#Download hg38 "Uniqness_map"
wget http://xinzhouneuroscience.org/wp-content/uploads/2019/05/Uniqness_map.tar.gz
```
## Dependencies for Github installation:
Aquila_stLFR2 utilizes <a href="https://www.python.org/downloads/">Python3 (+ numpy, pysam, sortedcontainers, and scipy)</a>, <a href="http://samtools.sourceforge.net/">SAMtools</a>, and <a href="https://github.com/lh3/minimap2">minimap2</a>. To be able to execute the above programs by typing their name on the command line, the program executables must be in one of the directories listed in the PATH environment variable (".bashrc"). <br />
Or you could just run "./install.sh" to check their availability and install them if not, but make sure you have installed "python3", "conda" and "wget" first. 

# Install through Github:
```
git clone https://github.com/maiziezhoulab/Aquila_stLFR2.git
cd Aquila_stLFR2
chmod +x install.sh
./install.sh
```


After running "./install.sh", two folders "Uniqness_map_hg19" and "Uniqness_map_hg38" would be downloaded, they include the uniqness map information for the human hG19/GRCh38 reference fasta file. 

## Running The Code:
Put the "Aquila_stLFR2/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or just use the fullpath of "**Aquila_stLFR_step1.py**" and "**Aquila_stLFR_step2.py**"

*We provide  <a href="https://github.com/maiziezhoulab/Aquila_stLFR2/blob/master/example_data/run_example_data.md">a small chromosome (chr21) example dataset</a> to run the whole pipeline before you try it into the large dataset. 


### Step 1: 
```
Aquila_stLFR2/bin/Aquila_stLFR_step1.py --fastq_file S12878.fastq --bam_file S12878.bam --vcf_file S12878_freebayes.vcf --sample_name S12878 --out_dir Assembly_results_S12878 --uniq_map_dir Aquila_stLFR2/Uniqness_map_hg19
```
#### *Required parameters
**--fastq_file:** "S12878.fastq" is the stLFR fastq file (with BX:Z:barcode at the header, you can use Aquila_stLFR/bin/Aquila_stLFR_fastq_preprocess.py to generate the input fastq file, <a href="https://github.com/maiziezhoulab/Aquila_stLFR2/blob/master/src/How_to_get_bam_and_vcf.md">check here for the processing details</a>)

**--bam_file:** "S12878.bam" is bam file generated from bwa-mem. How to get bam file, you can also check <a href="https://github.com/maiziezhoulab/Aquila_stLFR2/blob/master/src/How_to_get_bam_and_vcf.md">here</a>.

**--vcf_file:** "S12878_freebayes.vcf" is VCF file generated from variant caller like "FreeBayes". How to get vcf file, you can also check <a href="https://github.com/maiziezhoulab/Aquila_stLFR2/blob/master/src/How_to_get_bam_and_vcf.md">here</a>. 

**--sample_name:** "S12878" are the sample name you can define. 

**--uniq_map_dir:** "Aquila_stLFR2/Uniqness_map_hg19" is the uniqness file you can download by "./install.sh".

Note: use Uniqness_map_hg19 for BAM based on HG19 reference, Uniqness_map_hg38 for BAM based on HG38 reference

**--chr_number :** target chromosome. For example: using "--chr_number 1 "  will assemble chromosome 1. 

**--read_type :** PE for paired-end reads, SE for single-end reads 

#### *Optional parameters
**--out_dir:** default = ./Asssembly_results. You can define your own folder, for example "Assembly_results_S12878". 

**--block_threshold:** default = 200000 (200kb)
 
**--block_len_use:** default = 100000 (100kb)

**--num_threads:** default = 8. It's recommended not to change this setting unless large memory node could be used (2*memory capacity(it suggests for assembly below)), then try to use "--num_threads 12". 




#### Memory/Time Usage For Step 1
##### Running Step 1 for chromosomes parallelly on multiple(23) nodes

Coverage | Memory| Time for chr1 on a single node | 
--- | --- | --- | 
90X | 350GB | 1-10:20:10 |

Coverage | Memory| Time for chr21 on a single node | 
--- | --- | --- | 
90X | 100GB | 06:27:28 |





### Step 2: 
```
Aquila_stLFR2/bin/Aquila_stLFR_step2.py --out_dir Assembly_results_S12878 --num_threads 30 --reference <your reference file>
```
#### *Required parameters
**--reference:** reference file.
**--chr_number :** target chromosome. For example: using "--chr_number 1 "  will assemble chromosome 1. 
**--read_type :** PE for paired-end reads, SE for single-end reads 

#### *Optional parameters
**--out_dir:** default = ./Asssembly_results, make sure it's the same as "--out_dir" from ***Step1*** if you want to define your own output directory name.

**--num_threads:** default = 30, this determines the number of files assembled simultaneously by SPAdes.  

**--num_threads_spades:** default = 5, this is the "-t" for SPAdes. 

**--block_len_use:** default = 100000 (100kb)


### Step 3: 
```
Aquila_stLFR2/bin/Aquila_stLFR_step3.py --out_dir Assembly_results_S12878/ --num_threads 30 --reference <your reference file>
```
#### *Required parameters
**--reference:** reference file.
**--chr_number :** target chromosome. For example: using "--chr_number 1 "  will assemble chromosome 1. 
**--out_dir:** default = ./Asssembly_results, make sure it's the same as "--out_dir" from ***Step1*** and ***Step2*** if you want to define your own output directory name.




#### Memory/Time Usage For Step 2
##### Running Step 2 for chromosomes parallelly on multiple nodes
Coverage| Memory| Time for chr1 on a single node | --num_threads | --num_threads_spades|
--- | --- | --- | ---|---|
90X| 100GB | 11:13:19 |30 | 20|

Coverage| Memory| Time for chr21 on a single node | --num_threads | --num_threads_spades|
--- | --- | --- | ---|---|
90X| 100GB | 01:38:09 |30 | 20|



## Final Output:
**Assembly_Results_S12878/Assembly_Contigs_files:** Aquila_contig.fasta and Aquila_Contig_chr*.fasta 
```
Assembly_results_S12878
|
|-H5_for_molecules 
|   └-S12878_chr*_sorted.h5    --> (Fragment files for each chromosome including barcode, variants annotation (0: ref allele; 1: alt allele), coordinates for each fragment)
|
|-HighConf_file
|   └-chr*_global_track.p      --> (Pickle file for saving coordinates of high-confidence boundary points)
|
|-results_phased_probmodel
|   └-chr*.phased_final        --> (Phased fragment files)
|
|-phase_blocks_cut_highconf
|
|-Raw_fastqs
|   └-fastq_by_Chr_*           --> (fastq file for each chromosome)
|
|-ref_dir
|
|-Local_Assembly_by_chunks
|   └-chr*_files_cutPBHC
|       |-fastq_by_*_*_hp1.fastq                  --> (fastq file for a small phased chunk of haplotype 1)
|       |-fastq_by_*_*_hp2.fastq                  --> (fastq file for a small phased chunk of haplotype 2)
|       |-fastq_by_*_*_hp1_spades_assembly        --> (minicontigs: assembly results for the small chunk of haplotype 1) 
|       └-fastq_by_*_*_hp2_spades_assembly        --> (minicontigs: assembly results for the small chunk of haplotype 2)
|
└-Assembly_Contigs_files
    |-Aquila_cutPBHC_minicontig_chr*.fasta        --> (final minicontigs for each chromosome)
    |-Aquila_Contig_chr*.fasta                    --> (final contigs for each chromosome)
    └-Aquila_contig.fasta                         --> (final contigs for WGS)
```

## Final Output Format:
Aquila_stLFR outputs an overall contig file `Aquila_Contig_chr*.fasta` for each chromosome, and one contig file for each haplotype: `Aquila_Contig_chr*_hp1.fasta` and `Aquila_Contig_chr*_hp2.fasta`. For each contig, the header, for an instance, “>36_PS39049620:39149620_hp1” includes contig number “36”, phase block start coordinate “39049620”, phase block end coordinate “39149620”, and haplotype number “1”. Within the same phase block, the haplotype number “hp1” and “hp2” are arbitrary for maternal and paternal haplotypes. For some contigs from large phase blocks, the headers are much longer and complex, for an instance, “>56432_PS176969599:181582362_hp1_ merge177969599:178064599_hp1-177869599:177969599_hp1”. “56” denotes contig number, “176969599” denotes the start coordinate of the final big phase block, “181582362” denotes the end coordinate of the final big phase block, and “hp1” denotes the haplotype “1”. “177969599:178064599_hp1” and “177869599:177969599_hp1” mean that this contig is concatenated from minicontigs in small chunk (start coordinate: 177969599, end coordinate: 178064599, and haplotype: 1) and small chunk (start coordinate: 177869599, end coordinate: 177969599, and haplotype: 1). 

### Clean Data
##### If your hard drive storage is limited(Aquila_stLFR could generate a lot of intermediate files for local assembly), it is suggested to quily clean some data by running `Aquila_stLFR_clean.py`. Or you can keep them for some analysis (check the above output directory tree for details). 
```
Aquila_stLFR2/bin/Aquila_stLFR_clean.py --assembly_dir Assembly_results_S12878 
```

## Assembly Based Variants Calling and Phasing:
##### For example, you can use `Assemlby_results_S12878` as input directory to generate a VCF file which includes SNPs, small Indels and SVs. 
##### Please check <a href="https://github.com/maiziezhoulab/Aquila_stLFR2/blob/master/Assembly_based_variants_call/README.md/">Assembly_based_variants_call_and_phasing</a> for details. 



### Notes
#### For stLFR assembly, stLFR reads with barcode "0_0_0" are removed to get perfect diploid assembly.  


#### Please also check <a href="https://github.com/maiziezhoulab/Aquila">Aquila</a> here 


## Troubleshooting:
##### Please submit issues on the github page for <a href="https://github.com/maiziezhoulab/Aquila_stLFR2/issues">Aquila_stLFR2</a>. 
##### Or contact with me through <a href="can.luo@vanderbilt.edu">can.luo@vanderbilt.edu</a>

