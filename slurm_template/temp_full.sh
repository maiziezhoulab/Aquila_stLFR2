


#################################################

#           Parameters

#################################################

chr_num=6
reference=/data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa
bamfile=/lio/lfs/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350158968/align_EMA/V350158968_hg19_ema.bam 
dtype=PE
outdir=./Assembly_results_NA24385_chr${chr_num}/
fastq=   # single chromosome fastq
vcf=/data/maiziezhou_lab/ZhouLab_Projects/VandyCG_stLFR/Experiments/Constant_Use/Real_Data/V350158968/align_by_ema/gatk/gatk.vcf
t=30


/data/maiziezhou_lab/CanLuo/Software/Aquila_stLFR/bin/Aquila_stLFR_step1_test.py \
	--fastq_file $fastq \
	--bam_file  $bamfile  \
	--vcf_file $vcf \
	--sample_name NA24385 \
	--out_dir  $outdir \
	--uniq_map_dir /data/maiziezhou_lab/Maizie/Softwares/Aquila/Uniqness_map_hg19 \
	--chr_start ${chr_num}  --chr_end ${chr_num} --num_threads $t -r $dtype


python3 /data/maiziezhou_lab/CanLuo/Software/Aquila_stLFR/bin/Extract_qname_from_phased_molecule_cut_phase_blocks_v5.py \
	-i ./Assembly_results_NA24385_chr${chr_num}/ \
	-fq ${fastq} -rt PE -chr ${chr_num} -t 20 



/data/maiziezhou_lab/CanLuo/Software/Aquila_stLFR/bin/Aquila_stLFR_step2.py \
	--out_dir ${outdir} \
	--num_threads $t \
	--reference $reference \
	--chr_start ${chr_num} --chr_end ${chr_num} -r $dtype

#
####### original code vc and phasing 
#
#python3 /data/maiziezhou_lab/CanLuo/Software/Aquila_stLFR/bin/Aquila_stLFR_assembly_based_variants_call.py \
	#    --assembly_dir ${outdir} \
	#    --out_dir ${outdir}/Variant_Call \
	#    --ref_file ${reference} \
	#    --chr_start ${chr_num} --chr_end ${chr_num}

python3 /data/maiziezhou_lab/CanLuo/long_reads_project/bin/SVIM-asm_variant_caller_eval.py \
	-i ${outdir}/Assembly_Contigs_files/Aquila_contig.fasta \
	-o ${outdir}/Variant_Call_svim-asm/  -chr ${chr_num}

python3 /data/maiziezhou_lab/CanLuo/Software/Aquila_stLFR/volcanosv-vc/dippav_vc_sr.py \
	-contig ./Assembly_results_NA24385_chr${chr_num}/Assembly_Contigs_files/Aquila_Contig_chr${chr_num}.fasta \
	-ref /lio/lfs/maiziezhou_lab/maiziezhou_lab/CanLuo/long_reads_project/DipPAV2/hg19_ref_by_chr/hg19_chr${chr_num}.fa \
	-o volcanosv_chr${chr_num} \
	-chr ${chr_num} -t 10 

# python3 /data/maiziezhou_lab/CanLuo/Software/Aquila_stLFR/volcanosv-vc/check_del_fp.py -v volcanosv_chr6/final_vcf/dippav_variant_no_redundancy_filtered.vcf -r 30-100 -e /data/maiziezhou_lab/ZhouLab_Projects/VandyCG_stLFR/Experiments/2023_October/Aquila_stLFR_chr6/volcanosv_chr6/final_vcf/eval/Truvari_NGMLR_aligned_results_p0.1_P0.1_dist200_S0_O0_chr6/DEL_50_/ -o del_grid_search.csv


cov=`samtools depth  /lio/lfs/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350158968/align_EMA/V350158968_hg19_ema.bam  -r chr${chr_num}:20000000-40000000 |  awk '{sum+=$3} END { print sum/NR}'`


python3 /data/maiziezhou_lab/CanLuo/Software/Aquila_stLFR/volcanosv-vc/sv_filter.py \
	-v volcanosv_chr${chr_num}/final_vcf/dippav_variant_no_redundancy.vcf \
	-b $bamfile  \
	-cov $cov -t 30 -r 0.84

/data/maiziezhou_lab/CanLuo/long_reads_project/bin/truvari_eval.sh \
	${chr_num} ./volcanosv_chr${chr_num}/final_vcf/ \
	./volcanosv_chr${chr_num}/final_vcf/eval_ft dippav_variant_no_redundancy_filtered 500 0.5 0.5 30 0.01










