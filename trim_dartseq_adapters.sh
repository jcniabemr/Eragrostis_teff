#!/usr/bin/env bash
#SBATCH --partition medium 
#SBATCH -J trimgalore
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --mail-user=john.connell@niab.com
#SBATCH --mail-type=END,FAIL

####Trim adapter/collect stats 
####Larry exmaple #trim_galore -q 25 --length 75 --paired /filepath 
####Larry python script 
trim_barcodes=/mnt/shared/scratch/lperciva/allvall_JC/fix_barcodes/fix_barcodes.py
####Ref genome 
ref_genome=/home/jconnell/eragrostis/ref_genome/teff.Dabbi_50954_V3.fasta

####Activate env
source activate trimgalore

####Trim illumina seq adapter only and map to ref
# for x in $(ls /home/jconnell/eragrostis/2023_camb_teff/Dartseq_results/FASTQ_files/ | grep -v "targets_HJ7WCDMXY_2.csv"); do
# 	name=$(basename ${x} .FASTQ.gz)
# 	outdir=/home/jconnell/niab/eragrostis/trimmed_reads/${name}
# 	mkdir -p $outdir
# 	trim_galore -q 25 --length 50 --output_dir $outdir /home/jconnell/eragrostis/2023_camb_teff/Dartseq_results/FASTQ_files/${x}
# 	bwa mem -t 4 -M $ref_genome $outdir/"$name".FASTQ.gz_trimmed.fq.gz | samtools view -b - | samtools sort - -O 'BAM' -o $outdir/"$name"_alignment_sorted_1trim.bam
# done 

###Trim barcode+cutsite+seq adapter and map to ref
for x in $(ls /home/jconnell/eragrostis/2023_camb_teff/Dartseq_results/FASTQ_files/ | grep -v "targets_HJ7WCDMXY_2.csv"); do
	name=$(basename ${x} .FASTQ.gz)
	outdir=/home/jconnell/niab/eragrostis/trimmed_reads_barcodes/${name}
	mkdir -p $outdir
	python $trim_barcodes /home/jconnell/eragrostis/2023_camb_teff/Dartseq_results/FASTQ_files/${x} > $outdir/"$name"_inline-bc_trim.fastq
	trim_galore -q 25 --length 50 --output_dir $outdir $outdir/"$name"_inline-bc_trim.fastq
	bwa mem -t 4 -M $ref_genome $outdir/"$name"_inline-bc_trim_trimmed.fq | samtools view -b - | samtools sort - -O 'BAM' -o $outdir/"$name"_alignment_sorted_2trim.bam
done 

####Call variants 
conda deactivate 
source activate bcftools 

bcftools mpileup \
-Ou \
--bam-list <(for x in $(ls /home/jconnell/eragrostis/trimmed_reads_barcodes/*/*_alignment_sorted_2trim.bam); do echo $x; done) \
-q 10 \
-C 50 \
-a AD,DP \
-f ${ref_genome} | bcftools call -c -v -Ov > /home/jconnell/eragrostis/teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf 

####Replace bam file headders with two sed commands to strip /home/jconnell/eragrostis/trimmed_reads_barcodes/ and 3160705_alignment_sorted_2trim.bam leaving on the ID
cat /home/jconnell/eragrostis/teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf | sed 's/\/home\/jconnell\/eragrostis\/trimmed_reads_barcodes\///g; s/\/\([0-9]*\)_alignment_sorted_2trim\.bam//g' > edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf 
####Remove indels 
cat edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf | grep "#\|INDEL" > INDEL_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf
cat edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf | grep -v "INDEL" > SNP_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf
####Create bed file for INDEL SNP filter and depth/qual data for plotting 
python /home/jconnell/git_repos/niab_repos/eragrostis/extract_vcf_info.py -snp SNP_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf -indel INDEL_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf
####Exclude SNPs round INDELs
vcftools --vcf SNP_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf --exclude-bed snp_excluson_data.txt --recode --recode-INFO-all --out filtered_indel_SNP_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf
###Filter based on a quality score of >20 == 99% probability there is a variant at that site. 
bcftools filter -O v -o qual_filtered_indel_SNP_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf.recode.vcf -i 'QUAL >= 20' filtered_indel_SNP_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf.recode.vcf
####Remove bialelic SNP
cat qual_filtered_indel_SNP_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf.recode.vcf | awk -F '\t' '{split($5,a,","); if(length(a)<=1) print}' > bi_qual_filtered_indel_SNP_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq_recode.vcf
####Take lines were data is in both files  
python /home/jconnell/git_repos/niab_repos/eragrostis/find_similar_variants.py --dart bi_qual_filtered_indel_SNP_edited_teff.Dabbi_50954_V3_kora_asgori_dartseq_recode.vcf --seq unique_Asgori-Tank1_Kora-Tank1_filtered_combined.vcf.recode.vcf
####Merge dart and Kora_Asgori data 
####Before this you need to make sure the ## fields contain the contig and format fields, here this was just coppied from the other VCF 
bcftools sort Dart_common_SNP.vcf -Oz -o Dart_common_SNP.vcf.gz
bcftools sort Kora_Asgori_common_SNP.vcf -Oz -o Kora_Asgori_common_SNP.vcf.gz
####Create index
bcftools index Dart_common_SNP.vcf.gz
bcftools index Kora_Asgori_common_SNP.vcf.gz
bcftools merge Kora_Asgori_common_SNP.vcf.gz Dart_common_SNP.vcf.gz > combined_final_VCFDART_VCF_GS.vcf