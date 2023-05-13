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
bcftools index /home/jconnell/eragrostis/teff.Dabbi_50954_V3_kora_asgori_dartseq.vcf