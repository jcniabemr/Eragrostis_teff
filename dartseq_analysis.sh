#!/usr/bin/env bash 
#SBATCH -J dartseq_eragrostis
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4

####SNP call file 
snpcalls=/home/jconnell/niab/eragrostis/snpcalls

####Convert bcf to vcf 
bcftools convert \
-Ov \
-o ${snpcalls}/teff.Dabbi_50954_V3_kora_asgori.vcf \
${snpcalls}/teff.Dabbi_50954_V3_kora_asgori.bcf


####Sep INDELs from SNPs 

es teff.Dabbi_50954_V3_kora_asgori.vcf | grep "#\|INDEL" > INDEL_teff.Dabbi_50954_V3_kora_asgori.vcf
es teff.Dabbi_50954_V3_kora_asgori.vcf | grep -v "INDEL" > SNP_teff.Dabbi_50954_V3_kora_asgori.vcf


####Fitler vcf 
#For fitlering bcftools uses "&" to keep vairiants that meet all given criteria and "&&" to keep variants that match at least one of the given criteria


ALT=<ID=*,Description="Represents allele(s) other than observed.">
INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
INFO=<ID=RPBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)">
INFO=<ID=MQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)">
INFO=<ID=BQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)">
INFO=<ID=MQSBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)">
INFO=<ID=SCBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)">
INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">
FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
INFO=<ID=AF2,Number=1,Type=Float,Description="Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)">
INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
bcftools_callVersion=1.15+htslib-1.15
bcftools_callCommand=call -c -v -Ob; Date=Thu Mar  3 12:06:06 2022
bcftools_viewVersion=1.9+htslib-1.9
bcftools_viewCommand=view --header-only teff.Dabbi_50954_V3_kora_asgori.vcf; Date=Sat Apr 15 10:42:08 2023

bcftools filter -O v -o filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10 & INFO/DP <= 320 & MQBZ >=-2 && MQBZ <= 2 & RPBZ >= -7.5 && RPBZ <= 7.5 & SCBZ > 4' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf


https://www.htslib.org/workflow/filter.html

bcftools view -i 'QUAL >= 10 || DP <= 320 || MQBZ >=-2 || MQBZ <= 2 || RPBZ <= -7.5 || RPBZ >= 7.5 || SCBZ > 4' in.vcf

####Filter1 - based on quality, a score of >20 means that there is a 99% probability there is a variant at that site 







SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
1019991 
bcftools filter -O v -o Qual_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
994383
bcftools filter -O v -o Qual_Depth_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'INFO/DP <= 320' Qual_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
975206
bcftools filter -O v -o Qual_Depth_mqbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'MQBZ >=-2 && MQBZ <= 2'Qual_Depth_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
72387
bcftools filter -O v -o Qual_Depth_mqbz_rpbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'RPBZ >= -7.5 && RPBZ <= 7.5' Qual_Depth_mqbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
72212
bcftools filter -O v -o Qual_Depth_mqbz_rpbz_scbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'SCBZ > 4' Qual_Depth_mqbz_rpbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
55






bcftools filter -O v -o filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10 & INFO/DP <= 320 & MQBZ >=-5 && MQBZ <= 5 && MQSBZ >=-5 && MQSBZ <= 5 & RPBZ >= -7.5 && RPBZ <= 7.5 & SCBZ > 4' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf

bcftools filter -O v -o 2filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10 & INFO/DP <= 320 & MQSBZ >=-2 && MQSBZ <= 2 & RPBZ >= -7.5 && RPBZ <= 7.5 & SCBZ > 4' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
