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


####Filter 1 exclude SNPs round INDELs
vcftools --vcf SNP_teff.Dabbi_50954_V3_kora_asgori.vcf --out filtered_indel_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf --exclude-bed snp_excluson_data.txt














# ####Fitler vcf 
# #For fitlering bcftools uses "&" to keep vairiants that meet all given criteria and "&&" to keep variants that match at least one of the given criteria



# bcftools filter -O v -o filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10 & INFO/DP <= 320 & MQBZ >=-2 && MQBZ <= 2 & RPBZ >= -7.5 && RPBZ <= 7.5 & SCBZ > 4' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf


# https://www.htslib.org/workflow/filter.html

# bcftools view -i 'QUAL >= 10 || DP <= 320 || MQBZ >=-2 || MQBZ <= 2 || RPBZ <= -7.5 || RPBZ >= 7.5 || SCBZ > 4' in.vcf

# ####Filter1 - based on quality, a score of >20 means that there is a 99% probability there is a variant at that site 



# SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
# 1019991 
# bcftools filter -O v -o Qual_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
# 994383
# bcftools filter -O v -o Qual_Depth_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'INFO/DP <= 320' Qual_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
# 975206
# bcftools filter -O v -o Qual_Depth_mqbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'MQBZ >=-2 && MQBZ <= 2'Qual_Depth_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
# 72387
# bcftools filter -O v -o Qual_Depth_mqbz_rpbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'RPBZ >= -7.5 && RPBZ <= 7.5' Qual_Depth_mqbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
# 72212
# bcftools filter -O v -o Qual_Depth_mqbz_rpbz_scbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'SCBZ > 4' Qual_Depth_mqbz_rpbz_filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
# 55



# bcftools filter -O v -o filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10 & INFO/DP <= 320 & MQBZ >=-5 && MQBZ <= 5 && MQSBZ >=-5 && MQSBZ <= 5 & RPBZ >= -7.5 && RPBZ <= 7.5 & SCBZ > 4' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf

# bcftools filter -O v -o 2filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10 & INFO/DP <= 320 & MQSBZ >=-2 && MQSBZ <= 2 & RPBZ >= -7.5 && RPBZ <= 7.5 & SCBZ > 4' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf



