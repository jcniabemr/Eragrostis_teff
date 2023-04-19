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

# #For fitlering bcftools uses "&" to keep vairiants that meet all given criteria and "&&" to keep variants that match at least one of the given criteria
# https://www.htslib.org/workflow/filter.html
# ####Filter1 - based on quality, a score of >20 means that there is a 99% probability there is a variant at that site 
# bcftools filter -O v -o filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10 & INFO/DP <= 320 & MQBZ >=-5 && MQBZ <= 5 && MQSBZ >=-5 && MQSBZ <= 5 & RPBZ >= -7.5 && RPBZ <= 7.5 & SCBZ > 4' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf

