####Methods used to filter Eragrostis SNPs

####SNP call file 
snpcalls=/home/jconnell/niab/eragrostis/snpcalls

####Convert bcf to vcf 
bcftools convert \
-Ov \
-o ${snpcalls}/teff.Dabbi_50954_V3_kora_asgori.vcf \
${snpcalls}/teff.Dabbi_50954_V3_kora_asgori.bcf

####Sep INDELs from SNPs 
cat teff.Dabbi_50954_V3_kora_asgori.vcf | grep "#\|INDEL" > INDEL_teff.Dabbi_50954_V3_kora_asgori.vcf
cat	teff.Dabbi_50954_V3_kora_asgori.vcf | grep -v "INDEL" > SNP_teff.Dabbi_50954_V3_kora_asgori.vcf

####Create bed file for INDEL SNP filter and depth/qual data for plotting 
python /home/jconnell/git_repos/niab_repos/eragrostis/extract_vcf_info.py -snp SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -indel INDEL_teff.Dabbi_50954_V3_kora_asgori.vcf

####Exclude SNPs round INDELs
vcftools --vcf SNP_teff.Dabbi_50954_V3_kora_asgori.vcf --exclude-bed snp_excluson_data.txt --keep-INFO-all --recode --stdout > filtered_indel_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf

####Split up VCF by sample 
for x in $(echo $(bcftools query -l filtered_indel_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf)); do
	bcftools view -Ov -s $x -o "$x"_SNP_indelfiltered.vcf filtered_indel_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf
done

####Apply depth and qual filters to individual files 
#vcftools --vcf Asgori-Tank1_SNP_indelfiltered.vcf --minDP 25 --maxDP 140 --keep-INFO-all --recode --stdout > Asgori-Tank1_SNP_indelfiltered_depthfilter.vcf

bcftools filter -O v -o Asgori-Tank1_SNP_indelfiltered_depth_qualfilter.vcf -i 'QUAL >= 20 & INFO/DP <= 140 && INFO/DP >= 25' Asgori-Tank1_SNP_indelfiltered.vcf
bcftools filter -O v -o Kora-Tank1_SNP_indelfiltered_depth_qualfilter.vcf -i 'QUAL >= 20 & INFO/DP <= 150 && INFO/DP >= 40' Kora-Tank1_SNP_indelfiltered.vcf

# #For fitlering bcftools uses "&" to keep vairiants that meet all given criteria and "&&" to keep variants that match at least one of the given criteria
# https://www.htslib.org/workflow/filter.html
# ####Filter1 - based on quality, a score of >20 means that there is a 99% probability there is a variant at that site 
# bcftools filter -O v -o filter_SNP_teff.Dabbi_50954_V3_kora_asgori.vcf -i 'QUAL >= 10 & INFO/DP <= 320 & MQBZ >=-5 && MQBZ <= 5 && MQSBZ >=-5 && MQSBZ <= 5 & RPBZ >= -7.5 && RPBZ <= 7.5 & SCBZ > 4' SNP_teff.Dabbi_50954_V3_kora_asgori.vcf

