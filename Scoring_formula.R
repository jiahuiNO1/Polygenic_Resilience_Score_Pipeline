#Derive polyugenic scoring formula from GWAS meta-analysis summary statistics
#Only keep risk-residual genetic variants

#-->Summary statistics QC
sum = fread("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS/20211105_METAANALYSIS_DISCOVERY_meta_final_full_1.TBL", header = T)
colnames(sum)[[6]] = "Pvalue"
sum[sum$Allele1 == "a",]$Allele1 = "A"
sum[sum$Allele1 == "c",]$Allele1 = "C"
sum[sum$Allele1 == "t",]$Allele1 = "T"
sum[sum$Allele1 == "g",]$Allele1 = "G"

sum[sum$Allele2 == "a",]$Allele2 = "A"
sum[sum$Allele2 == "c",]$Allele2 = "C"
sum[sum$Allele2 == "t",]$Allele2 = "T"
sum[sum$Allele2 == "g",]$Allele2 = "G"

sum$OR = exp(sum$Effect)
nrow(sum) # 8321202

#-->Remove bad SNPs
badsnp = fread("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_sce1/risk/RemovetheseVariants_RiskVar_LD.txt")
sum = sum[!sum$MarkerName %in% badsnp$SNP,]
nrow(sum) #336225

#-->Remove high-LD region (MHC: chr6: 25:34 mb) (chr8: 7MB - 14MB)
out = fread("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/HRC_reference/MHC_CHR8_outlist.txt", header = T)

sum = sum[!sum$MarkerName %in% out$SNP,]
#nrow(sum) 
prs_stats = sum

#-->Remove SNPs with MAF < 0.05
prs_stats$A1_fre = prs_stats$A1_Tcounts/prs_stats$TOTAL_ALLELE
prs_stats$A2_fre = prs_stats$A2_Tcounts/prs_stats$TOTAL_ALLELE
prs_stats$minor_fre = ifelse(prs_stats$A1_fre < prs_stats$A2_fre, prs_stats$A1_fre, prs_stats$A2_fre)

outlist = prs_stats[prs_stats$minor_fre < 0.05, ]$MarkerName
prs_stats = prs_stats[!prs_stats$MarkerName %in% outlist, ]
#nrow(prs_stats)  

#-->Remove ambiguous SNPs
test1= prs_stats[prs_stats$Allele1 == "A" & prs_stats$Allele2 == "T",]
test2= prs_stats[prs_stats$Allele1 == "T" & prs_stats$Allele2 == "A",]
test3= prs_stats[prs_stats$Allele1 == "G" & prs_stats$Allele2 == "C",]
test4= prs_stats[prs_stats$Allele1 == "C" & prs_stats$Allele2 == "G",]
snp1 = test1$MarkerName
snp2 = test2$MarkerName
snp3 = test3$MarkerName
snp4 = test4$MarkerName
snplist = c(snp1, snp2, snp3, snp4)

prs_stats = prs_stats[!prs_stats$MarkerName %in% snplist, ]

#-->Remove duplicated SNPs
table(duplicated(prs_stats$MarkerName))
#FALSE 
#81194 

#-->Remove in-del SNPs
mullist1 = prs_stats[ nchar(prs_stats$Allele1) > 1,]$MarkerName
mullist2 = prs_stats[ nchar(prs_stats$Allele2) > 1,]$MarkerName
mullist = unique(c(mullist1,mullist2))
prs_stats = prs_stats[!prs_stats$MarkerName %in% mullist, ]
prs_stats$SE = prs_stats$StdErr
head(prs_stats)
nrow(prs_stats) # 81194

#-->Remove multiallelic variants
#all multiallelic SNPs generated in 20201114 scenario3 V3 scripts
rmlist = fread("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/HRC_reference/all_multiallelic_SNPs.txt", header = T)
table(prs_stats$MarkerName %in% rmlist$delete) # TRUE 8942
prs_stats = prs_stats[!prs_stats$MarkerName %in% rmlist$delete, ]
nrow(prs_stats) #72252

fwrite(prs_stats,
       file="I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score/20211105_resilience_Summary_statistics_with_allelefreqs_0.05_QCed_final_full.txt",
       quote=F,row.names=F,sep="\t")


#-->Clump summary markers with a MAF > 5% associated with resilience using Plink and reference panel 
cd /media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score
export PATH=~/Documents:$PATH
export PATH=~/Documents/plink:$PATH

dir=/media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/HRC_reference
plink --bfile ${dir}/ADNI_GO2_HRC_reference_noSAme \
      --clump 20211105_resilience_Summary_statistics_with_allelefreqs_0.05_QCed_final_full.txt \
      --clump-field Pvalue \
      --clump-snp-field MarkerName \
      --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.2 --clump-kb 250 \
      --out 20211105_AD_resilience_summary_clump_final_full
#-clump: 28808 clumps formed from 65779 top variants.


clumped = list.files("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score", full.names=T,pattern=".clumped")
clumped = clumped[grepl("20211105_", clumped)]
length(clumped)
read_clumped = lapply(clumped, function(x) fread(x,h=T,fill=T))
read_clumped = rbindlist(read_clumped)
read_clumped = data.frame(read_clumped)

prs_stats = fread("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score/20211105_resilience_Summary_statistics_with_allelefreqs_0.05_QCed_final_full.txt")
prs_stats = prs_stats[prs_stats$MarkerName %in% read_clumped$SNP, ]
nrow(prs_stats)
#28808

#-->Recode weights
prs_stats$WEIGHT = ifelse(prs_stats$Effect > 0, prs_stats$Effect, prs_stats$Effect * -1)
prs_stats$Resilience_allele = ifelse(prs_stats$Effect > 0, prs_stats$Allele1, prs_stats$Allele2)


#-->SNPs with study counts less than 5 were removed
prs_stats = prs_stats[!prs_stats$study_Tcounts < 5,]
nrow(prs_stats)
# 18723

hrc = fread("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/HRC_reference/ADNI_GO2_HRC_reference_noSAme.bim")

mer = merge(prs_stats,hrc, by.x = "MarkerName", by.y = "V2")
nrow(mer)

prs_stats = mer
prs_stats$Chromosome = prs_stats$V1
prs_stats$Position = prs_stats$V4

prs_stats = prs_stats[,colnames(prs_stats) %in% c("MarkerName","Resilience_allele","WEIGHT","Pvalue","Allele1","Allele2", "Chromosome", "Position"), with = F]

prs_stats = prs_stats[, c("MarkerName","Resilience_allele","WEIGHT","Pvalue","Allele1","Allele2", "Chromosome", "Position")]

fwrite(prs_stats,
       file="I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score/for_AD_resilience_score_formula_design1_age60.txt",
       quote=F,row.names=F, sep="\t")




