#GWAS on resilience 

#-->keep resilient samples from post-imputation genotypes and recode phenotypes 
cd /media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS
export PATH=~/Documents:$PATH
export PATH=~/Documents/plink_linux_i686_20201019:$PATH


ls /media/sf_ad_resilience_extend/HRC/post_imputation/pheno_for_update/*_phe.bim | sed '/GenADA\|kramer_merged\|ADNI_GO_2\|ADNI_GO2\|ADC7/d' > bim_files.txt
study_id=(ADC1 ADC2  ADC3 ADC4 ADC5 ADC6 VMUMMSSM goate2 mirw  ROSMAP  LOAD  ACT mayo Omni ADNI_1 kramernoA TGen rosmap2 whicap WASHU2 TARCC MTC)
path=/media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3

for id in "${study_id[@]}"; do
it=$(grep "$id" bim_files.txt | sed 's/.bim//g')
iter=`basename ${it}`
echo "working on" $iter


#keep samples and update pheno
plink --bfile ${it} \
--allow-no-sex \
--keep keep_resilient_discovery.txt \
--make-bed --out ${iter}_keep_discovery

plink --bfile ${iter}_keep_discovery \
--allow-no-sex \
--pheno  update_resilient_phe_discovery.txt \
--make-bed --out ${iter}_pheno_discovery

#allele frequency
plink --bfile ${iter}_pheno_discovery \
--allow-no-sex \
--chr 1-22   \
--maf 0.01 --geno 0.02 --mind 0.02 --hwe 1e-06 \
--freq counts  \
--out ${iter}_discovery

#association
  # == Run QC:
  # variant call rate > 98%
  # subject call rate > 98%
  # HWE p < 1e-06
  # allow no sex
  # maf > 1%
  # == Perform statistical analysis:
  # logistic regression
  # covariates: top 3 PCs

plink --bfile ${iter}_pheno_discovery \
--allow-no-sex \
--chr 1-22 \
--maf 0.01 --geno 0.02 --mind 0.02  --hwe 1e-06  \
--logistic beta hide-covar \
--ci 0.95 \
--covar ${path}/${id}_AD_NIAGADS_association_covar_g2.txt \
--out ${iter}_noCOV_discovery_g2

done




#-->Combine allele frequency file with association results file
cd /media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS
assoc=`ls *_noCOV_discovery_g2.assoc.logistic | sed '/ADNI_GO_2\|ADNI_GO2\|ADC7\|a_NG\|b_NG\|c_NG\|kramer_merged/d'`

for file in $assoc; do

study_id=`basename ${file%_noCOV_discovery_g2.assoc.logistic}`
echo ${study_id}
freq=`ls ${study_id}_discovery.frq.counts`

awk 'NR==FNR { n[$2] = $2; next } ($2 in n) {print $0}' ${file}  ${freq} > freq.temp

freq_temp="/media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS/freq.temp"

paste <(awk '{ gsub(/\r/,"", $12); print } ' ${file}) ${freq_temp}  | tail -n +2 > merge.temp

echo "Removing missing data from:" ${file}

grep -v NA merge.temp > ${study_id}_discovery_g2.logistic.merge

done



#-->Create file for meta-analysis
niagads=list.files(path="I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS", pattern="_discovery_g2.logistic.merge", full.names=T)
niagads = niagads[!grepl(paste(c("ADNI_GO_2","ADNI_GO2","ADC7","a_NG","b_NG","c_NG","ADC2"),sep="", collapse = "|"), niagads)]

file_list = niagads

for(x in 1:length(file_list)){
  
  basename = gsub("_discovery_g2.logistic.merge", "", basename(file_list[[x]]))
  # read merge file
  merge_dat = NA
  merge_dat = fread(file_list[[x]], sep = " ", fill = F, col.names = c('CHR_asso','SNP_asso','BP','A1_asso','TEST','NMISS','BETA','SE','L95','U95','STAT','P','CHR','SNP','a1','a2','C1','C2','G0'))
  
  merge_dat$check = ifelse(merge_dat$A1_asso == merge_dat$a1, "Y","N")
  print(table(merge_dat$check))
  merge_dat$A1 = ifelse(merge_dat$check == "Y", merge_dat$a1,merge_dat$a2)
  merge_dat$A2 = ifelse(merge_dat$check == "Y", merge_dat$a2,merge_dat$a1)
  merge_dat$A1_counts = ifelse(merge_dat$check == "Y", merge_dat$C1,merge_dat$C2)
  merge_dat$A2_counts = ifelse(merge_dat$check == "Y", merge_dat$C2,merge_dat$C1)
  merge_dat$total_counts = merge_dat$A1_counts + merge_dat$A2_counts
  merge_dat$study_label = "1"

  merge_dat$check2 = ifelse(merge_dat$A1_asso == merge_dat$A1, "Y","N")
  print(table(merge_dat$check2)) 
  
  fwrite(merge_dat,
         file=paste("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS/",basename,"_discovery_g2.formeta.assoc.txt",sep=""),
         quote=F, row.names=F, col.names = TRUE, sep= "\t")
}

#

#-->QQ plot of each study
niagads=list.files(path="I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS", pattern="_discovery_g2.formeta.assoc", full.names=T)
niagads = niagads[!grepl(paste(c("ADNI_GO_2","ADNI_GO2","ADC7","a_NG","b_NG","c_NG","ADC2"),sep="",collapse="|"), niagads)]

file_list = niagads

recordedplot = list()
for(x in 1:length(file_list)){
name = gsub("_discovery_g2.formeta.assoc.txt","",basename(file_list[[x]]))
name = gsub("POSTIMP_","",name)
name = gsub("_merge_d_phe","",name)

sum = fread(file_list[[x]])

observed <- sort(sum$P)
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

chisq <- qchisq(sum$P,1,lower.tail=F)
lambda = median(chisq)/qchisq(0.5,1) #1.750368

plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
text(x=2,y=6,paste("lambda =", format(lambda,digits=5),sep=""))
title(main = paste(name,"_QQ-plot", sep = ""), cex.main = 1)
recordedplot[[x]] <- recordPlot()
dev.off()

}

library(cowplot)
require(gridGraphics)
png("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/plots/20201126_eachStudy_qqplot_rerun_g2_part1.png", width=1200, height=1200)
plot_grid(plotlist=recordedplot[1:10], hjust = 0, vjust = 1 )
dev.off()
png("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/plots/20201126_eachStudy_qqplot_rerun_g2_part2.png", width=1200, height=1200)
plot_grid(plotlist=recordedplot[11:21], hjust = 0, vjust = 1 )
dev.off()



#-->Create a meta-analysis script for METAL
study_list1 = list.files(path = "I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS", pattern = "_discovery_g2.formeta.assoc.txt", full.names=FALSE)
study_list1 = study_list1[!grepl(paste(c("ADNI_GO_2","ADNI_GO2","ADC7","a_NG","b_NG","c_NG","ADC2"),sep="", collapse = "|"), study_list1)]

study_list = study_list1
#length(study_list)

header = paste("
SEPARATOR WHITESPACE
COLUMNCOUNTING LENIENT
MARKERLABEL SNP
ALLELELABELS A1 A2
PVALUELABEL P
EFFECT BETA
STDERR SE
SCHEME STDERR
GENOMICCONTROL ON
CUSTOMVARIABLE TOTAL_ALLELE
LABEL TOTAL_ALLELE as total_counts
CUSTOMVARIABLE A1_Tcounts
LABEL A1_Tcounts as A1_counts
CUSTOMVARIABLE A2_Tcounts
LABEL A2_Tcounts as A2_counts
CUSTOMVARIABLE study_Tcounts
LABEL study_Tcounts as study_label
OUTFILE METAANALYSIS_DISCOVERY_meta_g2_ .TBL")


proc_study = list()
for(x in 1:length(study_list)){
  proc_study[[x]] = paste("PROCESS ", study_list[[x]], sep= "")
}

full_script = c(paste(header), paste(proc_study), paste("\nANALYZE", sep=""))

write.table(full_script,
            file="I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS/metal_script_discovery_meta_g2.txt", quote=F,row.names=F,col.names=F)



#-->Run meta-analysis
export PATH=~/Documents/generic-metal/executables:$PATH
cd /media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS
metal metal_script_discovery_meta_g2.txt




