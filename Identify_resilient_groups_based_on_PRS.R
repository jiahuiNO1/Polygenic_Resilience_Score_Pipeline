#-->Calculate PRS (polygenic risk score)
ADGC = list.files("I:/ad_resilience_extend/HRC/post_imputation/pheno_for_update/", pattern ="_phe.bim", full.names=TRUE)

for(x in 1:length(ADGC)){
  cat("\rWorking on file:", x)
  BFILE = gsub(".bim", "", ADGC[[x]])
  cmd = paste("I:/ad_resilience_extend/HRC/learning_curve/plink_win64_20191130/plink --bfile ",BFILE, 
  	" --maf 0.05 --geno 0.02 --mind 0.02 --allow-no-sex --score I:/ad_resilience_extend/HRC/CHARGE_analysis/CHARGE_Script_V2/for_AD_PRS_formula_CHARGE_sce3.txt 3 6 5 include-cnt sum --q-score-range I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run/range_file_pres.txt I:/ad_resilience_extend/HRC/CHARGE_analysis/CHARGE_Script_V2/for_AD_PRS_formula_CHARGE_sce3.txt 3 4 --out I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/",basename(BFILE), sep="")
  system(cmd, show.output.on.console = FALSE)  
}


#-->Load risk scores and examine risk distribution in cases and controls 
library(viridis)
library(gridExtra)
require(RColorBrewer)
start = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3,  1e-2, 5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-01, 7e-01, 8e-01,9e-01, 1.0)
end = rep(0, length(start))

range_file = data.frame(names = paste("S",1:length(start), sep=""), end, start)

for(x in 1:nrow(range_file)){

#load PRSs
profile=paste(".S",x,".profile",sep="")  
pval = range_file[x,3]
risk_scores = list.files("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/", full.names=T, pattern = profile)
prefix= c("ADC1" ,"ADC2", "ADC3" ,"ADC4", "ADC5" ,"ADC6", "VMUMMSSM" ,"mayo", "kramernoA", "goate2", "ADNI_1", "ACT" ,"LOAD", "ROSMAP", "mirw" ,"Omni", "TGen" ,"MTC" ,"WASHU2", "whicap", "rosmap2" ,"TARCC" )

risk_scores = risk_scores[grepl(paste(prefix, sep="", collapse="|"), risk_scores)]
print(length(risk_scores))
read_scores = lapply(risk_scores, function(x) fread(x,header=TRUE))
names(read_scores) = gsub(profile, "", basename(risk_scores))
df = ldply(read_scores)
colnames(df)[1] = "StudyID"

#load PCs
pca = fread("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/20201122_PCA_separateByStudy_AD-NIAGADS.txt")
df = merge(df, pca, by = "IID")

#load APOE genotype
apoe = fread("H:/ad_resilience/niagads/pheno/pheno_of_studies/for_filter/All_local_APOE_covar.txt")
apoe$IID = as.character(apoe$IID)
df = merge(df, apoe, by = "IID")

#remove individuals with age < 60
df$affection[df$PHENO == 1] = "Control"
df$affection[df$PHENO == 2] = "AD"
df = df[!is.na(df$affection),]
df = df[!df$aaoaae2 < 60,]
df = df[!is.na(df$sex),]
df = df[!df$sex == "-9",]

#scaling PRSs
study_id = unique(df$StudyID)
tmp_list=list()
for(y in 1:length(study_id)){
  tmp = df[df$StudyID %in% study_id[[y]], ]
  tmp$PRS = scale(tmp$SCORE)
  tmp_list[[y]] = tmp
}
datAll = ldply(tmp_list)

#find max AD-PRS and 90th percentile within controls per study, and identify high-risk resilient groups
study_id = unique(datAll$StudyID)
tmp_list=list()
for(y in 1:length(study_id)){
  tmp = datAll[datAll$StudyID %in% study_id[[y]], ]
  upper_bn = max(tmp$PRS[tmp$PHENO == 1])
  lower_bn = quantile(tmp$PRS[tmp$PHENO == 1], 0.9)[[1]]
  tmp$upper_bn = upper_bn
  tmp$lower_bn = lower_bn
  tmp$HR_STATUS[tmp$PHENO == 1 & tmp$PRS >= lower_bn] = "High-risk-CN"
  tmp$HR_STATUS[tmp$PHENO == 2 & tmp$PRS >= lower_bn & tmp$PRS <= upper_bn] = "High-risk-AD"
  
  tmp_list[[y]] = tmp
}
datAll = ldply(tmp_list)


write.table(datAll, file=paste("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/AD_PRS-Scores_per_study_P",pval,".txt", sep=""),quote=F,row.names=F)

hr_status_n = table(datAll$HR_STATUS)

title = paste(names(hr_status_n)[[1]], " = ", hr_status_n[[1]], ", ", names(hr_status_n)[[2]], " = ", hr_status_n[[2]], sep= "")

#plot PRSs distributions in cases and controls
g = ggplot(datAll, aes(x = PRS, fill = affection)) +
  geom_density(col = 'black', lwd = 0.3, alpha = 0.6) +
  theme_classic() +
  geom_vline(data = datAll, aes(xintercept = lower_bn), color="red", lwd =0.75 , lty = 2) +
  geom_vline(data = datAll, aes(xintercept = upper_bn), color="red", lwd =0.75, lty = 2) +
  scale_fill_manual('Disease status', values = c("purple","gold")) +
  facet_wrap(~StudyID) +
  ggtitle(title) +
  ylab("Density distribution of subjects") +
  xlab("AD-PRS (z-standardized)") +
  theme(axis.text=element_text(size = 12, color="black"), axis.title=element_text(size = 12), panel.border=element_rect(size = 1, fill = NA))


png(paste("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/plots/PRS_resilient_sample_indiStudy_distribution_P",pval,"_histogram.png",sep=""),res=300,units="in",height=6,width=25)
print(g)
dev.off()

}

