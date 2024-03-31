#-->Calculate polygenic resilience score 
cd /media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score
export PATH=~/Documents:$PATH
export PATH=~/Documents/plink:$PATH

path="/media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_GWAS"

ls /media/sf_ad_resilience_extend/HRC/post_imputation/pheno_for_update/*_phe.bim | awk '/ADNIGO23|ANM|ADC7/ {print $0}' > bim_files.txt
bim_files=`cat bim_files.txt | sed 's/.bim//g'`

for it in ${bim_files}; do
iter=`basename $it`
echo "working on" $iter


#keep samples and update pheno
#plink --bfile ${it} \
#--allow-no-sex \
#--keep ${path}/keep_resilient_replication.txt \
#--make-bed --out ${iter}_keep_rep_final

#plink --bfile ${iter}_keep_rep_final \
#--allow-no-sex \
#--pheno  ${path}/update_resilient_phe_replication.txt \
#--make-bed --out ${iter}_pheno_rep_final


plink --bfile ${iter}_pheno_rep_final \
--maf 0.05 --geno 0.02 --mind 0.02 --allow-no-sex \
--score for_AD_resilience_score_formula_design1_age60.txt  1 2 3   include-cnt sum \
--q-score-range range_file_pres.txt  for_AD_resilience_score_formula_design1_age60.txt 1 4 \
--out ${iter}_all_rep_final_age60

done




#-->Test the association between polygenic resilience scores and resilient phenotype
#install.packages("fmsb")
require(fmsb)
require(broom)

mega_coefs = list()
meta_stats = list()
pooled_stats = list()

start = c(5e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-01, 7e-01, 8e-01,9e-01, 1.0)

end = rep(0, length(start))

range_file = data.frame(names = paste("S",1:length(start), sep=""), end, start)


for(x in 1:nrow(range_file)){

message("\r Running analysis of risk scores: ", range_file$start[x])

profile=paste(".S",x,".profile",sep="")  
pval = range_file[x,3]
res_scores = list.files("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score", full.names=T, pattern = profile)
res_scores = res_scores[grepl("all_rep_final_age60", res_scores)]
res_scores = res_scores[!grepl("all_rep_final_reduced", res_scores)]
print(length(res_scores))

n_files = length(res_scores)
studyID = c("ADNIGO23","ADC7","ANM")

#Read in  scores per study, and standardize scores to mean of 0 and sd of 1
dat_list = list()
coefs_study = list()

for(y in 1:length(studyID)){

res_id = res_scores[grepl(studyID[y],res_scores)]
res_read = fread(res_id[1])

res_read = res_read[,list(SCORE = sum(SCORESUM), Counts = sum(CNT2)), by=c("FID","IID","PHENO")]
#Normalization
res_read$PRS = (res_read$SCORE - mean(res_read$SCORE)) / (sd(res_read$SCORE))
res_read = data.frame(res_read)
res_read$CODE = ifelse(res_read$PHENO %in% "2", 1,0)

res_read$IID = as.character(res_read$IID)
res_read$StudyID = studyID[y]

#Load PCs
pca_df = fread(file=paste("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/",studyID[y],"_AD_NIAGADS_association_covar_g2.txt", sep =""))
pca_df$IID = as.character(pca_df$IID)

dat = merge(pca_df, res_read[,c("IID","PRS","CODE")], by.x = "IID", by.y = "IID")
dat = dat[!is.na(dat$aaoaae),]
dat = dat[!is.na(dat$sex),]

#-->Main model
dat1 = dat[,!colnames(dat) %in% c("FID","IID","CODE"), with = F]
fit = glm(dat$CODE ~ . , data=dat1 ,family=binomial())

#-->Reduced model
dat2 = dat[,!colnames(dat) %in% c("FID","IID","PRS","CODE"), with = F]
nullFit = glm(dat$CODE ~ . , data=dat2 , family=binomial())

#-->Extract variance explained
Rsq = fmsb::NagelkerkeR2(fit)$R2 - fmsb::NagelkerkeR2(nullFit)$R2

coefs = data.frame(coef(summary(fit)))
coefs$Rsq = Rsq
coefs$N_control = nrow(dat[dat$CODE == 0, ])
coefs$N_resilient = nrow(dat[dat$CODE == 1, ])
coefs$predictor = rownames(coefs)

coefs$profile = range_file$start[x]
coefs = coefs[coefs$predictor %in% "PRS", ]
coefs$StudyID = unique(res_read$StudyID)

coefs_study[[y]] = coefs
dat_list[[y]] = res_read

}

#Combine sum stats per study
pooled_coe = ldply(coefs_study)


#Combine individual-level scores and covariates across studies
dat_res = ldply(dat_list)
pca_df = fread(file = "I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/20201130_PCA_separateByStudy_ADNIGO23_ANM_ADC7.txt",sep="\t")
datMerge = merge(dat_res,pca_df, by = "IID")

covar = fread("H:/ad_resilience/niagads/pheno/pheno_of_studies/for_filter/All_local_APOE_covar.txt")
datMerge = merge(datMerge,covar, by = "IID")

datMerge = datMerge[!is.na(datMerge$aaoaae2),]
datMerge = datMerge[!datMerge$aaoaae2 == "-9",]
datMerge = datMerge[!is.na(datMerge$sex),]
datMerge$aaoaae2 = scale(datMerge$aaoaae2)
datMerge$PRS = scale(datMerge$PRS)

#-->Mega-analysis
fit = glm(CODE ~ PRS  + PC1 + PC2 + PC3 + PC4 + aaoaae2 + sex + StudyID, data= datMerge, family=binomial())
nullFit = glm(CODE ~ PC1 + PC2 + PC3 + PC4 + aaoaae2 + sex + StudyID, data = datMerge, family=binomial())

R2 = fmsb::NagelkerkeR2(fit)$R2 - fmsb::NagelkerkeR2(nullFit)$R2
coefs = data.frame(coef(summary(fit)))
coefs = coefs[rownames(coefs) %in% "PRS", ]
coefs$Rsq = R2
coefs$profile =  range_file$start[x]
coefs$N_control = sum(pooled_coe$N_control)
coefs$N_resilient = sum(pooled_coe$N_resilient)

mega_coefs[[x]] = coefs

#-->Run meta-analysis on study effects
require(metafor)
metafor = metafor::rma(yi = pooled_coe$Estimate, sei = pooled_coe$Std..Error,  weighted = TRUE, method="DL") #ok, im pretty sure the weights argument is not necessary
pooled_coe$weights = 1/(pooled_coe$Estimate + (pooled_coe$N_control + pooled_coe$N_resilient)*(pooled_coe$Std..Error)^2)
meta_stats[[x]] = data.frame(Beta = metafor$b,
                           SE = metafor$se,
                           pval = metafor$pval,
                           I2 = metafor$I2,
                           profile = range_file$start[x],
                           Ncontrols = sum(pooled_coe$N_control),
                           Ncases = sum(pooled_coe$N_resilient),
                           Rsq = weighted.mean(x = pooled_coe$Rsq, w = pooled_coe$weights))


colnames(pooled_coe)[colnames(pooled_coe) %in% "Estimate"] = "Beta"
colnames(pooled_coe)[colnames(pooled_coe) %in% "Std..Error"] = "SE"
colnames(pooled_coe)[colnames(pooled_coe) %in% "Pr...z.."] = "pval"
pooled_stats[[x]] = pooled_coe

}
all_stats = ldply(pooled_stats)
meta = ldply(meta_stats)
mega = ldply(mega_coefs)


colnames(mega)[colnames(mega) %in% "Estimate"] = "Beta"
colnames(mega)[grepl("Pr...z..", colnames(mega))] = "pval"
colnames(mega)[grepl("Std..Error", colnames(mega))] = "SE"



fwrite(meta, file="I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score/20211107_META_ADGC_studies_allMarker_rep_ADNIGO23_ANM_ADC7.txt",quote=F,row.names=F,sep="\t")
fwrite(mega, file="I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score/20211107_mega_ADGC_studies_allMarker_rep_ADNIGO23_ANM_ADC7.txt",quote=F,row.names=F,sep="\t")
fwrite(all_stats, file="I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/resilience_score/20211107_pooled_ADGC_studies_allMarker_rep_ADNIGO23_ANM_ADC7.txt",quote=F,row.names=F,sep="\t")



#-->Plot effect size and variance explained of each polygenic resilience score
require(ggplot2)
require(cowplot)
require(RColorBrewer)
cramp = colorRampPalette(brewer.pal(n=9, "Spectral"))(length(unique(all_stats$profile)))


gplt = ggplot(all_stats, aes(y = exp(as.numeric(Beta)), x = factor(profile))) +
geom_errorbar(aes(ymin =  exp(as.numeric(Beta) - (1.96*SE)), ymax = exp(as.numeric(Beta) + (1.96*SE)) ), color="dodgerblue3", width = 0.1) +
geom_point(pch = 18, size = 2.5, col = 'navy') +
facet_wrap(~StudyID) +
theme_grey() +
geom_text(aes(y = exp(as.numeric(Beta))), label = paste("p = ", round(all_stats$pval,4), sep=""), angle = 90, size = 3, nudge_x = 0.2) +
xlab(NULL) +
ylab("Odds Ratio") +
geom_hline(yintercept = 1.0, col = "black", lty = 2, lwd = 0.5) + 
theme(strip.background = element_rect(fill = "lightgrey"), panel.border=element_rect(size = 0.5, fill = NA), axis.text.x=element_text(size = 8 ,angle = 40, hjust = 1)) 

g1r = ggplot(all_stats, aes(x = factor(profile), fill = factor(profile), ifelse(exp(as.numeric(Beta)) < 1.0, -1*Rsq, Rsq))) +
geom_bar(stat='identity', color = "black", lwd = 0.1) +
facet_wrap(~StudyID) +
scale_fill_manual(NULL, values = cramp) +
guides(fill=FALSE) +
theme_grey() +
theme(axis.text.x=element_text(size = 8, angle = 40, hjust = 1), panel.border=element_rect(size=0.5,fill=NA)) +
xlab("P-value threshold for resilience scores") +
ylab(expression(paste("Variance explained, ",italic(R)^2)))

png(paste("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/plots/20211107_Indivi_ADGC_studies_allMarker_resilience_score_rep_ADNIGO23_ANM_ADC7.png",sep=""), res=300,units="in",height=8,width=8.5)
print(plot_grid(gplt, g1r, align='v', ncol=1))
dev.off()


g1 = ggplot(meta, aes(x = factor(profile), y = exp(as.numeric(Beta)))) +
  geom_hline(yintercept = 1.0, color="grey",lwd = 0.5, lty = 2) +
  geom_text(label = paste("p = ", round(as.numeric(meta$pval),4), sep=""),  angle = 90, nudge_x = 0.3) +
  geom_errorbar(aes(ymin = exp(as.numeric(Beta)- (1.96 * SE)) , ymax = exp(as.numeric(Beta)+ (1.96 * SE)) ), width=0.1, color="dodgerblue3") +
  geom_point(pch = 18, size=2.75, color="navy") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 8, angle = 45, hjust = 1)) +
  xlab("P-value threshold for resilience scores") +
  ylab("Odds Ratio") 

require(RColorBrewer)
cramp = colorRampPalette(brewer.pal(n=9, "Spectral"))(nrow(meta))

g2 = ggplot(meta, aes(x = factor(profile), ifelse(exp(as.numeric(Beta)) < 1.0, -1*Rsq, Rsq))) +
  geom_bar(stat='identity', color = "black", lwd = 0.1, fill = rev(cramp)) +
  theme_classic() +
  theme(axis.text.x=element_text(size = 8, angle = 45, hjust = 1)) +
  xlab("P-value threshold for resilience scores") +
  ylab(expression(paste("Variance explained, R"^2)))

require(gridExtra)

png(paste("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/plots/20211107_meta_coef_ADGC_studies_allMarker_resilience_score_rep_ADNIGO23_ANM_ADC7.png",sep=""),res=300,units="in",height=4,width=8)
grid.arrange(g1,g2,ncol=2)  
dev.off()


g1 = ggplot(mega, aes(x = factor(profile), y = exp(as.numeric(Beta)))) +
  geom_hline(yintercept = 1.0, color="grey",lwd = 0.5, lty = 2) +
  geom_text(label = paste("p = ", round(as.numeric(mega$pval),4), sep=""),  angle = 90, nudge_x = 0.3) +
  geom_errorbar(aes(ymin = exp(as.numeric(Beta)- (1.96 * SE)) , ymax = exp(as.numeric(Beta)+ (1.96 * SE)) ), width=0.1, color="dodgerblue3") +
  geom_point(pch = 18, size=2.75, color="navy") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 8, angle = 45, hjust = 1)) +
  xlab("P-value threshold for resilience scores") +
  ylab("Odds Ratio")

require(RColorBrewer)
cramp = colorRampPalette(brewer.pal(n=9, "Spectral"))(nrow(mega))

g2 = ggplot(mega, aes(x = factor(profile), ifelse(exp(as.numeric(Beta)) < 1.0, -1*Rsq, Rsq))) +
  geom_bar(stat='identity', color = "black", lwd = 0.1, fill = rev(cramp)) +
  theme_classic() +
  theme(axis.text.x=element_text(size = 8, angle = 45, hjust = 1)) +
  xlab("P-value threshold for resilience scores") +
  ylab(expression(paste("Variance explained, R"^2)))

require(gridExtra)

png(paste("I:/ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/plots/20211107_MEGA_coef_ADGC_studies_allMarker_resilience_score_rep_ADNIGO23_ANM_ADC7.png",sep=""),res=300,units="in",height=4,width=8)
grid.arrange(g1,g2,ncol=2)  
dev.off()

