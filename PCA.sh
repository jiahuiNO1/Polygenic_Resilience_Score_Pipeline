#To get PCs to account for population structure of genotypes

#
cd /media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3
export PATH=~/Documents:$PATH
export PATH=~/Documents/EIG-7.2.1/bin:$PATH
export PATH=~/Documents/plink_linux_i686_20201019:$PATH


ls /media/sf_ad_resilience_extend/preimp/*_initial_pca.bim | sed '/GenADA/d' > bim_files1.txt
bim_files=`cat bim_files1.txt | sed 's/.bim//g'`

new_pfx="FO"
prefix=( "Omni" )

path="/media/sf_ad_resilience_extend/HRC/whole_resilience_pipeline_1st_run_sce3/PCA_PHENO_UPDATE"

for it in "${prefix[@]}"; do
echo "working on" $it
iter=$(echo "${bim_files}" | grep "$it")
echo ${iter}

#-->Update sample ID
plink --bfile ${iter} \
--update-ids ${path}/update_ID_for_Omni.txt \
--make-bed --out ${it}_names

#-->Keep the samples in the phenotype file
plink --bfile ${it}_names \
--allow-no-sex \
--keep ${path}/${it}_keep.txt \
--make-bed --out ${it}_keep

#-->Update phenotype
plink --bfile ${it}_keep \
--allow-no-sex \
--pheno  ${path}/${it}_pheno.txt \
--make-bed --out ${it}_pheno

#-->LD prune
plink --bfile ${it}_pheno \
--filter-founders \
--maf 0.05 \
--indep-pairwise 100 50 0.1 \
--out pruned

#-->extract independent variants and samples and then PCA
plink --bfile ${it}_pheno \
--extract pruned.prune.in \
--make-bed --out ${it}_pruned

awk -f highLDregion_b37.awk ${it}_pruned.bim > highLDexcludes
awk '($1 < 1) || ($1 > 22) {print $2}' ${it}_pruned.bim > autosomeexcludes
cat highLDexcludes autosomeexcludes > highLD_and_autosomal_excludes  

plink --bfile ${it}_pruned \
--exclude highLD_and_autosomal_excludes \
--make-bed --out ${it}_pruned_two

#-->Run EIGENSOFT using LD-pruned binary
#Convert files to EIGENSOFT format using CONVERTF
#Requires par file to convert from packedped format to eigenstrat format
convertf -p <(printf "genotypename: "${it}"_pruned_two.bed
snpname: "${it}"_pruned_two.bim
indivname: "${it}"_pruned_two.fam
outputformat: EIGENSTRAT
genotypeoutname: "${it}".pop_strat.eigenstratgeno
snpoutname: "${it}".pop_strat.snp
indivoutname: "${it}".pop_strat.ind
pordercheck: NO
outlieroutname: "${it}"_out.outliers")


#-->Run SmartPCA again to remove outliers
smartpca.perl -i ${it}.pop_strat.eigenstratgeno \
-a ${it}.pop_strat.snp \
-b ${it}.pop_strat.ind \
-o ${it}.pop_strat_outliers.pca \
-p ${it}.pop_strat_outliers.plot \
-e ${it}.pop_strat_outliers.eval \
-l ${it}.pop_strat_outliers_smartpca.log \
-m 20 \
-t 1 \
-k 20 \
-s 6


#-->Extract outliers
awk '/REMOVED/ {print $3}' ${it}.pop_strat_outliers_smartpca.log | sed 's/:/ /g' > ${it}.pop_strat_outliers.outliers

plink --bfile ${it}_pruned_two \
--remove ${it}.pop_strat_outliers.outliers \
--make-bed --out ${it}_pruned_removed


#-->Re-run to assess which components to include as covariates in the final analysis
#Run ConvertF:
convertf -p <(printf "genotypename: "${it}"_pruned_removed.bed
snpname: "${it}"_pruned_removed.bim
indivname: "${it}"_pruned_removed.fam
outputformat: EIGENSTRAT
genotypeoutname: "${it}".PCS_for_covariates.eigenstratgeno
snpoutname: "${it}".PCS_for_covariates.snp
indivoutname: "${it}".PCS_for_covariates.ind
pordercheck: NO")

#Run SmartPCA:
smartpca.perl -i ${it}.PCS_for_covariates.eigenstratgeno \
-a ${it}.PCS_for_covariates.snp \
-b ${it}.PCS_for_covariates.ind \
-o ${it}.PCS_for_covariates.pca \
-p ${it}.PCS_for_covariates.plot \
-e ${it}.PCS_for_covariates.eval \
-l ${it}.PCS_for_covariates_smartpca.log \
-m 0 \
-t 20 \
-k 20 \
-s 6 

#Calculate association (short version):
sed -i -e 's/^[ \t]*//' -e 's/:/ /g' ${it}.PCS_for_covariates.pca.evec

R --file=PC-VS-OUTCOME_IN_R_SHORT.R --args ${it}.PCS_for_covariates ${path}/${it}_pheno.txt

done
