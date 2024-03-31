#!/user/bin/bash
cd /media/sf_ad_resilience_extend/preimp
export PATH=~/Documents:$PATH
export PATH=~/Documents/plink_linux_i686_20201019:$PATH

new_pfx="FO"
bim_files=`ls /media/sf_ad_resilience/niagads/genotype/out_for_HRC_imputation/*.bim`
path='/media/sf_ad_resilience/niagads/genotype'


for iter in ${bim_files}; do

echo $iter
it=`basename ${iter%%.*}`
echo ${it}

echo "working on:" $path/${it}  

#-->impute sex
#here is just for getting the sex check report
#plink --bfile $path/${it} --impute-sex --make-bed --out sex_check
#awk '{ if ($5 == "PROBLEM") {print $1, $2} }' sex_check.sexcheck > badsex.txt
#plink --bfile $path/${it} --remove badsex.txt --make-bed --out ${it}.sex 


#-->find and exclude duplicate variants
plink --bfile ${path}/${it} --list-duplicate-vars ids-only suppress-first --out dupvars 
plink --bfile ${path}/${it} --exclude dupvars.dupvar --make-bed --out ${new_pfx}_${it}_NOdupvars


#-->find and exclude missing variants and samples with too much missing variants calls
plink --bfile ${new_pfx}_${it}_NOdupvars --set-missing-var-ids @:#_$1_$2 --make-bed --out ${new_pfx}_${it}_setmissingvar
plink --bfile ${new_pfx}_${it}_setmissingvar --geno 0.02 --maf 0.01 -hwe 1e-06 --make-bed --out ${new_pfx}_${it}_preimpqc
plink --bfile ${new_pfx}_${it}_preimpqc --mind 0.02 --make-bed --out ${new_pfx}_${it}_preimpqc_mind


#-->limit to biallelic
plink --bfile ${new_pfx}_${it}_preimpqc_mind --biallelic-only strict --make-bed --out ${new_pfx}_${it}_preimpqc_mind_bi


rm FORCE_A1.txt
rm FLIP_STRAND.txt

#-->match allele with 1000 genome reference
echo '
# load pkg
require(data.table)
# import TARGET bim
target_bim = fread("'${new_pfx}_${it}_preimpqc_mind_bi'.bim")
ref_bim = fread("/media/sf_ad_resilience/niagads/reference/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01.bim")

# v2: variants(snp) name v4 : SNP position
ref_bim = ref_bim[match(target_bim$V2, ref_bim$V2)]
target_bim$V4 = ifelse(!is.na(ref_bim$V2), ref_bim$V4, target_bim$V4)

# Overwrite BIM file with correct SNP positions 
fwrite(target_bim, file = "'${new_pfx}_${it}_preimpqc_mind_bi'.bim",quote=F,row.names=F,sep="\t",col.names=FALSE)

# check that SNPs in target overlap reference
overlaps = intersect(target_bim$V2, ref_bim$V2)
if(length(overlaps) < 100){stop("Warning! <= 100 SNPs are common to the target and reference panel. Check that rsIDs are correct in the target.")}
target_bim = target_bim[target_bim$V2 %in% overlaps]
ref_bim = ref_bim[ref_bim$V2 %in% overlaps]
target_bim = target_bim[match(ref_bim$V2, target_bim$V2)] 

# check that REF/ALT alleles match exactly
a1_match = target_bim$V5 == ref_bim$V5
a2_match = target_bim$V6 == ref_bim$V6

# find potential 
missnp = target_bim$V2[which(a1_match == FALSE)]

a1_missnp = target_bim$V5[target_bim$V2 %in% missnp] == ref_bim$V6[ref_bim$V2 %in% missnp]
a2_missnp = target_bim$V6[target_bim$V2 %in% missnp] == ref_bim$V5[ref_bim$V2 %in% missnp]

#switch allele
switch_allele = missnp[which(a1_missnp == TRUE)]

if(length(switch_allele) > 1){
force_a1 = ref_bim[ref_bim$V2 %in% switch_allele, colnames(ref_bim) %in% c("V2","V5"),with=FALSE]
fwrite(force_a1, file="FORCE_A1.txt", quote=F,sep="\t",col.names=FALSE,row.names=FALSE)
}'  | R --vanilla



#-->update reference allele
FILE='FORCE_A1.txt'
if [ -f  $FILE ];
then
	plink --bfile  ${new_pfx}_${it}_preimpqc_mind_bi --reference-allele FORCE_A1.txt --make-bed --out ${new_pfx}_${it}_preimpqc_mind_bi_align
else
	plink --bfile  ${new_pfx}_${it}_preimpqc_mind_bi --make-bed --out ${new_pfx}_${it}_preimpqc_mind_bi_align
fi

plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align --list-duplicate-vars ids-only  --out dupvars2
plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align  --exclude dupvars2.dupvar --make-bed --out ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup


#-->relatedness check
#prune: include only founders; LD prune: window size(kb), variants counts for shifting window, r^2 threshold
plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup \
--filter-founders \
--maf 0.05 \
--indep-pairwise 100 50 0.1 \
--out prunedlist 

rm IBD_relcheck.genome

#mhc region (create for IBD computation)
plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup \
--chr 6 --from-mb 24 --to-mb 35 \
--write-snplist --out mhc

#include only independent samples, no mhc, only autosome, genome: IBD computation
plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup \
--extract prunedlist.prune.in \
--exclude mhc.snplist \
--chr 1-22 \
--genome \
--out IBD_relcheck

rm problem_rel.txt
awk '{ if ($10 >= 0.2) { print $1, $2} }' IBD_relcheck.genome > problem_rel.txt

FILE="problem_rel.txt"
if [ -f $FILE ]
then 
	plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup --remove problem_rel.txt --make-bed --out ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial
else
	plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup --make-bed --out ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial
fi


echo "finishing initial QC" # sex check, duplicate, missing, biallelic, relatedness

done


#-->PCA
bim_files=`ls /media/sf_ad_resilience/niagads/genotype/out_for_HRC_imputation/*.bim`

for iter in ${bim_files}; do
echo $iter
it=`basename ${iter%%.*}`

echo "working on:" ${it} 
#mhc region (create for prune)
plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial \
--chr 6 --from-mb 24 --to-mb 35 \
--write-snplist --out mhc
# chr 8 exclude(create for prune)
plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial \
--chr 8 --from-mb 8 --to-mb 14 \
--write-snplist --out chr8inv
#combine snplist
cat mhc.snplist chr8inv.snplist > excludesnplist.txt

echo "performing PCA for outliers detection"

#population outliers in PCA space
#LD prune
plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial \
--filter-founders \
--exclude excludesnplist.txt \
--maf 0.05 \
--indep-pairwise 100 50 0.1 \
--out pruned
#extract independent variants and samples and then PCA
plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial \
--extract pruned.prune.in \
--filter-founders \
--pca --out pca_founders
#edit pca header
echo '
require(data.table)
eigen = list.files( pattern=".eigenvec", full.names=T)
eigen = eigen[grepl("founders", eigen)]
for(x in 1:length(eigen)){
header = readLines(eigen[[x]],n=1)
headerCheck = header[grepl("PC", header)]
if(length(headerCheck) > 0) next
read = fread(eigen[[x]], h=F)
colnames(read)[1:2] = c("FID","IID")
colnames(read)[3:ncol(read)] = paste("PC", 1:(ncol(read)-2), sep= "")
fileout = gsub(".eigenvec", ".pca", eigen[[x]])
fwrite(read, file=fileout, quote=F,row.names=F,sep="\t")
} ' | R --vanilla

cat pca_founders.pca > ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial.pca 

rm OUTLIER*txt

echo '
require(data.table)
pca_files = list.files(pattern="founders.pca")
pca_files = pca_files[grepl("founders", pca_files)]
for(x in 1:length(pca_files)){
read_in = fread(pca_files[[x]],h=T)
read_in = as.data.frame(read_in)
foundOut = list()
for(y in 1:10){
pca_scores = read_in[,!colnames(read_in) %in% c("FID","IID")]
rownames(pca_scores) = read_in$IID
mean_pc = mean(pca_scores[,y])
var_pc = sd(pca_scores[,y])
thres = var_pc*6
foundOut[[y]] = rownames(pca_scores)[which(pca_scores[,y] > mean_pc + thres)]
message("\rDetected ", length(foundOut[[y]]), " outliers in PC",y)
}
foundOut = unlist(foundOut)
if(length(foundOut) > 0){
message("\rDetected ", length(unique(foundOut)), " unique population outliers!");
outlier_df = read_in[read_in$IID %in% foundOut, ];
filename = gsub(".pca", "", pca_files[[x]]);
fwrite(outlier_df[,colnames(outlier_df) %in% c("FID", "IID")], file=paste("OUTLIER_",filename,".txt", sep=""), col.names=F, quote=F,row.names=F,sep="\t")
}
}' | R --vanilla

cat OUTLIER*founder*txt > OUTLIER_all.txt

plink --bfile  ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial \
--remove OUTLIER_all.txt \
--make-bed --out ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial_pca

done



#-->converting files to vcf format
echo "converting plink to vcf format"
bim_files=`ls /media/sf_ad_resilience/niagads/genotype/out_for_HRC_imputation/*.bim`

for iter in ${bim_files}; do
echo $iter
it=`basename ${iter%%.*}`
echo "working on:" ${it}

for chrom in $(seq 1 22); do
plink --bfile ${new_pfx}_${it}_preimpqc_mind_bi_align_NOdup_initial_pca \
--chr ${chrom} --recode vcf-iid bgz \
--biallelic-only strict \
--out ${new_pfx}_${it}_${chrom}_preimpqc

tabix -f -p vcf ${new_pfx}_${it}_${chrom}_preimpqc.vcf.gz
done
done


#-->strand flipping
cd /media/sf_ad_resilience_extend/HRC/strand_flip_other_than_mirw_and_GenADA
bim_files=`ls /media/sf_ad_resilience/niagads/genotype/out_for_HRC_imputation/*.bim`
export PATH=~/Documents:$PATH
new_pfx="FO"
path1='/media/sf_ad_resilience_extend/preimp'
path2='/media/sf_ad_resilience/niagads/strand_flip_check'
for iter in ${bim_files}; do
echo $iter
it=`basename ${iter%%.*}`
echo "working on:" ${it}

for chrom in $(seq 1 22); do

java -jar ${path2}/conform-gt.24May16.cee.jar ref=${path2}/chr${chrom}.1kg.phase3.v5a.vcf.gz gt=${path1}/${new_pfx}_${it}_${chrom}_preimpqc.vcf.gz chrom=${chrom} match=POS out=con_${new_pfx}_${it}_${chrom}_preimpqc excludesamples=${path2}/non.eur.excl

done
done





