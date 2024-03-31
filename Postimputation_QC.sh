cd  /media/sf_ad_resilience_extend/HRC/post_imputation
prefix=("Omni" )

path="/media/sf_ad_resilience_extend/HRC/imputation"
export PATH=~/Documents:$PATH
export PATH=~/Documents/plink_linux_i686_20201019:$PATH

for it in "${prefix[@]}"; do
echo "working on" $it
for chrom in $(seq 1 22); do


plink --vcf  ${path}/${it}_chr${chrom}.dose.vcf.gz \
--double-id \
--biallelic-only strict \
--geno .02 --mind .02 --maf 1e-2 \
--make-bed --out s1


mv s1.log ${it}_chr${chrom}_s1.log


#-->Remove duplicates variant IDs from GWAS sample
echo "." > plink.dupvar 
#When multiple variants share the same bp coordinate and allele codes, it is likely that they are not actually distinct, and the duplicates should be merged or removed.
plink --bfile s1 --list-duplicate-vars ids-only suppress-first
plink --bfile s1 --exclude plink.dupvar --make-bed --out s3

mv s3.log ${it}_chr${chrom}_s3.log

#-->Find unique variants 
echo "." > OutInfoDupVars.txt

echo '
require(data.table)
readInfo = fread("s3.bim",header = F, sep = "\t")
cat("\rDetected",nrow(readInfo),"total markers")
snpid = readInfo$V2
snpid = snpid[duplicated(snpid)]
cat("\rDetected",length(snpid),"duplicate markers!")
if(length(snpid) > 0){
fwrite(data.table(snpid), file="OutInfoDupVars.txt",quote= F,row.names=F,col.names=F, sep = "\t")}
' | R --vanilla


plink --bfile s3 \
--exclude OutInfoDupVars.txt \
--make-bed --out s3r

#-->Remove duplicate variant IDs from qual score file
mv s3r.log ${it}_chr${chrom}_s3r.log

echo '
require(data.table)
require(R.utils)
readInfo = fread("/media/sf_ad_resilience_extend/HRC/imputation/'${it}_chr${chrom}.info'.gz") #can not delete the single quotes
print(head(readInfo))
cat("\rDetected",nrow(readInfo),"total markers")
snpid = readInfo$SNP
snpid = snpid[duplicated(snpid)]
readInfo = readInfo[!readInfo$SNP %in% snpid]
print(head(readInfo))
print(nrow(readInfo))
cat("\rDetected",length(snpid),"duplicate markers!")
if(nrow(readInfo) > 0){
fwrite(data.table(readInfo), file="INFO_dupvars.txt",quote= F,row.names=F,col.names=T, sep = "\t")}
' | R --vanilla


#-->Exclude markers with imputation score (INFO) < threshold
plink --bfile s3r \
--qual-scores INFO_dupvars.txt 7 1 1 \
--qual-threshold 0.8  \
--make-bed --out POSTIMP_${it}_${chrom}

done
done



#-->Update SNP name
cd  /media/sf_ad_resilience_extend/HRC/post_imputation
prefix=("Omni" )


path="/media/sf_ad_resilience_extend/post_imputation_qc"
export PATH=~/Documents:$PATH
export PATH=~/Documents/plink_linux_i686_20201019:$PATH

for it in "${prefix[@]}"; do
echo "working on" $it
for chrom in $(seq 1 22); do
#defalt: the new value is read from column 2 and the (old) variant ID from column 1
plink --bfile POSTIMP_${it}_${chrom} \
--update-name ${path}/snp_ForConvert_${chrom}.txt \
--make-bed --out POSTIMP_${it}_${chrom}_updated
done
done



#-->removes in-del variants 
cd  /media/sf_ad_resilience_extend/HRC/post_imputation
prefix=("Omni" )
export PATH=~/Documents:$PATH
export PATH=~/Documents/plink_linux_i686_20201019:$PATH

for it in "${prefix[@]}"; do
echo "working on" $it
for chrom in $(seq 1 22); do

echo '
require(data.table)
bim_files = list.files(pattern="'POSTIMP_${it}_${chrom}_updated'.bim")
head(bim_files)
bim = fread(bim_files)
bim_detect = bim[nchar(bim$V5) >= 2 | nchar(bim$V6) >= 2, ]
fwrite(data.table(snp = bim_detect$V2), file="'POSTIMP_${it}_${chrom}_indel'.txt", quote=F,sep="\t",col.names=FALSE,row.names=FALSE)
'  | R --vanilla

#removes all listed variants from the current analysis.
plink --bfile POSTIMP_${it}_${chrom}_updated \
--exclude POSTIMP_${it}_${chrom}_indel.txt \
--make-bed --out POSTIMP_${it}_${chrom}_updated_noINDEL
done
done



#-->Merge chromosome
cd  /media/sf_ad_resilience_extend/HRC/post_imputation
prefix=("Omni")
export PATH=~/Documents:$PATH
export PATH=~/Documents/plink_linux_i686_20201019:$PATH

for it in "${prefix[@]}"; do
echo "working on" $it

rm fileset.txt
for chrom in $(seq 1 22); do


#Write a list of the 22 autosomal filesets to be merged
echo POSTIMP_${it}_${chrom}_updated_noINDEL >> fileset.txt
done

#then merge them
plink --merge-list fileset.txt --make-bed --out POSTIMP_${it}_merged
done



#-->Exclude multiallelic SNPs
cd  /media/sf_ad_resilience_extend/HRC/post_imputation
prefix=("Omni")

export PATH=~/Documents:$PATH
export PATH=~/Documents/plink_linux_i686_20201019:$PATH

for it in "${prefix[@]}"; do
echo "working on" $it


echo '
require(data.table)
bim = fread("'POSTIMP_${it}_merged'.bim")
bim$chr_pos = paste(bim$V1, "_",bim$V4, sep = "")
bim_du = bim[duplicated(bim$chr_pos) == T,]
su = bim[bim$chr_pos %in% bim_du$chr_pos,]

fwrite(data.table(snp = su$V2), file="'${it}_sus_multi'.txt", quote=F,sep="\t",col.names=FALSE,row.names=FALSE)' | R --vanilla


plink --bfile POSTIMP_${it}_merged \
--exclude ${it}_sus_multi.txt \
--make-bed --out POSTIMP_${it}_merged_d
done


