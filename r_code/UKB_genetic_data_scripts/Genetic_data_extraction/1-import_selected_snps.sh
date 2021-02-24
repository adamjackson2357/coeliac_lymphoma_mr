module load plink

UKB=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018

for chr in $(seq 1 22)
do
echo $chr

plink --bfile $UKB/ukb_imp_chr$chr --fam plink_files/ukb_imp_19266.fam --extract snp_list.txt --make-bed --out plink_files/ukb_imp_chr$chr

done


