for P1 in results/*_P1_*.vcf
do
P2=${P1//_P1_/_P2_}
P1g=${P1//.vcf/.vcf.gz}
P2g=${P2//.vcf/.vcf.gz}
Out=${P1//_P1_/_}
Filt=${Out//.vcf/_filtered}
bgzip $P1
bgzip $P2
tabix -p vcf $P1g
tabix -p vcf $P2g
vcf-merge -R '0|0' $P1g $P2g > $Out
vcftools --vcf $Out --maf 0.1 --recode --recode-INFO-all --out $Filt
done
