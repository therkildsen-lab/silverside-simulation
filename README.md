# Simulation for the Atlantic Silverside system

Code written for SLiM 3.5 which specifies two highly connected populations, with set inversions arising at different time points

Output from Silverside_inversion_mutEnabled.slim (two vcf files) were merged with vcftools function vcf-merge
```
vcf-merge -R "0|0" <seed>_P1_SilverSide_Inversion.vcf.gz <seed>_P2_SilverSide_Inversion.vcf.gz > <seed>_Silverside_Inversion2.vcf
```