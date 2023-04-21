# SNP filtering

Look at missingness in unfiltdered vcf file:

```
vcftools --vcf unfiltered.vcf --missing-indv
```

Remove all individuals with <70% missing data using samplesname of samples I want to remove:

```
vcftools --vcf unfiltered.vcf --remove-indv <sample1> --remove-indv <sample2> --recode
```

Rename recoded vcf:

```
mv out.recode.vcf filtered_remove-indv.vcf
```

Filter for depth and allow SNPs with up to 20% of missing data

```
vcftools --vcf filtered_remove-indv.vcf --minDP 3 --maxDP 20 --max-missing 0.8 --recode
```

Rename recoded vcf:

```
mv out.recode.vcf fully_filtered.vcf
```


Fully filtered vcf file can now be used for analysis. 

vcftools was also used to subset the filtered dataset (e.g., for mainland/island birds or different populations)