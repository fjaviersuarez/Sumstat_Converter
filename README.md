# Sumstat_Converter
Sumstat_Converter is a R script to convert harmonized GWAS summary statistics into .assoc.logistic structure to be usable in PRS.


## Usage

```bash
Rscript Sumstat_Converter --base Base_File.assoc.logistic --target 
Summary_Statistics.txt --out Adapted_Summary_Statistics
```
### --base: 
It corresponds to the .assoc.logistic file that will be used by the script as a basis to rewrite the files.

### --target: 
It corresponds to the uncompressed [GWAS Catalog](https://www.ebi.ac.uk/gwas/home) file (with gunzip) you want to adapt to yours.

### --out:
The output **without extension**, the script already adds the necessary output to each file (.log, .prob, and .meta for the meta-analysis results, and .txt for the final output file generated with the corrected SNP and file).



## How it works

Sumstat_Converter selects the necessary and essential columns from the summary statistics for PRS, ordering and renaming them headers to be able to use the summary statistics in comparative studies. It also recreates the SNP column by concatenating the other columns and, through a meta-analysis using [plink1.9](https://www.cog-genomics.org/plink/) to generate a .prob file with problematic SNPs, corrects SNPs that do not match those in your base file. In this way, it solves SNPs with strand changes, alleles swapped, and with these two options together.
## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.
## Citation

If you find this project useful for your work, please consider citing it using the following BibTeX entry:

```bibtex
@misc{Sumstat_Converter,
  title = {Sumstat_Converter},
  author = {Suarez, F.},
  year = {2024},
  howpublished = {\url{https://github.com/fjaviersuarez/Sumstat_Converter}},
}
```
Thank you for acknowledging my work!

