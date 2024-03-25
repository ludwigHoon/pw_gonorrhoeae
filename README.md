# *N. Gonorrhoeae* SNP mining with minSNPs
Data source: [Pathogen**watch**](https://pathogen.watch/genomes/all?organismId=485)
- Assembled contigs: [https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Neisseria%20gonorrhoeae__fastas.zip](https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Neisseria%20gonorrhoeae__fastas.zip)
The zip file containing all the assembled contigs need to be unzipped with `7zz e -aou`, duplicated assembly can be found in the dataset, when duplicated, "_1" is added to the file name.
- Metadatas: 
    - [https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Neisseria%20gonorrhoeae__metadata.csv.gz](https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Neisseria%20gonorrhoeae__metadata.csv.gz)
    - [https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Neisseria%20gonorrhoeae__mlst.csv.gz](https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Neisseria%20gonorrhoeae__mlst.csv.gz)
    - [https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Neisseria%20gonorrhoeae__paarsnp.csv.gz](https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/Neisseria%20gonorrhoeae__paarsnp.csv.gz)


# Combining the metadata & creating SNP matrix
- See `Scripts/snippy_commands.sh` for commands that align each of samples with the reference genome and call the SNPs
- See `Scripts/process_metadata.R` for combining the metadata and filtering the samples
    - Run the resulting `snippy-core_AU.sh` and `snippy-core_global.sh`,and the following once the snippy-core is done:
    - awk '{print $2}' snippy_output/AU.tab > AU_SNP_POS.txt
    - awk '{print $2}' snippy_output/global.tab > global_SNP_POS.txt
- See `Scripts/process_snp_matrix.R` for processing and extracting unique sequences in SNP matrix.

Final Matrices can be accessed at [dataverse](https://doi.org/10.7910/DVN/K2AECH), note SNP matrix containing all samples is split into 5 files due to the size.
For SNP position reference, see `Data/FIN_AU_SNP_POS.csv` and `Data/FIN_global_SNP_POS.csv`

## Some statistics
Reference genome: TUM19854
Original number of samples: 38381
Actual number of samples: 31076

Australian samples:
- Actual samples: 2416
- Number of files: 4800
- Number of samples with duplicates: 2380

Global samples:
- Actual samples: 31076
- Number of files: 38287
- Number of samples with duplicates: 7211

# Generate phylogenetic tree
- AU only: `fasttree -nt AU_final_matrix.fasta > AU_only.tre`
- Global: `fasttree -nt global_final_matrix.fasta > global.tre`

# Running minSNPs
## Cefixime resistance
- AU only:
    - Script: `Scripts\AU_Cefixime.R` and `Scripts\AU_Cefixime2global.R`
    - Results: `Result/Cefixime_percent_AU_sample.tsv`, `Result/Cefixime_mcc_AU_sample.tsv`, `Result\Cefixime_AU2global.csv`
- Global sample:
    - Script: `Scripts\global_Cefixime.R` and `Scripts\global_Cefixime2All.R`
    - Results: `Result\Cefixime_mcc_global_<n>.tsv`, `Result\Cefixime_percent_global_<n>.tsv`, `Result\Cefixime_set2global.csv`, with `<n>` represent the subset ID (1-5) used to identify the SNP sets.
