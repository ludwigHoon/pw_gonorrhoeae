# AU Multidrug - Cefixime + Ciprofloxacin | Ciprofloxacin + Penicillin
library(minSNPs)
library(data.table)
library(BiocParallel)

# MCC + Percent
run_binary_analysis <- function(GOI, snp_matrix, target_name){
    BP <- MulticoreParam(workers = 120)

    print(paste("--- Running Percent --- at:", Sys.time()))
    # Percent_mode analysis
    percent_result <- find_optimised_snps(seqc = snp_matrix, metric = "percent", goi = GOI,
        number_of_result = 10, max_depth = 5, bp = BP,
        output_progress = TRUE)
    saveRDS(percent_result, paste0(target_name, "_percent_AU_sample.rds"))
    try(output_result(percent_result, view = "tsv", file_name = paste0(target_name, "_percent_AU_sample.tsv"), seqc = snp_matrix))
    
    print(paste("--- Running MCC --- at:", Sys.time()))
    # MCC mode analysis
    mcc_result <- find_optimised_snps(seqc = snp_matrix, metric = "mcc", goi = GOI,
        number_of_result = 10, max_depth = 5, bp = BP,
        output_progress = TRUE)
    saveRDS(mcc_result, paste0(target_name, "_mcc_AU_sample.rds"))
    try(output_result(mcc_result, view = "tsv", file_name = paste0(target_name, "_mcc_AU_sample.tsv"), seqc = snp_matrix))

    print(paste("--- Completed --- at:", Sys.time()))
}

# MCC_Multi + Simpson_by_group
run_multi_analysis <- function(metadata, snp_matrix, target_name){
    BP <- MulticoreParam(workers = 120)

    print(paste("--- Running mcc_multi --- at:", Sys.time()))
    # mcc_multi analysis
    percent_result <- find_optimised_snps(seqc = snp_matrix, metric = "mcc_multi",
        number_of_result = 10, max_depth = 5, meta = metadata, bp = BP, target = "target",
        output_progress = TRUE)
    saveRDS(percent_result, paste0(target_name, "_mcc_multi_AU_sample.rds"))
    try(output_result(percent_result, view = "tsv", meta = metadata, target = "target", file_name = paste0(target_name, "_mcc_multi_AU_sample.tsv"), seqc = snp_matrix))
    
    print(paste("--- Running simpson_by_group --- at:", Sys.time()))
    # simpson_by_group analysis
    mcc_result <- find_optimised_snps(seqc = snp_matrix, metric = "simpson_by_group",
        number_of_result = 10, max_depth = 5, meta = metadata, bp = BP, target = "target",
        output_progress = TRUE)
    saveRDS(mcc_result, paste0(target_name, "_simpson_by_group_AU_sample.rds"))
    try(output_result(mcc_result, view = "tsv", meta = metadata, target = "target", file_name = paste0(target_name, "_simpson_by_group_AU_sample.tsv"), seqc = snp_matrix))

    print(paste("--- Completed --- at:", Sys.time()))
}

AU_matrix <- read_fasta("AU_final_matrix.fasta")
metadata <- fread("no_dup_metadata.csv")
# Intermediate resistance is considered as resistance
for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
    metadata[[drug]] <- ifelse(metadata[[drug]] == "INTERMEDIATE", "RESISTANT", metadata[[drug]])
}
AU_metadata <- metadata[metadata$Country == "Australia" | Genome_Name == "Reference",]

# Multidrug - binary analysis
## Cefixime + Ciprofloxacin
GOI <- AU_metadata[AU_metadata$Cefixime == "RESISTANT" & AU_metadata$Ciprofloxacin == "RESISTANT"]$Genome_Name
run_binary_analysis(GOI, AU_matrix, "Cefixime_Ciprofloxacin")

## Ciprofloxacin + Penicillin
GOI <- AU_metadata[AU_metadata$Ciprofloxacin == "RESISTANT" & AU_metadata$Penicillin == "RESISTANT"]$Genome_Name
run_binary_analysis(GOI, AU_matrix, "Ciprofloxacin_Penicillin")

# Multidrug - multiclass analysis
for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
    AU_metadata[[drug]] <- ifelse(AU_metadata[[drug]] == "NOT_FOUND", "S", "R")
}
## Cefixime + Ciprofloxacin
mod_metadata <- AU_metadata
mod_metadata$target <- paste(mod_metadata$Cefixime, mod_metadata$Ciprofloxacin, sep = "_")
colnames(mod_metadata)[which(colnames(mod_metadata) == "Genome_Name")] <- "isolate"
run_multi_analysis(mod_metadata, AU_matrix, "Cefixime_Ciprofloxacin")

## Ciprofloxacin + Penicillin
mod_metadata <- AU_metadata
mod_metadata$target <- paste(mod_metadata$Ciprofloxacin, mod_metadata$Penicillin, sep = "_")
colnames(mod_metadata)[which(colnames(mod_metadata) == "Genome_Name")] <- "isolate"
run_multi_analysis(mod_metadata, AU_matrix, "Ciprofloxacin_Penicillin")
