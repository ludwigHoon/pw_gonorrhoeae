library(minSNPs)
library(data.table)
library(BiocParallel)

run_analysis <- function(metadata, snp_matrix){
    BP <- MulticoreParam(workers = 120)
    # Intermediate resistance is considered as resistance
    for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
        metadata[[drug]] <- ifelse(metadata[[drug]] == "INTERMEDIATE", "RESISTANT", metadata[[drug]])
    }

    GOI <- metadata[Cefixime == "RESISTANT"]$Genome_Name
    timestamp("--- Running Percent ---")
    # Percent_mode analysis for Cefixime
    percent_result <- find_optimised_snps(seqc = snp_matrix, metric = "percent", goi = GOI,
        number_of_result = 10, max_depth = 5, bp = BP,
        output_progress = TRUE)
    saveRDS(percent_result, "Cefixime_percent_AU_sample.rds")
    try(output_result(percent_result, view = "tsv", file_name = "Cefixime_percent_AU_sample.tsv", seqc = snp_matrix))
    
    timestamp("--- Running MCC ---")
    # MCC mode analysis for Cefixime
    mcc_result <- find_optimised_snps(seqc = snp_matrix, metric = "mcc", goi = GOI,
        number_of_result = 10, max_depth = 5, bp = BP,
        output_progress = TRUE)
    saveRDS(mcc_result, "Cefixime_mcc_AU_sample.rds")
    try(output_result(mcc_result, view = "tsv", file_name = "Cefixime_mcc_AU_sample.tsv", seqc = snp_matrix))

    timestamp("--- Completed ---")
}

AU_matrix <- read_fasta("AU_final_matrix.fasta")
metadata <- fread("no_dup_metadata.csv")
for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
    metadata[[drug]] <- ifelse(metadata[[drug]] == "INTERMEDIATE", "RESISTANT", metadata[[drug]])
}
AU_metadata <- metadata[metadata$Country == "Australia" | Genome_Name == "Reference",]
run_analysis(AU_metadata, AU_matrix)
