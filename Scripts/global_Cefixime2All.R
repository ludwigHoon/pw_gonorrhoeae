library(minSNPs)
library(data.table)

global_matrix_full <- read_fasta("global_final_matrix.fasta")

metadata <- fread("no_dup_metadata.csv")
for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
    metadata[[drug]] <- ifelse(metadata[[drug]] == "INTERMEDIATE", "RESISTANT", metadata[[drug]])
}
samples_in_set <- lapply(c(1:5), function(i){
    samples <- system(paste0("grep \">\" global_final_matrix_", i, ".fasta"), intern = TRUE)
    samples <- gsub(">", "", samples)
    return(samples)
})
metadata$set <- sapply(metadata$Genome_Name, function(s){
    inset <- which(sapply(samples_in_set, function(sset){s %in% sset}))
    return(paste0(inset, collapse = ","))
})

# Applicability of result in all samples, irrespective of whether the samples have been used to derive the SNP sets
all_results_f <- list.files(pattern = paste0("Cefixime_.*_global_.*.tsv"))
results <- lapply(all_results_f, process_result_file)
set_source <- sapply(all_results_f, function(r_f){
    metric <- unlist(strsplit(r_f, "_"))[2]
    set <- gsub(".tsv", "", unlist(strsplit(r_f, "_"))[4])
    return(paste0(metric, "_", set))
})
names(results) <- set_source
all_snp_sets <- unique(unlist(results))
SNPs_pos <- fread("FIN_global_SNP_POS.csv")

priority <- data.frame(target = c("GOI", "NON_GOI"), priority = c(1, 2))
colnames(metadata)[which(colnames(metadata) == "Cefixime")] <- "target"
colnames(metadata)[which(colnames(metadata) == "Genome_Name")] <- "isolate"
metadata$target <- ifelse(metadata$target == "RESISTANT", "GOI", "NON_GOI")

pb <- txtProgressBar(min = 0, max = length(set_source), style = 3)
all_summarised_results <- list()
for (ssource in set_source){
    metric <- unlist(strsplit(ssource, "_"))[1]
    oset <- unlist(strsplit(ssource, "_"))[2]

    O_metadata <- metadata[ set == oset,]
    N_metadata <- metadata[!isolate %in% O_metadata$isolate,]
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    
    O_seq <- global_matrix_full[names(global_matrix_full) %in% O_metadata$isolate]
    N_seq <- global_matrix_full[names(global_matrix_full) %in% N_metadata$isolate]
    stopifnot(length(O_seq) + length(N_seq) == length(global_matrix_full))

    sset <- results[[ssource]]
    sr_source <- summarise_result(snp_sets = sset,
        training_seqs = O_seq, validation_seqs = N_seq,
        training_metadata = O_metadata,
        validation_metadata = N_metadata,
        priority = priority, is_multi = FALSE, return_all_intermediate = TRUE,
        is_percent = if (metric == "percent") {
            TRUE
        } else {
            FALSE
        })
    sr_source$type <- metric
    sr_source$set <- oset
    all_summarised_results[[ssource]] <- sr_source
}
close(pb)
all_summarised_results <- rbindlist(all_summarised_results)
all_summarised_results$snp_sets <- sapply(all_summarised_results$snp_sets, function(s){
    fasta_positions <- strsplit(s, split = ", ")[[1]]
    reference_positions <- SNPs_pos[match(fasta_positions, SNPs_pos$Fasta_Position),]$Reference_Position
    return(paste0(reference_positions, collapse = ", "))
})
fwrite(all_summarised_results, "Cefixime_set2global.csv", row.names = FALSE)