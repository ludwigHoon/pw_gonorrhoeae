library(data.table)
library(minSNPs)
library(BiocParallel)

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
SNPs_pos <- fread("FIN_global_SNP_POS.csv")

# Applicability of Multidrugs resistant SNP sets result in all samples, irrespective of which subset of samples used to derive the set
run_binary_generalisation <- function(target_class){
    mod_metadata <- metadata
    all_mcc_result_f <- list.files(pattern = paste0(target_class, "_mcc_global_.*.tsv"))
    all_percent_result_f <- list.files(pattern = paste0(target_class, "_percent_global_.*.tsv"))
    results <- lapply(c(all_mcc_result_f, all_percent_result_f), process_result_file)
    result_type <- c(rep("mcc", length(all_mcc_result_f)), rep("percent", length(all_percent_result_f)))
    set_source <- gsub(".tsv", "", sapply(strsplit(c(all_mcc_result_f, all_percent_result_f), split = "_"), tail, n = 1))
    names(results) <- paste(result_type, set_source, sep = "_")
    
    priority <- data.frame(target = c("GOI", "NON_GOI"), priority = c(1, 2))
    all_snp_sets <- unique(unlist(results))
    drugs <- strsplit(target_class, split = "_")[[1]]
    mod_metadata$target <- ifelse(mod_metadata[[drugs[1]]] == "RESISTANT" & mod_metadata[[drugs[2]]] == "RESISTANT", "GOI", "NON_GOI")
    colnames(mod_metadata)[which(colnames(mod_metadata) == "Genome_Name")] <- "isolate"
    
    results_l <- list()
    pb <- txtProgressBar(min = 0, max = length(names(results)), style = 3)
    for (ssource in names(results)){
        setTxtProgressBar(pb, pb$getVal() + 1)
        metadata <- mod_metadata
        metric <- unlist(strsplit(ssource, "_"))[1]
        oset <- unlist(strsplit(ssource, "_"))[2]
        O_metadata <- metadata[ set == oset,]
        N_metadata <- metadata[!isolate %in% O_metadata$isolate,]
        O_seq <- global_matrix_full[names(global_matrix_full) %in% O_metadata$isolate]
        N_seq <- global_matrix_full[names(global_matrix_full) %in% N_metadata$isolate]
        stopifnot(length(O_seq) + length(N_seq) == length(global_matrix_full))
        sset <- results[[ssource]]
        result <- summarise_result(snp_sets = sset,
            training_seqs = O_seq, validation_seqs = N_seq,
            training_metadata = O_metadata,
            validation_metadata = N_metadata,
            priority = priority, is_multi = FALSE, return_all_intermediate = TRUE,
            is_percent = if (metric == "percent") {
                TRUE
            } else {
                FALSE
            }
        )
        result$type <- metric
        result$set <- oset
        results_l[[ssource]] <- result
    }
    close(pb)
    result_dt <- rbindlist(results_l, use.names = TRUE )
    result_dt$snp_sets <- sapply(result_dt$snp_sets, function(sset){
        sset <- as.numeric(unlist(strsplit(sset, split = ", ")))
        paste0(SNPs_pos[match(sset, SNPs_pos$Fasta_Position)]$Reference_Position, collapse = ", ")
    })
    fwrite(result_dt, paste0("binary_", target_class, "_global.csv"), row.names = FALSE)
}

## Binary class
### Ciprofloxacin + Penicillin
run_binary_generalisation("Ciprofloxacin_Penicillin")

### Cefixime + Ciprofloxacin
run_binary_generalisation("Cefixime_Ciprofloxacin")


run_multiclass_generalisation <- function(target_class){
    mod_metadata <- metadata
    for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
        mod_metadata[[drug]] <- ifelse(mod_metadata[[drug]] == "RESISTANT", "R", "S")
    }
    all_mcc_multi_result_f <- list.files(pattern = paste0(target_class, "_mcc_multi_global_.*.tsv"))
    all_simpson_by_group_result_f <- list.files(pattern = paste0(target_class, "_simpson_by_group_global_.*.tsv"))
    results <- lapply(c(all_mcc_multi_result_f, all_simpson_by_group_result_f), process_result_file)
    result_type <- c(rep("mcc_multi", length(all_mcc_multi_result_f)), rep("simpson_by_group", length(all_simpson_by_group_result_f)))
    set_source <- gsub(".tsv", "", sapply(strsplit(c(all_mcc_multi_result_f, all_simpson_by_group_result_f), split = "_"), tail, n = 1))
    names(results) <- paste(result_type, set_source, sep = "_")

    drugs <- strsplit(target_class, split = "_")[[1]]
    mod_metadata$target <- paste(mod_metadata[[drugs[1]]], mod_metadata[[drugs[2]]], sep = "_")
    colnames(mod_metadata)[which(colnames(mod_metadata) == "Genome_Name")] <- "isolate"
    
    results_l <- list()
    pb <- txtProgressBar(min = 0, max = length(names(results)), style = 3)
    for (ssource in names(results)){
        setTxtProgressBar(pb, pb$getVal() + 1)
        metadata <- mod_metadata
        split_ssource <- unlist(strsplit(ssource, "_"))
        oset <- tail(split_ssource, n = 1)
        metric <- paste(split_ssource[! split_ssource %in% oset], collapse = "_")
        
        O_metadata <- metadata[ set == oset,]
        N_metadata <- metadata[!isolate %in% O_metadata$isolate,]
        O_seq <- global_matrix_full[names(global_matrix_full) %in% O_metadata$isolate]
        N_seq <- global_matrix_full[names(global_matrix_full) %in% N_metadata$isolate]
        stopifnot(length(O_seq) + length(N_seq) == length(global_matrix_full))
        sset <- results[[ssource]]
        priority <- generate_prioritisation(O_metadata)
        result <- summarise_result(snp_sets = sset,
            training_seqs = O_seq, validation_seqs = N_seq,
            training_metadata = O_metadata,
            validation_metadata = N_metadata,
            priority = priority, is_multi = TRUE
        )
        result$type <- metric
        result$set <- oset
        results_l[[ssource]] <- result
    }
    close(pb)
    result_dt <- rbindlist(results_l, use.names = TRUE )
    result_dt$snp_sets <- sapply(result_dt$snp_sets, function(sset){
        sset <- as.numeric(unlist(strsplit(sset, split = ", ")))
        paste0(SNPs_pos[match(sset, SNPs_pos$Fasta_Position)]$Reference_Position, collapse = ", ")
    })
    fwrite(result_dt, paste0("multiclass_", target_class, "_global.csv"), row.names = FALSE)
}

## Multiple class
### Ciprofloxacin + Penicillin
run_multiclass_generalisation("Ciprofloxacin_Penicillin")


### Cefixime + Ciprofloxacin
run_multiclass_generalisation("Cefixime_Ciprofloxacin")


