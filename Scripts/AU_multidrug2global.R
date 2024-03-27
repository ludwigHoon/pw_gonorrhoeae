library(data.table)
library(minSNPs)
library(BiocParallel)

# Applicability of Multidrugs resistant SNP sets result in all samples, irrespective of country of origin

run_binary_generalisation <- function(target_class){
    percent_result <-paste0(target_class, "_percent_AU_sample.tsv")
    mcc_result <- paste0(target_class, "_mcc_AU_sample.tsv")
    results <- lapply(c(percent_result, mcc_result), process_result_file)
    names(results) <- c("percent", "mcc")

    # Convert SNPs to reference positions
    SNPs_pos <- fread("FIN_AU_SNP_POS.csv")
    REF_GENOME_pos <- data.table(selected_snps = unique(unlist(results)),
        reference_position = SNPs_pos[match(unique(unlist(results)), Fasta_Position)]$Reference_Position)[order(selected_snps)]
    # Scan each of the SNPs in the snippy temp output
    metadata <- fread("no_dup_metadata.csv")
    for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
        metadata[[drug]] <- ifelse(metadata[[drug]] == "INTERMEDIATE", "RESISTANT", metadata[[drug]])
    }

    res <- bplapply(metadata$Genome_Name, function(genome, positions){
        library(minSNPs)
        seq <- read_fasta(paste0("snippy_output/", genome, "/snps.aligned.fa"))
        pat <- generate_pattern(seq, positions)
        stopifnot(length(pat) == 1)
        return(unname(unlist(pat)))
    }, positions = REF_GENOME_pos$reference_position, BPPARAM =
        BatchtoolsParam(workers = 300, cluster="slurm", template="~/slurm_bc_template.tmpl",
                resources=list(walltime=60*60, ncpus=1), progress = TRUE))
    names(res) <- metadata$Genome_Name
    
    drugs <- strsplit(target_class, split = "_")[[1]]
    metadata$target <- ifelse(metadata[[drugs[1]]] == "RESISTANT" & metadata[[drugs[2]]] == "RESISTANT", "GOI", "NON_GOI")
    priority <- data.frame(target = c("GOI", "NON_GOI"), priority = c(1, 2))
    colnames(metadata)[which(colnames(metadata) == "Genome_Name")] <- "isolate"
    AU_metadata <- metadata[metadata$Country == "Australia" | isolate == "Reference",]
    VAL_metadata <- metadata[!isolate %in% AU_metadata$isolate,]
    AU_only <- res[names(res) %in% AU_metadata$isolate]
    VAL_seq <- res[!names(res) %in% names(AU_only)]
    VAL_seq <- setNames(strsplit(unlist(VAL_seq), split = ""), names(VAL_seq))
    AU_only <- setNames(strsplit(unlist(AU_only), split = ""), names(AU_only))
    stopifnot(all(names(VAL_seq) == VAL_metadata$isolate))

    # Test MCC SNP sets
    mcc_snp_sets <- lapply(results$mcc, match, REF_GENOME_pos$selected_snps)
    PERF_mcc_snp_sets <- summarise_result(snp_sets = mcc_snp_sets, training_seqs = AU_only, validation_seqs = VAL_seq, training_metadata = AU_metadata, validation_metadata = VAL_metadata, priority = priority, is_multi = FALSE, return_all_intermediate = TRUE, is_percent = FALSE)
    PERF_mcc_snp_sets$type <- "mcc"
    PERF_mcc_snp_sets$snp_sets <- sapply(PERF_mcc_snp_sets$snp_sets, function(sset){
        sset <- as.numeric(unlist(strsplit(sset, split = ", ")))
        paste0(REF_GENOME_pos[sset]$reference_position, collapse = ", ")
    })

    # Test percent SNP sets
    pct_snp_sets <- lapply(results$percent, match, REF_GENOME_pos$selected_snps)
    PERF_pct_snp_sets <- summarise_result(snp_sets = pct_snp_sets, training_seqs = AU_only, validation_seqs = VAL_seq, training_metadata = AU_metadata, validation_metadata = VAL_metadata, priority = priority, is_multi = FALSE, return_all_intermediate = TRUE, is_percent = TRUE)
    PERF_pct_snp_sets$type <- "percent"
    PERF_pct_snp_sets$snp_sets <- sapply(PERF_pct_snp_sets$snp_sets, function(sset){
        sset <- as.numeric(unlist(strsplit(sset, split = ", ")))
        paste0(REF_GENOME_pos[sset]$reference_position, collapse = ", ")
    })

    all_results <- rbind(PERF_mcc_snp_sets, PERF_pct_snp_sets)
    fwrite(all_results, paste0("binary_", target_class, "_AU2global.csv"), row.names = FALSE)
}


## Binary class
### Ciprofloxacin + Penicillin
run_binary_generalisation("Ciprofloxacin_Penicillin")

### Cefixime + Ciprofloxacin
run_binary_generalisation("Cefixime_Ciprofloxacin")


run_multiclass_generalisation <- function(target_class, get_subclass = FALSE){
    mcc_multi_result <-paste0(target_class, "_percent_AU_sample.tsv")
    simpson_by_group_result <- paste0(target_class, "_mcc_AU_sample.tsv")
    results <- lapply(c(mcc_multi_result, simpson_by_group_result), process_result_file)
    names(results) <- c("mcc_multi", "simpson_by_group")

    # Convert SNPs to reference positions
    SNPs_pos <- fread("FIN_AU_SNP_POS.csv")
    REF_GENOME_pos <- data.table(selected_snps = unique(unlist(results)),
        reference_position = SNPs_pos[match(unique(unlist(results)), Fasta_Position)]$Reference_Position)[order(selected_snps)]
    # Scan each of the SNPs in the snippy temp output
    metadata <- fread("no_dup_metadata.csv")
    for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
        metadata[[drug]] <- ifelse(metadata[[drug]] == "INTERMEDIATE", "RESISTANT", metadata[[drug]])
    }
    for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
        metadata[[drug]] <- ifelse(metadata[[drug]] == "RESISTANT", "R", "S")
    }

    res <- bplapply(metadata$Genome_Name, function(genome, positions){
        library(minSNPs)
        seq <- read_fasta(paste0("snippy_output/", genome, "/snps.aligned.fa"))
        pat <- generate_pattern(seq, positions)
        stopifnot(length(pat) == 1)
        return(unname(unlist(pat)))
    }, positions = REF_GENOME_pos$reference_position, BPPARAM =
        BatchtoolsParam(workers = 300, cluster="slurm", template="~/slurm_bc_template.tmpl",
                resources=list(walltime=60*60, ncpus=1), progress = TRUE))
    names(res) <- metadata$Genome_Name

    drugs <- strsplit(target_class, split = "_")[[1]]
    mod_metadata <- metadata
    mod_metadata$target <- paste(mod_metadata[[drugs[1]]], mod_metadata[[drugs[2]]], sep = "_")
    colnames(mod_metadata)[which(colnames(mod_metadata) == "Genome_Name")] <- "isolate"
    AU_metadata <- mod_metadata[mod_metadata$Country == "Australia" | isolate == "Reference",]
    VAL_metadata <- mod_metadata[!isolate %in% AU_metadata$isolate,]
    AU_only <- res[names(res) %in% AU_metadata$isolate]
    VAL_seq <- res[!names(res) %in% names(AU_only)]
    VAL_seq <- setNames(strsplit(unlist(VAL_seq), split = ""), names(VAL_seq))
    AU_only <- setNames(strsplit(unlist(AU_only), split = ""), names(AU_only))
    stopifnot(all(names(VAL_seq) == VAL_metadata$isolate))
    priority <- generate_prioritisation(AU_metadata)

    mcc_multi_snp_sets <- lapply(results$mcc_multi, match, REF_GENOME_pos$selected_snps)
    PERF_mcc_multi_snp_sets <- summarise_result(snp_sets = mcc_multi_snp_sets, training_seqs = AU_only,
        validation_seqs = VAL_seq, training_metadata = AU_metadata, validation_metadata = VAL_metadata,
        priority = priority, is_multi = (!get_subclass))
    PERF_mcc_multi_snp_sets$type <- "mcc_multi"
    PERF_mcc_multi_snp_sets$snp_sets <- sapply(PERF_mcc_multi_snp_sets$snp_sets, function(sset){
        sset <- as.numeric(unlist(strsplit(sset, split = ", ")))
        paste0(REF_GENOME_pos[sset]$reference_position, collapse = ", ")
    })

    simpson_by_group_snp_sets <- lapply(results$simpson_by_group, match, REF_GENOME_pos$selected_snps)
    PERF_simpson_by_group_snp_sets <- summarise_result(snp_sets = simpson_by_group_snp_sets,
        training_seqs = AU_only, validation_seqs = VAL_seq, training_metadata = AU_metadata,
        validation_metadata = VAL_metadata, priority = priority, is_multi = (!get_subclass))
    PERF_simpson_by_group_snp_sets$type <- "simpson_by_group"
    PERF_simpson_by_group_snp_sets$snp_sets <- sapply(PERF_simpson_by_group_snp_sets$snp_sets, function(sset){
        sset <- as.numeric(unlist(strsplit(sset, split = ", ")))
        paste0(REF_GENOME_pos[sset]$reference_position, collapse = ", ")
    })

    all_results <- rbind(PERF_mcc_multi_snp_sets, PERF_simpson_by_group_snp_sets)
    if (get_subclass){
        fwrite(all_results, paste0("multiclass_", target_class, "_subclass_AU2global.csv"), row.names = FALSE)
    } else {
        fwrite(all_results, paste0("multiclass_", target_class, "_AU2global.csv"), row.names = FALSE)
    }
}

## Multiple class
### Ciprofloxacin + Penicillin
run_multiclass_generalisation("Ciprofloxacin_Penicillin")


### Cefixime + Ciprofloxacin
run_multiclass_generalisation("Cefixime_Ciprofloxacin")

### Ciprofloxacin + Penicillin
run_multiclass_generalisation("Ciprofloxacin_Penicillin", get_subclass = TRUE)


### Cefixime + Ciprofloxacin
run_multiclass_generalisation("Cefixime_Ciprofloxacin", get_subclass = TRUE)