library(data.table)
library(minSNPs)
library(BiocParallel)

# Applicability of result in all samples, irrespective of country of origin
results <- lapply(c("Cefixime_percent_AU_sample.tsv", "Cefixime_mcc_AU_sample.tsv"), process_result_file)
names(results) <- c("percent", "mcc")
SNPs_pos <- fread("FIN_AU_SNP_POS.csv")
REF_GENOME_pos <- data.table(selected_snps = unique(unlist(results)),
    reference_position = SNPs_pos[match(unique(unlist(results)), Fasta_Position)]$Reference_Position)[order(selected_snps)]

# Copying reference fasta file to snippy_output to simplify analysis
system("mkdir snippy_output/Reference")
system("cp reference_TUM19854.fasta snippy_output/Reference/snps.aligned.fa")

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

# Check and ensure that all the samples in Australia SNP matrix is exactly the same
# as the SNPs pulled out from the snippy output
AU_matrix <- read_fasta("AU_final_matrix.fasta")
AU_only <- generate_pattern(AU_matrix, REF_GENOME_pos$selected_snps)
confirmation <- sapply(names(AU_only), function(x){
    return(AU_only[[x]] == res[[x]])
})
stopifnot(all(confirmation) == TRUE) # Something must have gone wrong if this fails


# Test the performance of the SNPs set in global, previously unseen samples
priority <- data.frame(target = c("GOI", "NON_GOI"), priority = c(1, 2))
colnames(metadata)[which(colnames(metadata) == "Cefixime")] <- "target"
colnames(metadata)[which(colnames(metadata) == "Genome_Name")] <- "isolate"
metadata$target <- ifelse(metadata$target == "RESISTANT", "GOI", "NON_GOI")
AU_metadata <- metadata[metadata$Country == "Australia" | isolate == "Reference",]
VAL_metadata <- metadata[!isolate %in% AU_metadata$isolate,]
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
fwrite(all_results, "Cefixime_AU2global.csv", row.names = FALSE)