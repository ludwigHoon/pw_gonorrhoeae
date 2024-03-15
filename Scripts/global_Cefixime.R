# Run mcc & percent for cefixime using each subset of global samples
library(BiocParallel)

successfulRUN <- bplapply(c(1:5), function(parameters){
    library(data.table)
    library(minSNPs)
    library(BiocParallel)

    BP <- MulticoreParam(workers = 80)
    file <- parameters
    metrics <- c("mcc", "percent")
    
    # Read in the global matrix
    global_matrix <- read_fasta(paste0("global_final_matrix_", file, ".fasta"))
    metadata <- fread("no_dup_metadata.csv")
    for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
        metadata[[drug]] <- ifelse(metadata[[drug]] == "INTERMEDIATE", "RESISTANT", metadata[[drug]])
    }
    GOI <- metadata[Cefixime == "RESISTANT"]$Genome_Name
    GOI <- GOI[GOI %in% names(global_matrix)]

    for (metric in metrics){
        print(paste("Running", metric, "for", file, "at:", Sys.time()))
        result <- find_optimised_snps(seqc = global_matrix, metric = metric, goi = GOI,
            number_of_result = 10, max_depth = 5, bp = BP,
            output_progress = TRUE)
        result_name <- paste0("Cefixime_", metric, "_global_", file)
        print(paste("Completed", metric, "for", file, "at:", Sys.time()))
        saveRDS(result, paste0(result_name,".rds"))
        try(output_result(result, view = "tsv", file_name = paste0(result_name, ".tsv"), seqc = global_matrix))
    }
    
    return(TRUE)
}, BPPARAM = BatchtoolsParam(workers = 5, cluster="slurm", template="~/slurm_bc_template.tmpl",
    resources=list(walltime=60*60*24*5, ncpus=32), progress = TRUE))
