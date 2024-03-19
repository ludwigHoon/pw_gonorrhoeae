# Run mcc & percent for cefixime using each subset of global samples
library(BiocParallel)

# Multidrug - binary analysis
global_binary <- bplapply(c(1:5), function(parameters){
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

    ## Cefixime + Ciprofloxacin
    ## Ciprofloxacin + Penicillin
    GOIS <- list(
        Cefixime_Ciprofloxacin = metadata[Cefixime == "RESISTANT" & Ciprofloxacin == "RESISTANT"]$Genome_Name,
        Ciprofloxacin_Penicillin = metadata[Ciprofloxacin == "RESISTANT" & Penicillin == "RESISTANT"]$Genome_Name
    )
    for (nGOI in names(GOIS)){
        for (metric in metrics){
            GOI <- GOIS[[nGOI]]
            GOI <- GOI[GOI %in% names(global_matrix)]
            print(paste("Running", metric, "for", file, "at:", Sys.time()))
            result <- find_optimised_snps(seqc = global_matrix, metric = metric, goi = GOI,
                number_of_result = 10, max_depth = 5, bp = BP,
                output_progress = TRUE)
            result_name <- paste0(nGOI, "_", metric, "_global_", file)
            print(paste("Completed", metric, "for", file, "at:", Sys.time()))
            saveRDS(result, paste0(result_name,".rds"))
            try(output_result(result, view = "tsv", file_name = paste0(result_name, ".tsv"), seqc = global_matrix))
        }
    }    
    return(TRUE)
}, BPPARAM = BatchtoolsParam(workers = 5, cluster="slurm", template="~/slurm_bc_template.tmpl",
    resources=list(walltime=60*60*24*5, ncpus=21), progress = TRUE))


# Multidrug - multiclass analysis
global_multiclass <- bplapply(c(1:5), function(parameters){
    library(data.table)
    library(minSNPs)
    library(BiocParallel)

    BP <- MulticoreParam(workers = 80)
    file <- parameters
    metrics <- c("mcc_multi", "simpson_by_group")
    
    # Read in the global matrix
    global_matrix <- read_fasta(paste0("global_final_matrix_", file, ".fasta"))
    metadata <- fread("no_dup_metadata.csv")
    for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
        metadata[[drug]] <- ifelse(metadata[[drug]] == "INTERMEDIATE", "RESISTANT", metadata[[drug]])
    }
    metadata <- metadata[Genome_Name %in% names(global_matrix)]
    for (drug in c("Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")){
        metadata[[drug]] <- ifelse(metadata[[drug]] == "NOT_FOUND", "S", "R")
    }
    ## Cefixime + Ciprofloxacin
    ## Ciprofloxacin + Penicillin
    nGOIS <- c(
        #"Cefixime_Ciprofloxacin",
        "Ciprofloxacin_Penicillin")
    for (nGOI in nGOIS){
        mod_metadata <- metadata
        if (nGOIS == "Cefixime_Ciprofloxacin"){
            mod_metadata$target <- paste(mod_metadata$Cefixime, mod_metadata$Ciprofloxacin, sep = "_")
        } else {
            mod_metadata$target <- paste(mod_metadata$Ciprofloxacin, mod_metadata$Penicillin, sep = "_")
        }
        colnames(mod_metadata)[which(colnames(mod_metadata) == "Genome_Name")] <- "isolate"
        for (metric in metrics){
            print(paste("Running", metric, "for", file, "at:", Sys.time()))
            result <- find_optimised_snps(seqc = global_matrix, metric = metric,
                number_of_result = 10, max_depth = 5, bp = BP, meta = mod_metadata, target = "target",
                output_progress = TRUE)
            result_name <- paste0(nGOI, "_", metric, "_global_", file)
            print(paste("Completed", metric, "for", file, "at:", Sys.time()))
            saveRDS(result, paste0(result_name,".rds"))
            try(output_result(result, view = "tsv", file_name = paste0(result_name, ".tsv"),
                meta = mod_metadata, target = "target", seqc = global_matrix))
        }
    }    
    return(TRUE)
}, BPPARAM = BatchtoolsParam(workers = 5, cluster="slurm", template="~/slurm_bc_template.tmpl",
    resources=list(walltime=60*60*24*5, ncpus=21), progress = TRUE))