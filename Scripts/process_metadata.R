library(data.table)
ng_paarsnp <- fread("Neisseria_gonorrhoeae__paarsnp.csv") # ng_paarsnp is metadata output from a pipeline predicting the resistance/susceptibility of Neisseria gonorrhoeae to antibiotics
ng_mlst <- fread("Neisseria_gonorrhoeae__mlst.csv") # ng_mlst is metadata output from a genotyping pipeline outputting the MLST of Neisseria gonorrhoeae
ng_meta <- fread("Neisseria_gonorrhoeae__metadata.csv") # ng_meta is Neisseria gonorrhoeae metadata input by user 

ng_merged_meta <- merge(ng_paarsnp, ng_mlst, by = "Genome ID")
ng_merged_meta <- merge(ng_merged_meta, ng_meta, by.x = "Genome ID", by.y = "id")

stopifnot(nrow(ng_merged_meta) == nrow(ng_paarsnp))
stopifnot(nrow(ng_merged_meta) == nrow(ng_mlst))
stopifnot(nrow(ng_merged_meta) == nrow(ng_meta))

colnames(ng_merged_meta)[which(colnames(ng_merged_meta) == "Genome Name.x")] <- "Genome_Name"

# Handling country metadata
ng_merged_meta[Country == "" & country != ""]$Country <- ng_merged_meta[Country == "" & country != ""]$country
# Fixing misnomers
print(sort(unique(ng_merged_meta$Country)))
ng_merged_meta[Country == "Brasil"]$Country <- "Brazil"
ng_merged_meta[Country == "HongKong"]$Country <- "Hong Kong"
ng_merged_meta[Country == "The Netherlands"]$Country <- "Netherlands"
ng_merged_meta[Country == "UK"]$Country <- "United Kingdom"
ng_merged_meta[Country == "USA"]$Country <- "United States"
ng_merged_meta[Country == "Viet Nam"]$Country <- "Vietnam"

umerged <- ng_merged_meta[,c("Genome_Name", "Country", "ST", "abcZ", "adk", "aroE", "fumC", "gdh", "pdhC", "pgm", "Cefixime", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline")]
umerged[is.na(as.numeric(umerged$ST))]$ST <- ""
SAMPLES_W_MULTIPLE_ENTRIES <- names(which(table(umerged$Genome_Name) > 1))

umerged$to_ignored <- FALSE
umerged[Genome_Name %in% SAMPLES_W_MULTIPLE_ENTRIES]$to_ignored <- TRUE

pb <- txtProgressBar(min = 0, max = length(SAMPLES_W_MULTIPLE_ENTRIES), style = 3)
# Remove those with conflicting MLSTs/antibiotic resistance predictions
for (sample_n in SAMPLES_W_MULTIPLE_ENTRIES){
    setTxtProgressBar(pb, pb$getVal() + 1)
    # Check if all the metadata are the same
    all_meta_same <- nrow(unique(umerged[Genome_Name == sample_n])) == 1
    # if not, remove the sample
    if (!all_meta_same){
        next
    } else {
        # if all metadata is the same, keep the first one
        umerged[Genome_Name == sample_n][1]$to_ignored <- FALSE
    }
}
close(pb)

no_dup <- umerged[to_ignored == FALSE]
print(paste("Final number of samples with no conflicts:", nrow(no_dup)))

# Merge with available FASTA files
library(BiocParallel)

ALL_FASTA_FILES <- gsub(".fasta", "", list.files("o_fastas", pattern = ".fasta"))
res <- bplapply(no_dup$Genome_Name, function(x){
    q <- paste0(x, "(_[0-9])?$")
    paste(grep(q, ALL_FASTA_FILES, value = TRUE), collapse = ", ")
}, BPPARAM = MulticoreParam(workers = 60, progress = TRUE))

no_dup$files <- unlist(res)


# Assemble Australia only SNP matrix
samples <- gsub(", ", " ", paste0(no_dup[Country == "Australia"]$files, collapse = " "))
command <- paste0("snippy-core --ref=reference_TUM19854.fasta --prefix=", "AU ", samples)
writeLines(command, "snippy-core_AU.sh")

# Assemble global SNP matrix
samples <- gsub(", ", " ", paste0(no_dup$files, collapse = " "))
command <- paste0("snippy-core --ref=reference_TUM19854.fasta --prefix=", "global ", samples)
writeLines(command, "snippy-core_global.sh")


no_dup <- rbind(no_dup, data.table(`Genome_Name` = "Reference", Country = "",
    ST = "7363", abcZ = "59", adk = "39", aroE = "67", fumC = "78", gdh = "148", pdhC = "153", pgm = "65",
    Cefixime = "RESISTANT", Azithromycin = "NOT_FOUND", Ceftriaxone = "RESISTANT",
    Ciprofloxacin = "RESISTANT", Penicillin = "RESISTANT", Sulfonamides = "RESISTANT",
    Spectinomycin = "NOT_FOUND", Tetracycline = "RESISTANT", to_ignored = FALSE, files = "reference_TUM19854"))

fwrite(no_dup, "no_dup_metadata.csv", row.names = FALSE)