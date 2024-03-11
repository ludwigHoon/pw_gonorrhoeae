# Process SNP matrix
library(minSNPs)
library(data.table)


metadata <- fread("no_dup_metadata.csv")

identify_conflict <- function(snp_matrix, dup_entry){
    library(BiocParallel)
    ml <- length(snp_matrix[[1]])
    position_to_exclude <- bplapply(1:nrow(dup_entry), function(entry) {
        matrix_s_names <- unlist(strsplit(dup_entry$files[entry], ", "))
        temp_exclude <- c()
        for (i in 1:ml){
            if (length(
                    unique(unlist(
                        generate_pattern(snp_matrix[matrix_s_names], i)
                    ))
                ) > 1){
                temp_exclude <- c(temp_exclude, i)
            }
        }
        return(temp_exclude)
    }, BPPARAM = MulticoreParam(workers = 120, progress = TRUE))
    names(position_to_exclude) <- dup_entry$Genome_Name
    return(position_to_exclude)
}



# Process AU matrix
AU_matrix <- read_fasta("snippy_output/AU.aln")
AU_SNP_ref <- fread("AU_SNP_POS.txt")
AU_dup_entry <- metadata[Country == "Australia" & grepl(", ", files)]
AU_conflict <- identify_conflict(AU_matrix, AU_dup_entry)
conflicts_positions <- data.table(cause = names(AU_conflict), positions = sapply(AU_conflict, paste0, collapse = ", "))
fwrite(conflicts_positions, "AU_conflict.csv", row.names = TRUE)
AU_conflict <- unique(unlist(AU_conflict))

No_conflict_pos <- 1:length(AU_matrix[[1]])
No_conflict_pos <- No_conflict_pos[!No_conflict_pos %in% AU_conflict]
AU_SNP_ref <- AU_SNP_ref[No_conflict_pos]
colnames(AU_SNP_ref) <- "Reference_Position"
AU_SNP_ref$Fasta_Position <- seq_len(nrow(AU_SNP_ref))
fwrite(AU_SNP_ref[,list(Fasta_Position, Reference_Position)], "FIN_AU_SNP_POS.csv", row.names = FALSE)

AU_final_matrix <- AU_matrix[names(AU_matrix) %in% c("Reference", metadata[Country == "Australia"]$Genome_Name)]
AU_final_matrix <- generate_pattern(AU_final_matrix, No_conflict_pos)

AU_final_matrix_no_dup_no_conflict <- strsplit(unlist(AU_final_matrix), split = "")
write_fasta(AU_final_matrix_no_dup_no_conflict, "AU_final_matrix.fasta")


# Process global matrix
global_matrix <- read_fasta("snippy_output/global.aln")
global_dup_entry <- metadata[grepl(", ", files)]
global_conflict <- identify_conflict(global_matrix, global_dup_entry)
conflicts_positions <- data.table(cause = names(global_conflict), positions = sapply(global_conflict, paste0, collapse = ", "))
fwrite(conflicts_positions, "global_conflict.csv", row.names = TRUE)
global_conflict <- unique(unlist(global_conflict))

No_conflict_pos <- 1:length(global_matrix[[1]])
No_conflict_pos <- No_conflict_pos[!No_conflict_pos %in% global_conflict]
global_SNP_ref <- fread("global_SNP_POS.txt")
global_SNP_ref <- global_SNP_ref[No_conflict_pos]
colnames(global_SNP_ref) <- "Reference_Position"
global_SNP_ref$Fasta_Position <- seq_len(nrow(global_SNP_ref))
fwrite(global_SNP_ref[,list(Fasta_Position, Reference_Position)], "FIN_global_SNP_POS.csv", row.names = FALSE)

global_final_matrix <- global_matrix[names(global_matrix) %in% c("Reference", metadata$Genome_Name)]
global_final_matrix <- generate_pattern(global_final_matrix, No_conflict_pos)

global_final_matrix_no_dup_no_conflict <- strsplit(unlist(global_final_matrix), split = "")
write_fasta(global_final_matrix, "global_final_matrix.fasta")
