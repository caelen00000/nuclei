library(scCustomize)
library(anndataR)
library(org.Dm.eg.db)
library(Seurat)
library(Matrix)

# change this to your own working_dir
setwd("C:\\Users\\camplain\\Documents\\snRNASeq_v2\\git_repo\\nuclei\\working_dir")

# Replaced alpha symbol with text, I think this is how the Fly Cell Atlas does it.
# TODO: replace with Flybase API calls?
manual_symbol_mapping <- data.frame(
  row.names = c("FBgn0289578", "FBgn0289579", "FBgn0289650", "FBgn0289678", "FBgn0289712", "FBgn0289725",    "FBgn0289743", "FBgn0289744", "FBgn0289745", "FBgn0289868", "FBgn0289933", "FBgn0289994", "FBgn0290109", "FBgn0013687"),
  symbol    = c("CG41994",     "Dld",         "jbug",        "pal",         "cpa",         "lncRNA:CR43126", "Cyp6a2",      "EndoG",       "Est17",       "Or67c",       "LanB1",       "Hsp27",       "eIF2alpha",   "mt:ori")
)

# For a given genotype and library list, this function converts the filtered .h5
# outputs from CellBender to .h5ad files for use in Python.
# Three .h5ad files are created per library, containing only the non-empty nuclei called by CellBender:
#   - One has the background counts removed
#   - One has the raw counts
#   - One has the raw counts with multimappers
# The ones with the raw counts may be useful if you don't trust CellBender's
# background removal or you want to see what was removed.
#
# At the same time, several useful things happen:
#   - FlyBase IDs are converted to gene symbols via the database in org.Dm.eg.db.
#     I've found it has missing symbols, so the mappings for these are added
#     from the manual_symbol_mapping dataframe. You'll get an error that tells 
#     you which IDs couldn't be mapped if they aren't in the dataframe either.
#     You'll have to look them up and add new entries.
#
#   - Per nucleus annotations are added, which will be useful later:
#     - library_id: e.g. "wildtype_a0_s1"
#     - genotype:   e.g. "wildtype"
#     - stage:      e.g. "a0"
#     - sample:     e.g. "s1"
#     - study:      always "PTTH". This is in distinction to "FCA", for when we
#                   integrate with the Fly Cell Atlas.
#
#     - {spliced, unspliced, ambiguous}_counts: number of reads with each splicing status
#     - fraction_unspliced_(unambiguous): fraction of unspliced reads vs either all reads or only unambiguous
#     - TODO: maybe include the entire Velocyto matrices
#
# AnnotationDbi::select and write_h5ad print some irrelevant suggestions so I suppressed them.
# If you have to debug this, consider deleting suppressMessages.
#
# Arguments:
#   cellbender_dir: the directory with your CellBender outputs
#   star_dir:       the STAR output directory
#   genotype:       either "wildtype" or "mutant"
#   library_list:   a list of the libraries you want to convert following my usual naming scheme but without the genotype, e.g. c("a0_s1", "a0_s2")
#   out_dir:        the directory where the new .h5ad files will be written


cellbender_h5_to_h5ad <- function(cellbender_dir,
                                  star_dir,
                                  genotype,
                                  library_list,
                                  out_dir){
  
  for(library in library_list){
    library_id <- paste0(genotype, "_", library)
    
    print(paste0("Processing ", library_id))
    
    stage = strsplit(library, "_")[[1]][1]
    sample = strsplit(library, "_")[[1]][2]
    
    matrix_cellbender <- Read_CellBender_h5_Mat(paste0(cellbender_dir, "\\", library_id, "\\cellbender_filtered.h5"))
    matrix_star       <- Read10X(paste0(star_dir, "\\", library_id, "\\Solo.out\\GeneFull_Ex50pAS\\raw"))
    matrix_star_multi <- Read10X(paste0(star_dir, "\\", library_id, "\\Solo.out\\GeneFull_Ex50pAS\\raw_multi_floor"))
    
    matrix_star       <- matrix_star[, colnames(matrix_star) %in% colnames(matrix_cellbender)]
    matrix_star_multi <- matrix_star_multi[, colnames(matrix_star_multi) %in% colnames(matrix_cellbender)]
    
    velocyto_dir <- paste0(star_dir, "\\", library_id, "\\Solo.out\\Velocyto\\raw\\")
    
    velocyto_barcodes <- read.delim(paste0(velocyto_dir, "barcodes.tsv"),
                                    header = FALSE)$V1
    
    spliced_matrix   <- readMM(paste0(velocyto_dir, "spliced.mtx"))
    unspliced_matrix <- readMM(paste0(velocyto_dir, "unspliced.mtx"))
    ambiguous_matrix <- readMM(paste0(velocyto_dir, "ambiguous.mtx"))
    
    colnames(spliced_matrix) <- velocyto_barcodes
    colnames(unspliced_matrix) <- velocyto_barcodes
    colnames(ambiguous_matrix) <- velocyto_barcodes
    
    spliced_matrix   <- spliced_matrix[, colnames(spliced_matrix) %in% colnames(matrix_cellbender)]
    unspliced_matrix <- unspliced_matrix[, colnames(unspliced_matrix) %in% colnames(matrix_cellbender)]
    ambiguous_matrix <- ambiguous_matrix[, colnames(ambiguous_matrix) %in% colnames(matrix_cellbender)]
    
    spliced_counts   <- colSums(spliced_matrix)
    unspliced_counts <- colSums(unspliced_matrix)
    ambiguous_counts <- colSums(ambiguous_matrix)
    
    # QClus includes ambiguous counts in the total, dropletQC doesn't. Why?
    total_counts_ambiguous   <- spliced_counts + unspliced_counts + ambiguous_counts
    total_counts_unambiguous <- spliced_counts + unspliced_counts
    
    fraction_unspliced_ambiguous <- unspliced_counts / total_counts_ambiguous
    fraction_unspliced_ambiguous[is.nan(fraction_unspliced_ambiguous)] <- 0
    
    fraction_unspliced_unambiguous <- unspliced_counts / total_counts_unambiguous
    fraction_unspliced_unambiguous[is.nan(fraction_unspliced_unambiguous)] <- 0
    
    write_anndata <- function(matrix, out_subdir){
      # transpose the matrix so obs and var aren't backwards
      anndata <- AnnData(t(matrix))
      
      suppressMessages(
        symbols <- AnnotationDbi::select(org.Dm.eg.db,
                                         keys = anndata$var_names,
                                         keytype = "FLYBASE",
                                         columns = c("SYMBOL"))$SYMBOL
      )
      
      unmapped_IDs <- anndata$var_names[is.na(symbols)]
      
      manual_mapped_symbols <- manual_symbol_mapping[unmapped_IDs, "symbol"]
      
      if(anyNA(manual_mapped_symbols)){
        print("Not all unmapped IDs could be manually mapped:")
        print(unmapped_IDs[is.na(manual_mapped_symbols)])
        
        stop()
      }
      
      symbols[is.na(symbols)] <- manual_mapped_symbols
      
      anndata$var_names <- symbols
      
      anndata$obs_names <- paste0(anndata$obs_names, "_", library_id)
      
      anndata$obs$library_id <- library_id
      anndata$obs$genotype   <- genotype
      anndata$obs$stage      <- stage
      anndata$obs$sample     <- sample
      anndata$obs$study      <- "PTTH"
      
      anndata$obs$spliced_counts   <- spliced_counts
      anndata$obs$unspliced_counts <- unspliced_counts
      anndata$obs$ambiguous_counts <- ambiguous_counts
      
      anndata$obs$fraction_unspliced_ambiguous <- fraction_unspliced_ambiguous
      anndata$obs$fraction_unspliced_unambiguous <- fraction_unspliced_unambiguous
      
      dir.create(paste0(out_dir, "\\", out_subdir),
                 recursive = TRUE,
                 showWarnings = FALSE)
      
      suppressMessages(
        anndata$write_h5ad(paste0(out_dir, "\\", out_subdir, "\\", library_id, ".h5ad"))
      )
    }
    
    write_anndata(matrix_cellbender, "background_removed")
    write_anndata(matrix_star, "calling_only")
    write_anndata(matrix_star_multi, "calling_only_multi")
  }
}

cellbender_h5_to_h5ad(cellbender_dir = "cellbender\\concat_all_multimapper_em_1st_pass_expect_10k_tuned",
                      star_dir       = "Z:\\Caelen\\snRNAseq_v2\\concat_all_multimapper_em_1st_pass",
                      genotype       = "wildtype",
                      library_list   = c("a0_s1", "a0_s2",
                                         "p1_s1", "p1_s2",
                                         "p2_s1", "p2_s2",
                                         "p3_s1", "p3_s2",
                                         "p4_s1", "p4_s2",
                                         "wp_s1", "wp_s2"),
                      out_dir        = "cellbender_h5ad\\expect_10k_tuned")




# Converts the lists of rRNA and mtRNA downloaded from FlyBase to lists of gene symbols
# No need to run this, this is just for reference.
flybase_gene_list_to_symbols <- function(){
  # Reads a list of rRNA genes obtained here: https://flybase.org/cgi-bin/cvreport.pl?cvterm=SO:0001637
  # I removed the bobbed genes, as well as the pre-rRNA.
  # I checked and unsurprisingly the pre-rRNA is expressed in almost every nucleus, so I figured it wouldn't be informative as a QC metric.
  rRNA <- read.delim("../preprocessing_and_cell_calling/flybase_rRNA_no_pre.txt", sep = "\n", header = FALSE)$V1
  
  rRNA <- AnnotationDbi::select(org.Dm.eg.db,
                                keys = rRNA,
                                keytype = "FLYBASE",
                                columns = c("SYMBOL"))$SYMBOL
  
  write.table(rRNA, "../preprocessing_and_cell_calling/symbol_rRNA_no_pre.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Reads a list of mtRNA genes obtained here: https://flybase.org/cgi-bin/cvreport.pl?cvterm=SO:0000088
  mtRNA <- read.delim("../preprocessing_and_cell_calling/flybase_mtRNA.txt", sep = "\n", header = FALSE)$V1
  
  mtRNA <- AnnotationDbi::select(org.Dm.eg.db,
                                 keys = mtRNA,
                                 keytype = "FLYBASE",
                                 columns = c("SYMBOL"))$SYMBOL
  
  write.table(mtRNA, "../preprocessing_and_cell_calling/symbol_mtRNA.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Reads a list of mitochondrial ribosomal protein coding genes obtained here (with mt:srRNA and mt:lrRNA excluded): https://flybase.org/reports/FBgg0000059.html
  # Note none of these are in the mitochondrial genome
  mito_ribo_protein <- read.delim("../preprocessing_and_cell_calling/flybase_mito_ribo_protein.txt", sep = "\n", header = FALSE)$V1
  
  mito_ribo_protein <- AnnotationDbi::select(org.Dm.eg.db,
                                             keys = mito_ribo_protein,
                                             keytype = "FLYBASE",
                                             columns = c("SYMBOL"))$SYMBOL
  
  write.table(mito_ribo_protein, "../preprocessing_and_cell_calling/symbol_mito_ribo_protein.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Reads a list of cytoplasmic ribosomal protein coding genes obtained here: https://flybase.org/reports/FBgg0000141.html
  ribo_protein <- read.delim("../preprocessing_and_cell_calling/flybase_ribo_protein.txt", sep = "\n", header = FALSE)$V1
  
  ribo_protein <- AnnotationDbi::select(org.Dm.eg.db,
                                        keys = ribo_protein,
                                        keytype = "FLYBASE",
                                        columns = c("SYMBOL"))$SYMBOL
  
  write.table(ribo_protein, "../preprocessing_and_cell_calling/symbol_ribo_protein.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
