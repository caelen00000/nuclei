library(scCustomize)
library(anndataR)
library(org.Dm.eg.db)

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
# AnnotationDbi::select and write_h5ad print some irrelevant suggestions so I suppressed them.
# If you have to debug this, consider deleting suppressMessages.

# Arguments:
#   cellbender_dir: the directory with your CellBender outputs
#   genotype:       either "wildtype" or "mutant"
#   library_list:   a list of the libraries you want to convert following my usual naming scheme but without the genotype, e.g. c("a0_s1", "a0_s2")
#   out_dir:        the directory where the new .h5ad files will be written
cellbender_h5_to_h5ad <- function(cellbender_dir,
                                  genotype,
                                  library_list,
                                  out_dir){
  dir.create(out_dir, showWarnings = FALSE)
  
  for(library in library_list){
    library_id <- paste0(genotype, "_", library)
    stage = strsplit(library, "_")[[1]][1]
    sample = strsplit(library, "_")[[1]][2]
    
    matrix <- Read_CellBender_h5_Mat(paste0(cellbender_dir, "\\", library_id, "\\cellbender_filtered.h5"))
    
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
    
    suppressMessages(
      anndata$write_h5ad(paste0(out_dir, "\\", library_id, ".h5ad"))
    )
  }
}

cellbender_h5_to_h5ad(cellbender_dir = "cellbender\\concat_all_multimapper_em_1st_pass_auto",
                      genotype       = "wildtype",
                      library_list   = c("a0_s1", "a0_s2",
                                         "p1_s1", "p1_s2",
                                         "p2_s1", "p2_s2",
                                         "p3_s1", "p3_s2",
                                         "p4_s1", "p4_s2",
                                         "wp_s1", "wp_s2"),
                      out_dir        = "raw_h5ad")



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
}