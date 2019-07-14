
# ---- metadata ---

metadata_fn <- "../mapping.tsv"
r_data_dir <- "../data"

taxtab_fn <- paste0(r_data_dir, "/", "taxtab.rds")
seqtab_nochim_fn <- paste0(r_data_dir, "/", "seqtab_nochim.rds")
fitGTR_fn <- paste0(r_data_dir, "/", "fitGTR.rds")
meta_fn <- paste0(r_data_dir, "/", "metadata.rds")

seqtab_nochim <- readRDS(seqtab_nochim_fn)
taxtab <- readRDS(taxtab_fn)
fitGTR <- readRDS(fitGTR_fn)

# Here we have two sources of metadata "typical" source with several factors related
# to the experiment, like treatment and other attributes related to the study like gender
# And the second source is food frequecy questionaries (FFQ), which should also be considered
# From memory there were some samples that had done FFQs but haven't been sequenced. Basically I
# need to account for that and filter both metadata files after merging. For downstream analysis
# phyloseq object needs to be "balanced"

meta_data <- meta_data %>% select(-BarcodeSequence,
                                  -LinkerPrimerSequence,
                                  -Description,
                                  -Sample_Well)

samdf <- meta_data %>%
          column_to_rownames(var = 'sample') %>%
          as.data.frame()

ps <- phyloseq(tax_table(taxtab),
               sample_data(samdf),
               otu_table(seqtab_nochim, taxa_are_rows = FALSE),
               phy_tree(fitGTR$tree))

ps

ps_fn_rds <- paste0(r_data_dir, "/", "ps.rds")
ps_fn_rdata <- paste0(r_data_dir, "/", "ps.RData")

if(!file.exists(ps_fn_rdata)) {
  saveRDS(object = ps, file = ps_fn_rds)
  save(ps, file = ps_fn_rdata)
} else {
  load(ps_fn_rdata)
}

if(!file_test("-d", "../images")) {
  dir.create("../images", recursive = T)
}
