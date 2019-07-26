
# ---- metadata ---

mk_phyloseq <- function(taxtab, seqtab_nochim, fitGTR, meta_data) {

  samdf <- read_tsv(meta_data) %>%
  #samdf <- read_csv(meta_data) %>%
            mutate(sample_id = sample) %>%
            column_to_rownames(var = 'sample') %>%
            as.data.frame()

  ps <- phyloseq(tax_table(taxtab),
                 sample_data(samdf),
                 otu_table(seqtab_nochim, taxa_are_rows = FALSE),
                 phy_tree(fitGTR$tree))

  r_data_dir <- "data"
  ps_fn_rds <- paste0(r_data_dir, "/", "ps.rds")
  ps_fn_rdata <- paste0(r_data_dir, "/", "ps.RData")

  if(!file.exists(ps_fn_rdata)) {
    saveRDS(object = ps, file = ps_fn_rds)
    save(ps, file = ps_fn_rdata)
  } else {
    load(ps_fn_rdata)
  }

  if(!file_test("-d", "images")) {
    dir.create("images", recursive = T)
  }

  return(ps)
}

# ---- phylo_tree ----

filt_phyloseq <- function(ps) {

  p_tree1 <- plot_tree(ps, method = "treeonly", title = "Raw tree")

  ps_tree <- phy_tree(ps)
  m <- ps_tree$edge.length %>% mean
  s <- ps_tree$edge.length %>% sd

  ps_filt <- prune_taxa(ps_tree$tip.label[ps_tree$edge[ps_tree$edge.length < m+(3*s), 2]], ps)
  p_tree2 <- plot_tree(ps_filt, method = "treeonly", title = "Filtred tree")

  outlier_seqs <-ps_tree$tip.label[ps_tree$edge[ps_tree$edge.length > m+(3*s), 2]]

  r_data_dir <- "data"
  ps_filt_fn <- paste0(r_data_dir, "/", "ps_filt.RData")

  if(!file.exists(ps_filt_fn)) {
    save(ps_filt, file = ps_filt_fn)
  } else {
    load(ps_filt_fn)
  }

  return(list("raw_tree" = p_tree1,
	      "filt_tree" = p_tree2,
	      "outlier_seqs" = outlier_seqs,
	      "ps_filt" = ps_filt))
}
