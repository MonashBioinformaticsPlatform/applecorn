
# ---- metadata ---

mk_phyloseq <- function(taxtab,
                        seqtab_nochim,
                        fitGTR,
                        meta_data,
                        r_data_dir = "data") {

  samdf <- read_tsv(meta_data) %>%
  #samdf <- read_csv(meta_data) %>%
            mutate(sample_id = sample) %>%
            column_to_rownames(var = 'sample') %>%
            as.data.frame()

  ps <- phyloseq(tax_table(taxtab),
                 sample_data(samdf),
                 otu_table(seqtab_nochim, taxa_are_rows = FALSE),
                 phy_tree(fitGTR$tree))

  ps_fn_rds <- paste0(r_data_dir, "/", "ps.rds")
  ps_fn_rdata <- paste0(r_data_dir, "/", "ps.RData")

  if(!file.exists(ps_fn_rdata)) {
    saveRDS(object = ps, file = ps_fn_rds)
    save(ps, file = ps_fn_rdata)
  } else {
    load(ps_fn_rdata)
    ps <- readRDS(ps_fn_rds)
  }

  return(ps)
}

# ---- phylo_tree ----

mk_tree <- function(ps, plt_title = "") {

  p_tree <- plot_tree(ps, method = "treeonly", title = plt_title)

  ps_tree <- phy_tree(ps)
  m <- ps_tree$edge.length %>% mean
  s <- ps_tree$edge.length %>% sd

  return(list("tree" = p_tree,
              "mean" = m,
              "sd" = s))
}

filt_phyloseq <- function(ps, m, s) {

  ps_tree <- phy_tree(ps)
  ps_filt <- prune_taxa(ps_tree$tip.label[ps_tree$edge[ps_tree$edge.length < m+(3*s), 2]], ps)

  outlier_seqs <- ps_tree$tip.label[ps_tree$edge[ps_tree$edge.length > m+(3*s), 2]]

  return(list("outlier_seqs" = outlier_seqs,
              "ps_filt" = ps_filt))
}
