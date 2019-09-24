
# ---- alpha1 ----

mk_alpha <- function(ps_filt, r_data_dir, color = NULL) {

  rarefy_depth <- ps_filt %>% sample_sums %>% min

  ps_filt_rarefied <- rarefy_even_depth(ps_filt,
                                        sample.size = rarefy_depth,
                                        replace = FALSE,
                                        rngseed = 13)

  ps_filt_rarefied <- prune_taxa(taxa_sums(ps_filt_rarefied) > 0, ps_filt_rarefied)

  ps_filt_rarefied_fn <- paste0(r_data_dir, "/", "ps_filt_rarefied.RData")

  if(!file.exists(ps_filt_rarefied_fn)) {
    save(ps_filt_rarefied, file = ps_filt_rarefied_fn)
  } else {
    load(ps_filt_rarefied_fn)
  }

  p_alpha1 <- plot_richness(ps_filt_rarefied,
                            x="sample",
                            #measures=c("Observed", "Chao1", "Shannon", "Simpson"),
                            measures=c("Observed", "Chao1", "Shannon"),
                            color=color)

  return(list("ps_filt_rarefied" = ps_filt_rarefied,
              "plot" = p_alpha1,
              "min_depth" = rarefy_depth))
}
