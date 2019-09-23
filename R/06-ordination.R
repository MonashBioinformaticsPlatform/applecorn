# ---- ordination_init ----

mk_ordination <- function(ps_filt, dist = "wunifrac") {

  #options :
  # - wunifrac
  # - unifrac
  # - jacard
  # - bray
  pslog_filt <- transform_sample_counts(ps_filt, function(x) log(1 + x))

  # ---- do_filt ----

  #Note need to think here about the data.
  # the transform_sample_counts() function under the hood simply
  # runs apply with a given function on otu_table. It feels a bit simplified since we are
  # not normalising for the library sizes here. I think distance plots PCA/PCoA
  # would be more informative when libs had been normalised for depth..?

  # this can transform otu_table into workable form, depending on what you want todo
  # you might not even need to spread there.

  ord_log <- ordinate(pslog_filt, method = "PCoA", distance = dist)
  plt <- plot_ordination(pslog_filt,
                         ord_log,
                         type = "samples",
                         color = "treat",
                         #shape = "cage",
                         axes=1:2) +
                          geom_point(size = 4) +
                          #geom_text(aes(label = sample_id),
                          #          hjust = 1.3,
                          #          angle = -15,
                          #          size = 3) +
                          ggtitle(dist)
                          #scale_color_manual(values = c("#94641F", "#C5944E", "#008000", "#7CFC00"))# +

  return(list("ord_plt" = plt))

}

#mk_plts <- FALSE
#
#if(mk_plts) {
#  lapply(ord_plts %>% names, function(p) {ggsave(plot = ord_plts[[p]],
#                                                 filename = p,
#                                                 device = "jpg",
#                                                 height = 4.4,
#                                                 width = 6.05,
#                                                 scale = 1,
#                                                 dpi = 600,
#                                                 units = "in")})
#}

# ---- permanova_pairwise ----

#'
#' DO pairwise PERMANOVA tests.
#'
#' @param ps phyloseq object.
#' @param test_var column name that has group names, e.g treat vs control
#' @return results data.frame (R2, F.Model, p.value, FWER)
#' @examples
#' permanova <- permanova_pairwise(ps_filt, treat)
#'
permanova_pairwise <- function(ps) {

  samdf <- ps %>% sample_data

  set.seed(15653)

  facts <- samdf$treat %>% unique
  co <- combn(unique(as.character(facts)),2)
  res <- list()

  for(i in 1:ncol(co)) {
    comp_name <- paste(co[,i], collapse = " vs ")
    #TODO do a check that test_var actually exists in the samdf

    samdf_pair <- samdf %>%
                    dplyr::filter(treat %in% co[,i])

    filt <- samdf_pair %>%
              dplyr::select(sample_id) %>%
              unlist

    #d <- phyloseq::distance(ps, method = "bray")
    d <- phyloseq::distance(ps, method = "wunifrac")
    d_pair <- d %>% as.matrix()
    d_pair <- d_pair[filt, filt] %>% as.dist

    perm_pair <- adonis(d_pair ~ treat, data=samdf_pair)$aov.tab
    # the order: Name, F.Model, R2, p.value
    df <- data.frame("R2" = perm_pair[1,5],
                     "F.Model" = perm_pair[1,4],
                     "p.value" = perm_pair[1,6])

    res[[comp_name]] <- df
  }

  df <- do.call(rbind, res) %>%
          rownames_to_column(var = "comps") %>%
          as_tibble

  fwer <- p.adjust(df$p.value, method = 'bonferroni')
  df[['fwer']] <- fwer

  return(df %>% arrange(fwer))
}


#permanova <- permanova_pairwise(ps_filt_rarefied)
#permanova
#
#permanova_fn <- paste0(r_data_dir, "/", "permanova.csv")
#
#if(!file.exists(permanova_fn)) {
#  permanova %>% write_tsv(permanova_fn)
#}
