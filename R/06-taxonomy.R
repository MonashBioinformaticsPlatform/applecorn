
# ---- taxonomy ----

mk_barplt <- function(main_df,
                      r_data_dir,
                      sample_order = NULL,
                      n_taxa = 12,
                      legend = "right",
                      test_var = "treat") {

  # legend = "none" to turn it off
  # n_taxa = NULL to turn it off

  if(is.null(n_taxa)) {
    n_taxa <- main_df %>%
               select(name) %>%
               distinct %>%
               unlist %>%
               length
  }

  top_taxa <- main_df %>%
                group_by(name) %>%
                summarise(n = sum(counts)) %>%
                arrange(-n) %>%
                head(n = n_taxa) %>%
                dplyr::select(name) %>%
                unlist %>%
                unname

  df_filt <- NULL

  if(!is.null(sample_order)) {
    df_filt <- main_df %>% mutate(sample = factor(sample, level = sample_order))
  } else {
    df_filt <- main_df
  }

  df_filt <- df_filt %>%
               dplyr::filter(name %in% top_taxa) %>%
               separate(name, c('phylum', 'class', 'order', 'family', 'genus', 'md5'), sep = ";") %>%
               unite(name, c('phylum', 'class', 'order', 'family', 'genus', 'md5'), sep = "\n")

  f <- as.formula(paste("~", test_var))

  p_taxa <- df_filt %>%
              ggplot(aes(sample, cpm, fill = name)) +
                geom_bar(stat = "identity") +
                facet_wrap(f, scales = "free_x") +
                ggtitle(paste0(n_taxa, " most abundant taxa")) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position = legend,
                      axis.title.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.title.y=element_blank())

  if(n_taxa <= 12) {
    p_taxa <- p_taxa + scale_fill_brewer(palette="Set3")
  }

  return(p_taxa)

}
