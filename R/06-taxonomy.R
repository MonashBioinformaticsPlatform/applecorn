
# ---- taxonomy ----

unknown_yet <- function() {
  r_data_dir <- "../data"
  main_df_fn <- paste0(r_data_dir, "/", "main_df.RData")
  load(file = main_df_fn)

  n_taxa <- 12

  top_taxa <- main_df %>%
                group_by(name) %>%
                summarise(n = sum(counts)) %>%
                arrange(-n) %>%
                head(n = n_taxa) %>%
                dplyr::select(name) %>%
                unlist %>%
                unname

  p_taxa <- main_df %>%
              dplyr::filter(name %in% top_taxa) %>%
              unite(sample2, sample, cage, sep = "-") %>%
              #ggplot(aes(sample, counts, fill = name)) +
              ggplot(aes(sample2, counts, fill = name)) +
                geom_bar(stat = "identity") +
                scale_fill_brewer(palette="Set3") +
                facet_wrap(~treat, scales = "free_x") +
                ggtitle(paste0(n_taxa, " most abundant taxa")) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))


  p_taxa
}
