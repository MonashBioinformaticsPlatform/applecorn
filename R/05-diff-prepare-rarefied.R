
mk_main_rarefied_df <- function(ps_filt_rarefied, r_dara_dir, k = 1e6) {

  df_cnts <- ps_filt_rarefied %>%
                otu_table() %>%
                as.data.frame() %>%
                rownames_to_column(var = "sample") %>%
                as_tibble() %>%
                gather(seq, counts, -sample)

  df_taxa <- ps_filt_rarefied %>%
                tax_table() %>%
                as.data.frame() %>%
                rownames_to_column(var = 'seq') %>%
                as_tibble() %>%
                rowwise %>%
                mutate(md5 = substr(digest(seq, algo = 'md5', seed = 333),1,6)) %>%
                ungroup %>%
                #mutate(sp = ifelse(Species %>% is.na(), md5, Species)) %>%
                #unite(name, Phylum, Genus, md5, sep = "-", remove = F) %>%
                unite(name, Phylum, Class, Order, Family, Genus, md5, sep = ';') %>%
                #mutate(phylum = substr(Phylum, 1, 3)) %>%
                #unite(name, phylum, Family, Genus, md5, sep = "-") %>%
                dplyr::select(name, seq)

  meta_data <- ps_filt_rarefied %>%
                  sample_data() %>%
                  data.frame() %>%
                  rownames_to_column(var = 'sample') %>%
                  as_tibble()

  #aggregating cpm normalised counts with taxa
  # each row is unique taxa in a given sample
  df <- full_join(df_cnts, df_taxa, by = "seq")
  #here I'm adding metadata information to that data.frame
  df2 <- left_join(df, meta_data, by = "sample")
  # keeping df2 for legacy purposes, but this is main data.frame since
  # alot of the downstream analysis will rotate around that data.frame
  main_df <- df2

  main_df_fn <- paste0(r_data_dir, "/", "main_df_rarefied.RData")

  if(!file.exists(main_df_fn)) {
    save(main_df, file = main_df_fn)
  } else {
    load(main_df_fn)
  }

  return(main_df)
}
