
permanova_pairwise <- function(ps, test_var, dist_type) {

  #'
  #' DO pairwise PERMANOVA tests.
  #'
  #' @param ps phyloseq object.
  #' @param test_var column name that has group names, e.g treat vs control
  #' @return results data.frame (R2, F.Model, p.value, FWER)
  #' @examples
  #' permanova <- permanova_pairwise(ps_filt, treat)
  #'

  samdf <- ps %>%
            sample_data() %>%
            as.data.frame()

  set.seed(15653)

  facts <- samdf[[test_var]] %>%
            unique %>%
            as.character()

  co <- combn(facts, 2)

  res <- list()

  for(i in 1:ncol(co)) {
    comp_name <- paste(co[,i], collapse = " vs ")
    #TODO do a check that test_var actually exists in the samdf
    samdf_pair <- samdf %>% dplyr::filter(get(test_var) %in% co[,i])

    filt <- samdf_pair$sample_id
    # this is fix in case sample ids are numbers. laters of d_pair[filt, filt] does subsample correctly
    filt <- filt %>% as.character

    d <- phyloseq::distance(ps, method = dist_type)
    d_pair <- d %>% as.matrix()
    d_pair <- d_pair[filt, filt] %>% as.dist

    f <- as.formula(paste("d_pair", test_var, sep = "~"))

    perm_pair <- adonis(f, data=samdf_pair)$aov.tab
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

  df_out <- df %>% arrange(fwer)

  return(df_out)
}
