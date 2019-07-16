
applecorn <- drake_plan(
  dada2_res = run_dada2(raw_data),
  taxtab = do_taxa_ann(dada2_res[["seqtab_nochim"]], classifier),
  tree = mk_tree(dada2_res[["seqtab_nochim"]], taxatab),
  ps = mk_phyloseq(taxtab, dada2_res[["seqtab_nochim"]], tree, metadata),
  ps_filt = filt_phyloseq(ps),
  alpha = mk_alpha(ps_filt[["ps_filt"]]),
  rare_curve = do_rare_curve(dada2_res[["info"]], alpha[["ps_filt_rarefied"]]),
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE)
)
