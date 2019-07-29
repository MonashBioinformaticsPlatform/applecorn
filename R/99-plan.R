
dada2_plan <- drake_plan(
                  dada2_res = run_dada2(raw_data),
                  taxtab = do_taxa_ann(dada2_res[["seqtab_nochim"]], classifier),
                  tree = mk_tree(dada2_res[["seqtab_nochim"]], taxatab)
              )

kraken_filt_plan <- drake_plan(
                        otus_fasta = mk_otus_fasta(dada2_res[["seqtab_nochim"]], taxtab),
                        kraken_filt = do_taxtab_filt(config$kraken_contams_fn, otus_fasta[["df_seq"]]),
                        kraken_filt2 = do_seqtab_filt(config$kraken_contams_fn, df_seq, dada2_res[["seqtab_nochim"]])
                      )

ps_plan <- drake_plan(
               ps = mk_phyloseq(taxtab, dada2_res[["seqtab_nochim"]], tree, metadata),
               ps_filt = filt_phyloseq(ps),
               alpha = mk_alpha(ps_filt[["ps_filt"]])
             )

rare_plan <- drake_plan(
                 rare_curve = do_rare_curve(dada2_res[["info"]], ps_filt[["ps_filt"]])
             )

#wunifrac = mk_ordination(ps_filt[['ps_filt']], "wunifrac")
#unifrac = mk_ordination(ps_filt[['ps_filt']], dist = "unifrac")

if(config$kraken_filt %>% as.logical) {
  applecorn = bind_plans(dada2_plan,
                         kraken_filt_plan,
                         ps_plan)
} else {
  applecorn = bind_plans(dada2_plan,
                         ps_plan)
}


#applecorn = drake_plan(
#  steps,
#  report = rmarkdown::render(
#    knitr_in("report.Rmd"),
#    output_file = file_out("report.html"),
#    quiet = TRUE)
#)
