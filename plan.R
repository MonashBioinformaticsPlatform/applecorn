
if(debug) {
  print(config)
}

kraken <- config$kraken_filt %>% as.logical
kraken_contams_fn <- config$kraken_contams_fn

rarefy <- config$rarefy %>% as.logical

applecorn <- drake_plan(
                        dada2_res = run_dada2(raw_data),
                        taxtab = do_taxa_ann(dada2_res[["seqtab_nochim"]], classifier),
                        tree = mk_tree(dada2_res[["seqtab_nochim"]], taxatab),

                        rare_curve = target(command = do_rare_curve(dada2_res[["info"]], ps_filt[["ps_filt"]]),
                                            trigger = trigger(condition = rarefy)),

                        otus_fasta = target(command = mk_otus_fasta(dada2_res[["seqtab_nochim"]], taxtab),
                                            trigger = trigger(condition = kraken)),
                        kraken_filt = target(command = do_taxtab_filt(kraken_contams_fn, otus_fasta[["df_seq"]]),
                                             trigger = trigger(condition = kraken)),
                        kraken_filt2 = target(command = do_seqtab_filt(kraken_contams_fn, otus_fasta[["df_seq"]], dada2_res[["seqtab_nochim"]]),
                                              trigger = trigger(condition = kraken)),

                        ps = mk_phyloseq(taxtab, dada2_res[["seqtab_nochim"]], tree, metadata),
                        ps_filt = filt_phyloseq(ps),
                        alpha = mk_alpha(ps_filt[["ps_filt"]]),

                        wunifrac = mk_ordination(ps_filt[['ps_filt']], "wunifrac"),
                        unifrac = mk_ordination(ps_filt[['ps_filt']], dist = "unifrac"),

                        report = rmarkdown::render(
                          knitr_in("report.Rmd"),
                          output_file = file_out("report.html"),
                          quiet = TRUE)

                        )
                        #kraken = target(
                        #          command = (
                        #                     otus_fasta = mk_otus_fasta(dada2_res[["seqtab_nochim"]], taxtab),
                        #                     kraken_filt = do_taxtab_filt(kraken_contams_fn, otus_fasta[["df_seq"]]),
                        #                     kraken_filt2 = do_seqtab_filt(kraken_contams_fn, otus_fasta[["df_seq"]], dada2_res[["seqtab_nochim"]])
                        #                     ),
                        #          trigger = trigger(condition = kraken)
                        #          )
                        #)

