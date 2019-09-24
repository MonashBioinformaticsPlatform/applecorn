prog <- unlist(strsplit(commandArgs()[4], split = "--file="))[2]
base_prog <- basename(prog)
origin <- dirname(prog)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {

  config_fn <- "config.yml"
  config_full_fn <- file.path(origin, config_fn)

  config_fn <- gsub(".yml", "_example.yml", config_fn)

  if(!file.exists(config_fn)) {
    file.copy(config_full_fn, config_fn)
  }

  msg <- paste0("\n\n",
                "USAGE: ",
                base_prog,
                " <CONFIG>",
                "\n\n")

  stop(msg)
}

libs <- list.files(file.path(origin, "R"), pattern = ".R", full.names = TRUE)

config_fn <- normalizePath(args[1])
config <- yaml::read_yaml(config_fn)

debug <- config$debug

if(debug) {
  print(libs)
}
lapply(libs, source)

raw_data <- config$raw_data
r_data_dir <- config$r_data_dir
image_dir <- config$image_dir
filt_path <- config$filt_path

report_res <- config$report_fn
report_fn <- "report.Rmd"
report_full_fn <- file.path(origin, report_fn)
file.copy(report_full_fn, report_res)

multiple <- FALSE
if(length(raw_data) > 1) {
  multiple <- TRUE
}

seqtabs <- c()
dada <- NULL

if(multiple) {

  for(i in 1:length(raw_data)) {

    raw_dir <- basename(raw_data[[i]])
    print(raw_dir)

    fq_dir <- normalizePath(raw_data[[i]])
    filt_dir <- paste0(filt_path, "/", raw_dir)
    res_dir <- paste0(r_data_dir, "/", raw_dir)

    dada <- get_otus(fq_dir,
                     suffix = config$suffix,
                     split = config$split,
                     extn = config$extn,
                     trim_left = config$trim_left,
                     trunc_len = config$trunc_len,
                     #trunc_q = config$trunc_q,
                     max_ee = config$max_ee,
                     r_data_dir = res_dir,
                     filt_path = filt_dir)

    seqtabs <- c(seqtabs, dada[['seqtab_nochim']])

  }

}

seqtab_nochim <- NULL

if(!is.null(seqtab_nochim) & length(seqtabs) > 1) {
  seqtab_nochim <- mergeSequenceTables(seqtabs)
} else {
  dada <- get_otus(raw_data,
                   suffix = config$suffix,
                   split = config$split,
                   extn = config$extn,
                   trim_left = config$trim_left,
                   trunc_len = config$trunc_len,
                   max_ee = config$max_ee,
                   r_data_dir = r_data_dir,
                   filt_path = filt_path)

  seqtab_nochim <- dada[['seqtab_nochim']]

}

metadata <- normalizePath(config$metadata)
classifier <- normalizePath(config$taxa_db)

taxtab <- do_taxa_ann(seqtab_nochim,
                      classifier,
                      r_data_dir)

fitGTR <- fit_tree(seqtab_nochim,
                   taxatab,
                   r_data_dir = r_data_dir)

ps <- mk_phyloseq(taxtab,
                  seqtab_nochim,
                  fitGTR,
                  metadata,
                  r_data_dir = r_data_dir)

tree <- mk_tree(ps, plt_title = "Raw tree")
ps_filt <- filt_phyloseq(ps, tree[["mean"]], tree[["sd"]])
tree_filt <- mk_tree(ps_filt[["ps_filt"]], plt_title = "Filtered tree")

main_ps <- ps

if(config$ps_filt) {
  main_ps <- ps_filt[["ps_filt"]]
}

test_var <- config$test_var

plts <- list()

alpha <- mk_alpha(main_ps,
                  r_data_dir,
                  color = test_var)

alpha_fn <- paste0(image_dir, "/alpha.jpg")
plts[[alpha_fn]] <- alpha[["plot"]]

rare_curve <- NULL

if(config$rarefy) {
  rare_curve <- do_rare_curve(dada[["info"]], ps_filt[["ps_filt"]])
}

dist_binary <- config$dist_binary
dist_abund <- config$dist_abund

plt_binary <- mk_ordination(main_ps, dist_binary)
plt_binary_fn <- paste0(image_dir, "/", dist_binary, ".jpg")
plts[[plt_binary_fn]] <- plt_binary

plt_abund <- mk_ordination(main_ps, dist_abund)
plt_abund_fn <- paste0(image_dir, "/", dist_abund, ".jpg")
plts[[plt_abund_fn]] <- plt_abund

permanova <- permanova_pairwise(main_ps, test_var, dist_abund)

permanova_fn <- paste0(r_data_dir, "/", "permanova.csv")
permanova %>% write_tsv(permanova_fn)

main_df <- mk_main_df(main_ps, r_data_dir)

n_taxa <- config$n_taxa %>% as.numeric
legend <- config$legend

barplt <- mk_barplt(main_df,
                    r_data_dir,
                    n_taxa = n_taxa,
                    legend = legend,
                    test_var = test_var)

barplt_fn <- paste0(image_dir, "/", "barplot.jpg")
plts[[barplt_fn]] <- barplt

opts <- list(origin = origin,
             dada = dada,
             taxtab = taxtab,
             ps = list("raw" = ps, "filt" = main_ps),
	         tree = tree,
	         tree_filt = tree_filt,
	         alpha = alpha,
	         rare_curve = rare_curve,
             permanova = permanova,
             plt_binary = plt_binary,
             plt_abund = plt_abund,
             barplt = barplt)

if(config$mk_plots) {

  lapply(ord_plts %>% names, function(p) {ggsave(plot = ord_plts[[p]],
                                                 filename = p,
                                                 device = "jpg",
                                                 height = 4.4,
                                                 width = 6.05,
                                                 scale = 1,
                                                 dpi = 600,
                                                 units = "in")})
}

rmarkdown::render(input = report_res,
                  params = opts)

