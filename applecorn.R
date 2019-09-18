prog <- unlist(strsplit(commandArgs()[4], split = "--file="))[2]
base_prog <- basename(prog)
origin <- dirname(prog)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {

  config_fn <- "config.yml"
  config_full_fn <- file.path(origin, config_fn)
  file.copy(config_full_fn, config_fn)

  msg <- paste0("\n\n",
                "USAGE: ",
                base_prog,
                " <CONFIG>",
                "\n\n")

  stop(msg)
}

libs <- list.files(file.path(origin, "R"), full.names = TRUE)

config_fn <- normalizePath(args[1])
config <- yaml::read_yaml(config_fn)

debug <- config$debug

if(debug) {
  print(libs)
}
lapply(libs, source)

raw_data <- config$raw_data
r_data_dir <- config$r_data_dir
filt_path <- config$filt_path

report_fn <- config$report_fn
report_full_fn <- file.path(origin, report_fn)
file.copy(report_full_fn, report_fn)

report_out <- gsub(".Rmd", ".html", report_fn)
report_full_out <- file.path(origin, report_out)
file.copy(report_full_out, report_out)

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

opts <- list(origin = origin,
             dada = dada,
             taxtab = taxtab,
             ps = ps)

rmarkdown::render(input = report_fn,
                  params = opts)

