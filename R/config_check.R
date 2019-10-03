# TODO
# want to check test_var variable and some other I guess as well
# to somehow make sure user have looked at them and changed

chk_fn <- function(fn) {

  if(!file.exists(fn)) {
    msg <- paste0("MSG: File doesn't exists ", fn, "\n")
    cat(msg)
    stop()
  }

}

conf_chk <- function(config) {

  dist_binary <- config$dist_binary
  dist_abund <- config$dist_abund
  fit_tree <- config$fit_tree

  if(!fit_tree & (dist_binary == "unifrac") ) {
    msg <- paste0("MSG:",
                  "\n",
                  "Can't use UniFrac distance without phylogenetic tree.",
                  "\n",
                  "Check fit_tree option in the config OR change distance e.g jaccard",
                  "\n")
    cat(msg)
    stop()
  }

  if(!fit_tree & (dist_abund == "wunifrac") ) {
    msg <- paste0("MSG:",
                  "\n",
                  "Can't use weighted UniFrac distance without phylogenetic tree.",
                  "\n",
                  "Check fit_tree option in the config OR change distance e.g bray",
                  "\n")
    cat(msg)
  }

  raw_data <- config$raw_data
  r_data_dir <- config$r_data_dir
  images_dir <- config$images_dir

  metadata <- normalizePath(config$metadata)
  classifier <- normalizePath(config$taxa_db)

  chks <- c(raw_data, metadata, classifier)
  lapply(chks, chk_fn)

}
