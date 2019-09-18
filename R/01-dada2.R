# ---- get_seqtab ----

get_otus <- function(fq_dir,
                     suffix = "_R1,_R2",
                     split = "_R",
                     extn = "fastq.gz",
                     trim_left = c(23, 23),
                     trunc_len = c(245, 200),
                     max_ee = c(2, 2),
                     trunc_q = 2,
                     r_data_dir = "data",
                     filt_path = "filtered",
                     large = FALSE,
                     multithread = TRUE,
                     omega_a = 1e-40,
                     omega_c = 1e-40) {

  #' Denoising fastq and OTU/ASV picking
  #'
  #' @param fq_dir path to fastq directory
  #' @param suffix
  #' @param split substring to split on in order to get sample name
  #' @param extn file extension to look for in fq_dir
  #' @param trim_left is a filterAndTrim() function paramenter, for more info ?dada2::filterAndTrim
  #' @param trunc_len is a filterAndTrim() function paramenter, for more info ?dada2::filterAndTrim
  #' @param r_data_dir directory to store RData and .rds files as well as other resulting files
  #' @param filt_path directory for filtered fastq files
  #' @param threads [4] number of threads to use for parallel processing
  #' @return named list with two items otus = seqtab_nochim and info = df_info
  #' @examples
  #' add(1, 1)
  chk <- strsplit(suffix, split = ",") %>% unlist

  paired = FALSE
  r1 <- NULL
  r2 <- NULL

  if(length(chk) == 2) {
    paired = TRUE
    r2 <- chk[2]
  }
  else if(length(chk) > 2) {
    #TODO add more information to the error message
    stop("ERROR: Wrong suffix value")
  }
  else if(length(chk) <= 0) {
    stop("ERROR: Wrong suffix value")
  }

  r1 <- chk[1]

  if(!file_test("-d", r_data_dir)) {
    dir.create(r_data_dir, recursive = T)
  }

  if(!file_test("-d", filt_path)) {
    dir.create(filt_path, recursive = T)
  }

  fns <- sort(list.files(fq_dir, pattern = extn, full.names = TRUE))

  fnFs <- NULL
  fnRs <- NULL
  filtFs <- NULL
  filtRs <- NULL

  if(paired) {

    fnFs <- fns[grepl(r1, fns)]
    fnRs <- fns[grepl(r2, fns)]

    filtFs <- file.path(filt_path, basename(fnFs))
    filtRs <- file.path(filt_path, basename(fnRs))

  } else {
    fnFs <- fns[grepl(r1, fns)]
    filtFs <- file.path(filt_path, basename(fnFs))
  }

  sample_names <- sapply(strsplit(basename(filtFs), split), `[`, 1)

  filt_n_trim <- list(compress = TRUE,
                      trimLeft = trim_left,
                      truncLen = trunc_len,
                      truncQ = trunc_q,
                      trimRight = 0,
                      maxLen = Inf,
                      minLen = 20,
                      maxN = 0,
                      minQ = 0,
                      maxEE = max_ee,
                      rm.phix = TRUE,
                      orient.fwd = NULL,
                      matchIDs = FALSE,
                      id.sep = "\\s",
                      id.field = NULL,
                      n = 1e+05,
                      OMP = TRUE,
                      verbose = TRUE)

  trimed_dada2_fn = paste0(r_data_dir, "/", "trimed_dada2.rds")
  trimed_dada2 <- NULL

  if(!file.exists(trimed_dada2_fn)) {

    threads <- 1

    if(multithread) {
      threads <- detectCores()
    }

    if(paired) {
      trimed_dada2 <- mcmapply(fastqPairedFilter,
                               mapply(c,
                                      fnFs,
                                      fnRs,
                                      SIMPLIFY = FALSE),
                               mapply(c,
                                      filtFs,
                                      filtRs,
                                      SIMPLIFY = FALSE),
                               MoreArgs = filt_n_trim,
                               mc.cores = threads,
                               mc.silent = TRUE)
    }
    else {
      trimed_dada2 <- mcmapply(fastqFilter,
                               fnFs,
                               filtFs,
                               MoreArgs = filt_n_trim,
                               mc.cores = threads,
                               mc.silent = TRUE)
    }

    if(is.null(trimed_dada2)) {
      stop("ERROR: Can't happen")
    }

    trimed_dada2 <- trimed_dada2 %>% t
    rownames(trimed_dada2) <- sample_names

    saveRDS(object = trimed_dada2, file = trimed_dada2_fn)

  }
  else {
    trimed_dada2 <- readRDS(trimed_dada2_fn)
  }

  names(filtFs) <- sample_names
  names(filtRs) <- sample_names

  seqtab_fn <- paste0(r_data_dir, "/", "seqtab.rds")
  denoised_fn <- paste0(r_data_dir, "/", "denoised_merged_fqs.rds")
  dada_forward_fn <- paste0(r_data_dir, "/", "data_forward.rds")
  dada_reverse_fn <- paste0(r_data_dir, "/", "data_reverse.rds")
  err_forward_fn <- paste0(r_data_dir, "/", "err_forward.rds")
  err_reverse_fn <- paste0(r_data_dir, "/", "err_reverse.rds")

  seqtab <- NULL
  dadaFs <- NULL
  dadaRs <- NULL

  if(!file.exists(denoised_fn) & paired) {

    errF <- NULL
    errR <- NULL

    if(!file.exists(err_forward_fn) & !file.exists(err_reverse_fn)) {

      errF <- learnErrors(filtFs,
                          randomize = TRUE,
                          verbose = TRUE,
                          MAX_CONSIST = 20,
                          multithread=multithread)

      saveRDS(object = errF, file = err_forward_fn)

      errR <- learnErrors(filtRs,
                          randomize = TRUE,
                          verbose = TRUE,
                          MAX_CONSIST = 20,
                          multithread=multithread)

      saveRDS(object = errR, file = err_reverse_fn)
    }
    else {

      errF <- readRDS(err_forward_fn)
      errR <- readRDS(err_reverse_fn)

    }

    if(large) {

      mergers <- vector("list", length(sample_names))
      names(mergers) <- sample_names

      dadaFs <- vector("list", length(sample_names))
      names(dadaFs) <- sample_names

      dadaRs <- vector("list", length(sample_names))
      names(dadaRs) <- sample_names

      for(s in sample_names) {

        derepFs <- derepFastq(filtFs[[s]])

        dadaF <- dada(derepFs,
                      err=errF,
                      verbose = TRUE,
                      OMEGA_C = omega_c,
                      OMEGA_A = omega_a,
                      MIN_FOLD = 1,
                      multithread = multithread)

        derepRs <- derepFastq(filtRs[[s]])

        dadaR <- dada(derepRs,
                      err=errR,
                      verbose = TRUE,
                      OMEGA_C = omega_c,
                      OMEGA_A = omega_a,
                      MIN_FOLD = 1,
                      multithread = multithread)

        merger <- mergePairs(dadaF, derepFs, dadaR, derepRs)
        mergers[[s]] <- merger

        dadaFs[[s]] <- dadaF
        dadaRs[[s]] <- dadaR

      }
      rm(derepFs)
      rm(derepRs)

      saveRDS(object = dadaFs, file = dada_forward_fn)
      saveRDS(object = dadaRs, file = dada_reverse_fn)

      seqtab <- makeSequenceTable(mergers)
      saveRDS(object = seqtab, file = seqtab_fn)
      saveRDS(object = mergers, file = denoised_fn)

    } else {

      derepFs <- derepFastq(filtFs)
      names(derepFs) <- sample_names

      derepRs <- derepFastq(filtRs)
      names(derepRs) <- sample_names

      dadaFs <- dada(derepFs,
                     err=errF,
                     pool=TRUE,
                     verbose = TRUE,
                     OMEGA_C = omega_c,
                     OMEGA_A = omega_a,
                     MIN_FOLD = 1,
                     multithread = multithread)

      dadaRs <- dada(derepRs,
                     err=errR,
                     pool=TRUE,
                     verbose = TRUE,
                     OMEGA_C = omega_c,
                     OMEGA_A = omega_a,
                     MIN_FOLD = 1,
                     multithread = multithread)

      saveRDS(object = dadaFs, file = dada_forward_fn)
      saveRDS(object = dadaRs, file = dada_reverse_fn)

      mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
      seqtab <- makeSequenceTable(mergers)

      saveRDS(object = seqtab, file = seqtab_fn)
      saveRDS(object = mergers, file = denoised_fn)
    }

  } else if(!file.exists(seqtab_fn) & !paired) {

    if(large) {
      stop("MSG: Not implemented")
      #dds <- vector("list", length(sample_names))
      #names(dds) <- sample_names
      #for(sam in sample_names) {
      #  cat("Processing:", sam, "\n")
      #  derep <- derepFastq(filts[[sam]])
      #  ddF <- dada(derepF, err=errF, multithread=TRUE)
      #  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
      #}
    }
    else {

      derepFs <- derepFastq(filtFs)
      names(derepFs) <- sample_names

      errF <- learnErrors(filtFs,
                          randomize = TRUE,
                          verbose = TRUE,
                          MAX_CONSIST = 20,
                          multithread=TRUE)

      dadaFs <- dada(derepFs,
                     err=errF,
                     pool=TRUE,
                     verbose = TRUE,
                     OMEGA_C = omega_c,
                     OMEGA_A = omega_a,
                     MIN_FOLD = 1,
                     multithread = multithread)

      saveRDS(object = dadaFs, file = dada_forward_fn)

      seqtab <- makeSequenceTable(dadaFs)
      saveRDS(object = seqtab, file = seqtab_fn)
    }

  }

  else {
    if(paired) {
      dadaFs <- readRDS(dada_forward_fn)
      dadaRs <- readRDS(dada_reverse_fn)
      mergers <- readRDS(denoised_fn)
      seqtab <- readRDS(seqtab_fn)
    }
    else {
      dadaFs <- readRDS(dada_forward_fn)
      seqtab <- readRDS(seqtab_fn)
    }
  }

  seqtab_nochim_fn <- paste0(r_data_dir, "/", "seqtab_nochim.rds")
  seqtab_nochim <- NULL

  if(!file.exists(seqtab_nochim_fn)) {
    seqtab_nochim <- removeBimeraDenovo(seqtab,
                                        multithread=multithread,
                                        #method = "consensus",
                                        method = "pooled",
                                        verbose=TRUE)

    saveRDS(object = seqtab_nochim, file = seqtab_nochim_fn)
  } else {
    seqtab_nochim <- readRDS(seqtab_nochim_fn)
  }

  # ---- sanity_check1.1 ----

  chimeras_rate <- 1 - sum(seqtab_nochim)/sum(seqtab)

  # ---- lib_size_info1 ----

  get_uniq <- function(x) {
    s <- getUniques(x) %>% sum
    return(s)
  }

  track <- NULL
  if(paired) {

    track <- cbind(trimed_dada2,
                   sapply(dadaFs, get_uniq),
                   sapply(dadaRs, get_uniq),
                   sapply(mergers, get_uniq),
                   rowSums(seqtab_nochim))

    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("lib_size",
                         "filtered",
                         "denoisedF",
                         "denoisedR",
                         "merged",
                         "nonchim")
  }
  else {

    track <- cbind(trimed_dada2,
                   sapply(dadaFs, get_uniq),
                   rowSums(seqtab),
                   rowSums(seqtab_nochim))

    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    #TODO merged column name isn't right here, but leaving here for now for simplicity
    # column name should be raw OTU's or something like that
    colnames(track) <- c("lib_size",
                         "filtered",
                         "denoisedF",
                         "merged",
                         "nonchim")
  }

  df_info <- track %>%
              as.data.frame() %>%
              rownames_to_column(var = "sample_id") %>%
              as_tibble() %>%
              arrange(lib_size)

  df_info_fn = paste0(r_data_dir, "/", "lib_size_info.tsv")

  if(!file_test("-f", df_info_fn)) {
      df_info %>% write_tsv(df_info_fn)
  }

  #df_info %>% arrange(-nonchim) %>% knitr::kable()

  return(list("info" = df_info,
              "seqtab" = seqtab,
              "seqtab_nochim" = seqtab_nochim,
              "chim_rate" = chimeras_rate))
}
