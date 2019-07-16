# ---- get_seqtab ----

#fqs <- "../raw-data"

run_dada2 <- function(fqs) {

  miseq_path <- fqs
  fns <- sort(list.files(miseq_path, full.names = TRUE))
  fnFs <- fns[grepl("_R1_trimed", fns)]
  fnRs <- fns[grepl("_R2_trimed", fns)]

  r_data_dir <- "data"
  if(!file_test("-d", r_data_dir)) {
    dir.create(r_data_dir, recursive = T)
  }

  filt_path <- "filtered"
  if(!file_test("-d", filt_path)) {
    dir.create(filt_path, recursive = T)
  }

  # check the quality
  #ii <- sample(length(fnFs), 3)
  #for(i in ii) {
  #  print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd"))
  #}
  #for(i in ii) {
  #  print(plotQualityProfile(fnRs[i]) + ggtitle("Rev"))
  #}

  filtFs <- file.path(filt_path, basename(fnFs))
  filtRs <- file.path(filt_path, basename(fnRs))
  sample_names <- sapply(strsplit(basename(filtFs), "_R"), `[`, 1)

  trimed_dada2_fn = paste0(r_data_dir, "/", "trimed_dada2.rds")

  trimed_dada2 <- NULL

  # Trim
  if(!file.exists(trimed_dada2_fn)) {
    trimed_dada2 <- filterAndTrim(fnFs,
                                  filtFs,
                                  fnRs,
                                  filtRs,
                                  trimLeft=10,
                                  truncLen=c(245, 200),
                                  maxN=0,
                                  maxEE=c(2, 2),
                                  truncQ=2,
                                  rm.phix=TRUE,
                                  compress=TRUE,
                                  multithread=F)

    saveRDS(object = trimed_dada2, file = trimed_dada2_fn)
  } else {
    trimed_dada2 <- readRDS(trimed_dada2_fn)
  }

  denoised_fn = paste0(r_data_dir, "/", "denoised_merged_fqs.rds")
  dada_forward_fn = paste0(r_data_dir, "/", "data_forward.rds")
  dada_reverse_fn = paste0(r_data_dir, "/", "data_reverse.rds")

  if(!file.exists(denoised_fn)) {

    # Dereplicate
    derepFs <- derepFastq(filtFs)
    derepRs <- derepFastq(filtRs)

    names(derepFs) <- sample_names
    names(derepRs) <- sample_names

    # Error estimation
    n_samples <- fnFs %>% length
    #ddF <- dada(derepFs[1:n_samples], err=NULL, selfConsist=TRUE, multithread = TRUE)
    #ddR <- dada(derepRs[1:n_samples], err=NULL, selfConsist=TRUE, multithread = TRUE)

    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtFs, multithread=TRUE)

    #plotErrors(ddF)
    #plotErrors(ddR)

    #NOTE: this step is slow minutes to hours, expect at least 30 minutes
    dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread = TRUE)
    dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread = TRUE)

    saveRDS(object = dadaFs, file = dada_forward_fn)
    saveRDS(object = dadaRs, file = dada_reverse_fn)

    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
    saveRDS(object = mergers, file = denoised_fn)

  } else {
    # mergers is an R class list that hold denoised sequences
    # from all samples. The lenght (the number of items in the list)
    # equals to the number of samples i.e merged fastq files
    # in this case we had 69 sample, there fore mergers %>% length == 69
    mergers <- readRDS(denoised_fn)
    # each item in the list is a data.frame that has set number of columns - nine
    # and variable number of rows, each row is a "unique" feature found in that sample
    # sequence:  char
    # abundance: int
    # forward:   int
    # reverse:   int
    # nmatch:    int
    # nmismatch: int
    # nindel:    int
    # prefer:    num
    # accept:    logi
    dadaFs <- readRDS(dada_forward_fn)
    dadaRs <- readRDS(dada_reverse_fn)
  }

  # make table of counts
  seqtab <- makeSequenceTable(mergers)
  # remove chimeras
  #NOTE: this was rather confusing matrix object
  # basically rownames are sample names i.e each row is a single sample
  # whereas colnames are individual sequences found in that sample
  # and the cell is the actual raw counts for that sequence in that sample
  # I guess the confusion was due to the fact that I've used my matricies in the
  # other form where columns are the sample names and rows are the individual features
  # seqtab %>% colnames
  # seqtab %>% rownames
  #WARNING: slow step, 20-30 minutes

  seqtab_fn <- paste0(r_data_dir, "/", "seqtab.rds")
  seqtab_nochim_fn <- paste0(r_data_dir, "/", "seqtab_nochim.rds")

  if(!file.exists(seqtab_nochim_fn)) {
    seqtab_nochim <- removeBimeraDenovo(seqtab)
    saveRDS(object = seqtab, file = seqtab_fn)
    saveRDS(object = seqtab_nochim, file = seqtab_nochim_fn)
  } else {
    trimed_dada2 <- readRDS(trimed_dada2_fn)
    seqtab <- readRDS(seqtab_fn)
    seqtab_nochim <- readRDS(seqtab_nochim_fn)
  }

  # ---- sanity_check1.1 ----

  chimeras_rate <- 1 - sum(seqtab_nochim)/sum(seqtab)

  # ---- lib_size_info1 ----

  get_uniq <- function(x) {
    s <- getUniques(x) %>% sum
    return(s)
  }

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
  rownames(track) <- sample_names


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
	      "seqtab_nochim" = seqtab_nochim,
	      "chim_rate" = chimeras_rate))
}
