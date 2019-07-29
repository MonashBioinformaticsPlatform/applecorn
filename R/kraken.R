
mk_otus_fasta <- function(seqtab_nochim, taxtab) {

  r_data_dir <- "data"

  df_seq <- taxtab %>%
              as.data.frame() %>%
              rownames_to_column(var = "seq") %>%
              as_tibble() %>%
              select(seq) %>%
              rowwise %>%
              #mutate(md5 = substr(digest(seq, algo = 'md5', seed = 333),1,6)) %>%
              mutate(md5 = digest(seq, algo = 'md5', seed = 333)) %>%
              ungroup
              #unite(name, Phylum, Family, Genus, Species, sep = ":") %>%
              #select(name, seq)

  seqs <- DNAStringSet(df_seq$seq)
  names(seqs) <- df_seq$md5

  otus_fn <- paste0(r_data_dir, "/", "seqtab_nochim.fa")
  if(!file.exists(otus_fn)) {
    writeXStringSet(seqs, otus_fn)
  }

  return(list("df_seq" = df_seq,
              "otus_fn" = otus_fn))

}

do_taxtab_filt <- function(kraken_contams_fn, df_seq) {

  #run kraken manually to get contams.txt file
  contams <- read_tsv(kraken_contams_fn, col_names = F) %>% unlist %>% unname

  r_data_dir <- "data"
  alignment_filt_fn <- paste0(r_data_dir, "/", "alignment_filt.rds")

  if(!file.exists(alignment_filt_fn)) {
    df_seq2 <- df_seq %>% filter(!(md5 %in% contams))

    seqs <- df_seq2$seq
    names(seqs) <- df_seq2$seq

    alignment_filt <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
    saveRDS(object = alignment_filt, file = alignment_filt_fn)

  } else {
    alignment_filt <- readRDS(alignment_filt_fn)
  }

  fitGTR_fn <- paste0(r_data_dir, "/", "fitGTR.rds")

  if(!file.exists(fitGTR_fn)) {

    phang.align <- phyDat(as(alignment_filt, "matrix"), type="DNA")
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm) # Note, tip order != sequence order
    fit <- pml(treeNJ, data=phang.align)

    fitGTR <- update(fit, k=4, inv=0.2)
    fitGTR <- optim.pml(fitGTR,
                        model="GTR",
                        optInv=TRUE,
                        optGamma=TRUE,
                        rearrangement = "NNI",
                        control = pml.control(trace = 0))

    saveRDS(object = fitGTR, file = fitGTR_fn)

  } else {
    fitGTR <- readRDS(fitGTR_fn)
  }

  detach("package:phangorn", unload=TRUE)

  return(list("fitGTR" = fitGTR))

}

do_seqtab_filt <- function(kraken_contams_fn, df_seq, seqtab_nochim) {

  contams <- read_tsv(kraken_contams_fn, col_names = F) %>% unlist %>% unname

  chk <- df_seq %>%
            filter(md5 %in% contams)

  seqtab_nochim_filt_fn <- paste0(r_data_dir, "/", "seqtab_nochim_filt.rds")
  if(!file.exists(seqtab_nochim_filt_fn)) {
    seqtab_nochim_filt <- seqtab_nochim[, !(colnames(seqtab_nochim) %in% chk$seq)]
    saveRDS(object = seqtab_nochim_filt, file = seqtab_nochim_filt_fn)
  } else {
    seqtab_nochim_filt <- readRDS(seqtab_nochim_filt_fn)
  }

  chk2 <- seqtab_nochim[,chk$seq] %>%
            as.data.frame() %>%
            rownames_to_column(var = "sample") %>%
            as_tibble()

  chk3 <- chk2 %>%
            gather(seq, v, -sample) %>%
            filter(v != 0) %>%
            group_by(sample) %>%
            mutate(tot = sum(v)) %>%
            ungroup %>%
            unite(v2, v, tot, sep = ":") %>%
            group_by(seq) %>%
            mutate(name = paste0(sample, collapse = "|"),
                   v3 = paste0(v2, collapse = "|")) %>%
            ungroup %>%
            unite(name2, name, v3, sep = "=") %>%
            select(name2, seq) %>%
            distinct()

  seqs <- DNAStringSet(chk3$seq)
  names(seqs) <- chk3$name2
  contams_fn <- paste0(r_data_dir, "/", "contams.fa")

  if(!file.exists(otus_fn)) {
    writeXStringSet(seqs, contams_fn)
  }

  return(list("seqtab_nochim_filt" = seqtab_nochim_filt))
}
