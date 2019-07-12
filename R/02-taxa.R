
# ---- taxonomy_assignment ----

ptm <- proc.time()

r_data_dir <- "../data"

ref_fasta <- "~/taxonomic_classifiers/rdp_train_set_16.fa.gz"
#ref_fasta <- "~/taxonomic_classifiers/gg_13_8_train_set_97.fa.gz"

seqtab_nochim_fn <- paste0(r_data_dir, "/", "seqtab_nochim.rds")

if(!file.exists(seqtab_nochim_fn)) {
  stop("ERROR: This can't happen")
}

seqtab_nochim <- readRDS(seqtab_nochim_fn)

# total number of OTUs
seqtab_nochim %>% dim
tot_seqs <- seqtab_nochim %>% ncol

seqs_len <- seqtab_nochim %>%
              colnames() %>%
              lapply(nchar) %>%
              unlist

m <- seqs_len %>% mean
s <- seqs_len %>% sd

msg <- paste0("MSG: Total number of unique OTUs identified, ",tot_seqs,
              " mean ", m,
              " and standard deviation ", s,
              " for the length of the sequence")
cat(msg, sep = '\n')

taxtab_fn <- paste0(r_data_dir, "/", "taxtab.rds")

if(!file.exists(taxtab_fn)) {
  
  msg <- paste0("MSG: Starting taxonomic annotation using, ", ref_fasta %>% basename, " database")
  cat(msg, sep = '\n')
  
  #WARNING: slow step 5-10 minutes
  taxtab <- assignTaxonomy(seqtab_nochim,
                           refFasta = ref_fasta,
                           tryRC = TRUE,
                           multithread = TRUE)

  colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  saveRDS(object = taxtab, file = taxtab_fn)
  
  ptm <- proc.time() - ptm
  
  msg <- paste0("MSG: Taxonomic annotation took, ", round(ptm[['elapsed']]/60, 3), " minutes")
  cat(msg, sep = '\n')

} else {
  taxtab <- readRDS(taxtab_fn)
}

alignment_fn <- paste0(r_data_dir, "/", "alignment.rds")

if(!file.exists(alignment_fn)) {
  
  # Multiple sequence alignment and phylotree
  seqs <- getSequences(seqtab_nochim)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  
  msg <- paste0("MSG: Starting multiple sequence alignment with ", tot_seqs, " sequnces")
  cat(msg, sep = '\n')
  
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
  saveRDS(object = alignment, file = alignment_fn)

  ptm <- proc.time() - ptm
  msg <- paste0("MSG: Multiple sequence alignment took ", round(ptm[['elapsed']]/60, 3), " minutes")
  cat(msg, sep = '\n')
  
} else {
  alignment <- readRDS(alignment_fn)
}
  
fitGTR_fn <- paste0(r_data_dir, "/", "fitGTR.rds")

if(!file.exists(fitGTR_fn)) {
  
  msg <- paste0("MSG: Building phylogenetic tree using NJ method")
  cat(msg, sep = '\n')
  
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  
  ptm <- proc.time() - ptm
  msg <- paste0("MSG: NJ clustering took ", round(ptm[['elapsed']]/60, 3), " minutes")
  cat(msg, sep = '\n')
  
  msg <- paste0("MSG: Starting pml fit")
  cat(msg, sep = '\n')
  
  fit = pml(treeNJ, data=phang.align)

  ptm <- proc.time() - ptm
  msg <- paste0("MSG: Fitting pml took ", round(ptm[['elapsed']]/60, 3), " minutes")
  cat(msg, sep = '\n')
  
  ## negative edges length changed to 0!

  msg <- paste0("MSG: Starting pml optimisation ",
                 "This is typically very slow step ",
                 "Can be hours", sep = "\n")
  
  cat(msg, sep = '\n')
  
  fitGTR <- update(fit, k=4, inv=0.2)
  #NOTE: this step is by far the longest step over an hour
  fitGTR <- optim.pml(fitGTR,
                      model="GTR",
                      optInv=TRUE,
                      optGamma=TRUE,
                      rearrangement = "stochastic",
                      control = pml.control(trace = 0))
  
  #h <- 5 
  #lim <- h * 3600
  #
  #withTimeout(fitGTR <- optim.pml(fitGTR,
  #                    model="GTR",
  #                    optInv=TRUE,
  #                    optGamma=TRUE,
  #                    rearrangement = "stochastic",
  #                    control = pml.control(trace = 0)), 
  #            timeout = lim,
  #            onTimeout = "error")

  ptm <- proc.time() - ptm
  msg <- paste0("MSG: Tree optimisation took ", round(ptm[['elapsed']]/60, 3), " minutes")
  cat(msg, sep = '\n')
  
  saveRDS(object = fitGTR, file = fitGTR_fn)
  
} else {
  fitGTR <- readRDS(fitGTR_fn)
}

detach("package:phangorn", unload=TRUE)