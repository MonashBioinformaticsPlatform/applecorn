
# ---- taxonomy_assignment ----

do_taxa_ann <- function(seqtab_nochim, ref_fasta) {

  # total number of OTUs
  seqtab_nochim %>% dim
  tot_seqs <- seqtab_nochim %>% ncol

  seqs_len <- seqtab_nochim %>%
                colnames() %>%
                lapply(nchar) %>%
                unlist

  m <- seqs_len %>% mean %>% round(2)
  s <- seqs_len %>% sd %>% round(2)

  msg <- paste0("MSG: Total number of unique OTUs identified ",tot_seqs,
                "\n",
                "MSG: Mean length ", m,
                " and standard deviation ", s)
  cat(msg, sep = '\n')

  r_data_dir <- "data"
  taxtab_fn <- paste0(r_data_dir, "/", "taxtab.rds")

  if(!file.exists(taxtab_fn)) {

    msg <- paste0("MSG: Starting taxonomic annotation using, ", ref_fasta %>% basename, " database")
    cat(msg, sep = '\n')

    #WARNING: slow step 5-10 minutes
    taxtab <- assignTaxonomy(seqtab_nochim,
                             refFasta = ref_fasta,
                             tryRC = TRUE,
                             multithread = TRUE)

    #colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    saveRDS(object = taxtab, file = taxtab_fn)

  } else {
    taxtab <- readRDS(taxtab_fn)
  }

  return(taxtab)
}

mk_tree <- function(seqtab_nochim, taxatab) {

  r_data_dir <- "data"
  alignment_fn <- paste0(r_data_dir, "/", "alignment.rds")

  if(!file.exists(alignment_fn)) {

    seqs <- getSequences(seqtab_nochim)
    names(seqs) <- seqs # This propagates to the tip labels of the tree
    alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
    saveRDS(object = alignment, file = alignment_fn)

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

    msg <- paste0("MSG: Starting pml fit")
    cat(msg, sep = '\n')

    fit <- pml(treeNJ, data = phang.align)

    msg <- paste0("MSG: Starting pml optimisation ",
                   "This is typically very slow step ",
                   "Can be hours", sep = "\n")

    cat(msg, sep = '\n')

    fitGTR <- update(fit, k=4, inv=0.2)
    #NOTE: this step is by far the longest step over an hour
    #NOTE Okay I can massively speed things up here by using NNI rearrangement rather then stochastic
    # Qiime2 uses FastTree and my understanding is that FastTree does the same i.e NNI
    # This is from phangorn docs https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.pdf
    # For larger trees the NNI rearrangements often get stuck in local maxima.
    # We added two stochatic algorithm to improve topology search.
    # So I think for my purpose it is fine to use NNI, I don't particular care about optimal topology
    # although note sure if topology has an effect on UniFrac i.e how taxa are related..
    #TODO? The only other thing that I should investigate is to try to filter branches at treeNJ stage before optimasing
    # this in theory should help ..?
    fitGTR <- optim.pml(fitGTR,
                        model="GTR",
                        optInv=TRUE,
                        optGamma=TRUE,
                        #rearrangement = "stochastic",
                        rearrangement = "NNI",
                        control = pml.control(trace = 0))

    saveRDS(object = fitGTR, file = fitGTR_fn)

  } else {
    fitGTR <- readRDS(fitGTR_fn)
  }

  detach("package:phangorn", unload=TRUE)

  return(fitGTR)
}
