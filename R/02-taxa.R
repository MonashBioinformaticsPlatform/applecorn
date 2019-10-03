
# ---- taxonomy_assignment ----

do_taxa_ann <- function(seqtab_nochim,
                        ref_fasta,
                        r_data_dir = "data") {

  taxtab_fn <- paste0(r_data_dir, "/", "taxtab.rds")

  if(!file.exists(taxtab_fn)) {

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

fit_tree <- function(seqtab_nochim,
                     taxatab,
                     r_data_dir = "data") {

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

    phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
    dm <- phangorn::dist.ml(phang.align)
    treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order

    fit <- phangorn::pml(treeNJ, data = phang.align)

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
    fitGTR <- phangorn::optim.pml(fitGTR,
                                  model="GTR",
                                  optInv=TRUE,
                                  optGamma=TRUE,
                                  #rearrangement = "stochastic",
                                  rearrangement = "NNI",
                                  control = phangorn::pml.control(trace = 0))

    saveRDS(object = fitGTR, file = fitGTR_fn)

  } else {
    fitGTR <- readRDS(fitGTR_fn)
  }

  #detach("package:phangorn", unload=TRUE)

  return(fitGTR)
}
