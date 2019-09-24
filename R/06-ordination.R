# ---- ordination_init ----

mk_ordination <- function(ps_filt, dist = "wunifrac") {

  #options :
  # - wunifrac
  # - unifrac
  # - jacard
  # - bray
  pslog_filt <- transform_sample_counts(ps_filt, function(x) log(1 + x))

  # ---- do_filt ----

  #Note need to think here about the data.
  # the transform_sample_counts() function under the hood simply
  # runs apply with a given function on otu_table. It feels a bit simplified since we are
  # not normalising for the library sizes here. I think distance plots PCA/PCoA
  # would be more informative when libs had been normalised for depth..?

  # this can transform otu_table into workable form, depending on what you want todo
  # you might not even need to spread there.

  ord_log <- ordinate(pslog_filt, method = "PCoA", distance = dist)
  plt <- plot_ordination(pslog_filt,
                         ord_log,
                         type = "samples",
                         color = "treat",
                         #shape = "cage",
                         axes=1:2) +
                          geom_point(size = 4) +
                          #geom_text(aes(label = sample_id),
                          #          hjust = 1.3,
                          #          angle = -15,
                          #          size = 3) +
                          ggtitle(dist)
                          #scale_color_manual(values = c("#94641F", "#C5944E", "#008000", "#7CFC00"))# +

  return(list("ord_plt" = plt))

}
