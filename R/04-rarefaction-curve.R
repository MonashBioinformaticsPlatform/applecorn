# ---- rare_curve1 ----

# rarefaction code from here https://github.com/joey711/phyloseq/issues/143

#NOTE in future ideally want to drop plyr and reshap2 packages in favour tidyverse
# although not sure if ldply has change and/or part of tidyverse
# this is mainly to keep the code

est_rare_richness <- function(ps_dat, measures, depth) {

  ps_dat_filt <- prune_samples(sample_sums(ps_dat) >= depth, ps_dat)

  rarified_ps_dat <- rarefy_even_depth(physeq = ps_dat_filt,
                                       sample.size = depth,
                                       replace = FALSE,
                                       verbose = FALSE)

  alpha_diversity <- estimate_richness(rarified_ps_dat, measures = measures)

  res <- alpha_diversity %>%
            rownames_to_column(var = "sample") %>%
            as_tibble() %>%
            gather(measure, "Alpha_diversity", -sample)

  return(res)
}

mk_rare_curves <- function(ps_dat, measures, depths, parallel = TRUE) {

  # set parallel options if libraryd
  if(parallel) {
    paropts  <- list(.packages=c("phyloseq", "reshape2"))
  } else {
    paropts  <- NULL
  }

  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply

  rare_curve <- purrr::map_dfr(depths,
                               est_rare_richness,
                               ps_dat = ps_dat,
                               measures = measures,
                               .id = "depth") %>%
                       mutate(depth = depth %>% as.numeric)

  return(rare_curve)
}

# NOTE think about smarter way to set the sampling sizes
# sampling size should be dynamically identified from the data itself
# since the range will depend on the maxmimum reads depth from all the samples
# in this particular study
# > df_info %>% select(nonchim) %>% max
# [1] 30369
# therefore sampling_sizes should go upto 35000 for example
#df_info data frame should be availibale here, perhaps use that
#
# Just want to add to the comments above. I tested it our with 35K
# and it doesn't appear to work, which makes sense since we can't draw 35k sequences
# from a pool of 31k

mk_range <- function(d_min, d_max) {

  if(d_min > 1000) {
    d_min <- 1000
  }

  #TODO assert(d_min >= 0)

  #adjust binary search to d_min
  d <- d_max-d_min
  # binary search, worth case senario. I simply want to get even number of steps
  n_steps <- floor(log2(d))
  steps <- floor(d/n_steps)

  range <- seq(d_min, d_max, steps)

  return(range)
}

do_rare_curve <- function(ps_filt, n_reps = 10) {

  set.seed(42)

  d_max <- ps_filt %>% sample_sums %>% max
  d_min <- ps_filt %>% sample_sums %>% min

  #depths <- rep(sampling_sizes, each = 10)

  sizes <- mk_range(d_min, d_max)
  depths <- rep(sizes, each = n_reps)

  rare_curve <- mk_rare_curves(ps_filt,
                               c('Observed', 'Shannon'),
                               depths = depths)

  # mean in this case referes to the average number of OTUs detected from
  # random sampling. Remember that we do 10 random trials, at each depth.
  # if we were to look at the distribution of all 10 trials, it should be pretty narrow
  # I guess we have sd there as well.
  # and so doing this
  # > rare_curve_sum %>% as_tibble() %>% select(Depth, Alpha_diversity_sd) %>% plot
  # tells us (as expected) that with bigger sample sd decreases

  rare_curve_sum <- rare_curve %>%
                        group_by(sample, depth, measure) %>%
                        summarise(alpha_mean = mean(Alpha_diversity),
                                  alpha_sd = sd(Alpha_diversity))

  # ---- rare_curve2 ----

  #NOTE this is more general info about rarefaction curves.
  # The purpose of the such curve is to see the discoveribility of
  # taxa (OTUs/ASVs) at different sequencing depth. We simulate those
  # sequencing depths by simply taking a random sample from our total pool
  # and counting how many OTUs we've discovered at that depth.
  # However typically our data is grouped in several different ways
  # by treatment types, or by gender or someother way. Rarefaction curvers
  # are build in the context of a particular group.
  # Once again for the sake of the pipeline one can either grab some default
  # or provide a shiny app.

  p_rare <- rare_curve_sum %>%
              mutate(depth_fct = depth %>% as.factor) %>%
              ggplot(aes(depth,
                         alpha_mean,
                         ymin = alpha_mean - alpha_sd,
                         ymax = alpha_mean + alpha_sd,
                         colour = sample,
                         group = sample)) +
                  geom_line(show.legend = F) +
                  geom_pointrange(show.legend = F) +
                  scale_x_continuous(breaks = sizes,
                                     limits = c(min(sizes), max(sizes))) +
                  facet_wrap(~measure, scales = 'free_y') +
                  ylab(paste0("Alpha diversity, mean across ", n_reps, " reps")) +
                  ggtitle("Rarefaction curves") +
                  theme(axis.text.x = element_text(angle = 90, hjust = 1))

  return(p_rare)
}

# Perhaps rarefaction curve most of (by default) should always
# be per sample, i.e which samples appear to "under sample" OTUs and
# which samples seem to "loose (gain?) power" at a greater seq depth

# jsut FYI https://github.com/gauravsk/ranacapa
# probably should use that in the future? need to check the code first
# also note vegan:rarecurve() function
