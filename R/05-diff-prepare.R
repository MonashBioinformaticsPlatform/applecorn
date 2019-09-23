# ---- main_df_setup ----

r_data_dir <- "data"
#r_data_dir <- "data_day20"
ps_filt_fn <- paste0(r_data_dir, "/", "ps_filt.RData")
load(file = ps_filt_fn)

#ps_filt_rarefied_fn <- paste0(r_data_dir, "/", "ps_filt_rarefied.RData")
#load(file = ps_filt_rarefied_fn)

k <- 1000000

df_cnts <- ps_filt %>%
            otu_table() %>%
            as.data.frame() %>%
            rownames_to_column(var = "sample") %>%
            as_tibble() %>%
            gather(seq, counts, -sample) %>%
            group_by(sample) %>%
            mutate(tot = sum(counts)) %>%
            ungroup() %>%
            #mutate(freq = counts/tot)
            mutate(cpm = (counts/tot)*k)

#df_cnts <- ps_filt_rarefied %>% 
#              otu_table() %>% 
#              as.data.frame() %>%
#              rownames_to_column(var = "sample") %>%
#              as_tibble() %>%
#              gather(seq, counts, -sample)

df_taxa <- ps_filt %>%
#df_taxa <- ps_filt_rarefied %>%
              tax_table() %>%
              as.data.frame() %>%
              rownames_to_column(var = 'seq') %>%
              as_tibble() %>%
              rowwise %>%
              mutate(md5 = substr(digest(seq, algo = 'md5', seed = 333),1,6)) %>%
              ungroup %>%
              #mutate(sp = ifelse(Species %>% is.na(), md5, Species)) %>% 
              #unite(name, Phylum, Genus, md5, sep = "-", remove = F) %>%
              unite(name, Phylum, Class, Order, Family, Genus, md5, sep = ';') %>% 
              #mutate(phylum = substr(Phylum, 1, 3)) %>%
              #unite(name, phylum, Family, Genus, md5, sep = "-") %>%
              dplyr::select(name, seq)

meta_data <- ps_filt %>%
#meta_data <- ps_filt_rarefied %>%
                sample_data() %>%
                data.frame() %>%
                rownames_to_column(var = 'sample') %>%
                as_tibble()

#aggregating cpm normalised counts with taxa
# each row is unique taxa in a given sample
df <- full_join(df_cnts, df_taxa, by = "seq")
#here I'm adding metadata information to that data.frame
df2 <- left_join(df, meta_data, by = "sample")
# keeping df2 for legacy purposes, but this is main data.frame since
# alot of the downstream analysis will rotate around that data.frame
main_df <- df2

main_df_fn <- paste0(r_data_dir, "/", "main_df.RData")

if(!file.exists(main_df_fn)) {
  save(main_df, file = main_df_fn)
} else {
  load(main_df_fn)
}

# ---- corr1 ----

n_ranked <- 10

top_ranked <- df2 %>%
                select(name, freq) %>%
                group_by(name) %>%
                summarise(freq = sum(freq)) %>%
                arrange(-freq) %>%
                head(n_ranked) %>%
                select(name) %>%
                unlist %>%
                unname

p_taxa <- df2 %>%
            #select(sample, freq, name, contains("BP"), -contains("test")) %>%
            mutate(s4fd_chk = ifelse(S4FD > 0.57, "upper", "lower")) %>%
            filter(!is.na(S4FD)) %>%
            #mutate(s4fd_chk = ifelse(NightHR > 69.79, "upper", "lower")) %>%
            #filter(!is.na(NightHR)) %>%
            mutate(most_abundant = ifelse(name %in% top_ranked, name, "other")) %>%
            filter(most_abundant != "other") %>%
            ggplot(aes(sample, freq, fill = most_abundant)) +
              geom_bar(stat = "identity") +
              scale_fill_brewer(type = "div") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              facet_wrap(~s4fd_chk, scales = "free_x")
              #facet_wrap(~hypertension, scales = "free_x")
p_taxa

#NOTE
# df2 table contains a few (166 to be exact) taxa that cotain zero counts.
# The reason that there is a taxa with zero counts is because those are rare
# types and needed to occur in at least one sample. I tracked down one specific
# sequence and it was associate with sample 70 and it had 3 counts. That sample
# was drop due to the fact that I didn't have either metadata or FFQ.
# Therefore those taxa that have zero counts in df2 are safe to drop for further analysyis

df_filt <- df2 %>%
            group_by(name) %>%
            mutate(cnts_chk = sum(counts)) %>%
            ungroup %>%
            filter(cnts_chk != 0)

cor_dat <- df_filt %>%
              select(sample, freq, name, contains("BP"), -contains("test")) # %>%
              spread(name, freq) %>%
              column_to_rownames(var = "sample")

#cor_mat <- cor(cor_dat)
cor_mat <- cor(cor_dat[,1:8], cor_dat[,9:ncol(cor_dat)], method = "pearson")

cor_df <- cor_mat %>%
            as.data.frame() %>%
            rownames_to_column(var = "measures") %>%
            as_tibble %>%
            gather(name, r, -measures) %>%
            mutate(bp_type2 = ifelse(grepl("SBP", x = measures), "systolic", "dystolic")) %>%
            mutate(bp_type = ifelse(grepl("Day", x = measures), "day",
                                     ifelse(grepl("Night", x = measures), "night",
                                            ifelse(grepl("Office", x = measures), "office", "overall"))))

p1_cor <- cor_df %>%
            #filter(r > 0.35, grepl("Night", x = measures)) %>%
            #filter(r > 0.25) %>%
            ggplot(aes(name, r, color = bp_type)) +
              geom_point() +
              #scale_color_brewer(type = "div") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              facet_wrap(~bp_type2, ncol = 1)
p1_cor

approx_chk <- function(r, chk) {

  r <- round(r, digits = 3)
  chk <- round(chk, digits = 3)

  k <- 0.02

  #print(r-k)
  #print(r+k)

  if(r == chk) {
    return(r)
  }
  else if(chk > r-k & chk < r+k) {
    return(chk)
  }
  else {
    return(0)
  }
}

cor_df2 <- cor_df %>%
            #filter(r > 0) %>%
            filter(r > 0.25) %>%
            #filter(r < -0.25) %>%
            #filter(r < -0.25 | r > 0.25) %>%
            #filter(measures == "NightDBP") %>%
            group_by(measures) %>%
            mutate(chk = sort(table(r), decreasing = T)[1] %>% names %>% as.numeric()) %>%
            #NOTE this nesting and unnesting has problems of duplicaing rows
            # essentially rows will be duplicates by the length of the list
            # going to drop this for now, and gonna stick with the top one for now
            # but this needs to be mentioned upfront
            #mutate(chk = list(sort(table(r), decreasing = T)[1:5] %>% names %>% as.numeric())) %>%
            #unnest() %>%
            ungroup %>%
            rowwise %>%
            mutate(r_cln = approx_chk(r, chk)) %>%
            ungroup

p2_cor <- cor_df2 %>%
            filter(r_cln != 0) %>%
            #filter(bp_type2 == "night") %>%
            ggplot(aes(name, r_cln, color = bp_type)) +
              geom_point() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              facet_wrap(~bp_type2, ncol = 1)

p2_cor

dyst_taxa <- cor_df2 %>%
                filter(bp_type2 == "dystolic",
                       bp_type == "night") %>%
                select(name) %>%
                unlist %>%
                unname

syst_taxa <- cor_df2 %>%
                filter(bp_type2 == "systolic",
                      bp_type == "night") %>%
                select(name) %>%
                unlist %>%
                unname

night_targets <- intersect(dyst_taxa, syst_taxa)

p_taxa2 <- df2 %>%
            select(sample, freq, name, contains("BP"), -contains("test")) %>%
            mutate(most_abundant = ifelse(name %in% night_targets, name, "other")) %>%
            filter(most_abundant != "other") %>%
            ggplot(aes(sample, freq, fill = most_abundant)) +
              geom_bar(stat = "identity", show.legend = FALSE) +
              #scale_fill_brewer(type = "div") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_taxa2

