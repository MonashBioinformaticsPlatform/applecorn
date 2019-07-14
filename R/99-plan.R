
applecorn <- drake_plan(
  dada2_res = run_dada2(raw_data),
  taxa_res = run_taxa(),
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE)
)
