
applecorn <- drake_plan(
  report = knit(knitr_in("report.Rmd"), file_out("report.md"), quite = TRUE),
  dada2 = dada2("raw-data")
)
