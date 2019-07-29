#!/usr/bin/env Rscript

start_time <- Sys.time()
msg <- paste0("Starting applecorn run: ", start_time)
cat(msg, sep = "\n")

prog <- unlist(strsplit(commandArgs()[4], split = "--file="))[2]
base_prog <- basename(prog)
origin <- dirname(prog)

config_fn <- "config.yml"
if(!file.exists(config_fn)) {

  msg <- paste0("Initiating config file. Please review ", config_fn, " file and restart ", base_prog)
  cat(msg, sep = "\n")

  config_full_fn <- file.path(origin, config_fn)
  file.copy(config_full_fn, config_fn)
}

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  msg <- paste0("\n\n",
                "USAGE: ",
                base_prog,
                " <CONFIG>",
                "\n\n")
  stop(msg)
}

config_fn <- normalizePath(args[1])
config <- yaml::read_yaml(config_fn)

debug <- as.logical(config$debug)

if(debug) {
  print(origin)
  print(config)
}

raw_data <- normalizePath(config$raw_data)
metadata <- normalizePath(config$metadata)
classifier <- normalizePath(config$classifier)
res_dir <- config$output_dir

if(!file_test("-d", res_dir)) {
   dir.create(res_dir, recursive = T)
}

msg <- paste0("MSG: Setting working to ", res_dir)
cat(msg, sep = "\n")

#TODO: there is a smarter way of doing this, probably just use config list?
f <- c(raw_data, metadata, classifier)

file_chk <- function(f) {
  if(!file.exists(f)) {
    msg <- paste0("\n\n",
                  "File doesn't exist. Check your path",
                  "\n--> ",
                  f,
                  "\n\n")
    stop(msg)
  }
}

quite <- lapply(f, file_chk)

setwd(res_dir)
if(debug) {
  getwd()
}

libs <- list.files(file.path(origin, "R"), full.names = TRUE)
if(debug) {
  print(libs)
}
lapply(libs, source)

report_fn <- "report.Rmd"
report_fn_full <- file.path(origin, report_fn)

report_dat <- readLines(report_fn_full)
report_dat <- gsub("DATABASE_PLACE_HOLDER", basename(classifier), report_dat)
writeLines(report_dat, report_fn)

# Now, your functions and workflow plan should be in your environment.
ls()

plan <- file.path(origin, "plan.R")
source(plan)

if(debug) {
  print(applecorn)
}

make(applecorn)

end_time <- Sys.time()
msg <- paste0("Finished applecorn run: ", end_time)
cat(msg, sep = "\n")

wall_time <- end_time - start_time
wall_time
