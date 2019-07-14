#!/usr/bin/env Rscript

start_time <- Sys.time()
msg <- paste0("Starting applecorn run: ", start_time)
cat(msg, sep = "\n")

prog <- unlist(strsplit(commandArgs()[4], split = "--file="))[2]
origin <- dirname(prog)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  msg <- paste0("\n\n",
			   	"USAGE: ",
			   	basename(prog),
			   	" <FASTQ_DIR>",
				"\n\n")
  stop(msg)
}

raw_data <- normalizePath(args[1])
if(!file_test("-d", raw_data)) {
  msg <- paste0("\n\n",
			   	"Directory doesn't exist. Check your path",
			   	"\n--> ",
			   	raw_data,
			   	"\n\n")
  stop(msg)
}

libs <- c("R/00-libraries.R",
          "R/01-dada2.R",
          "R/02-taxa.R",
          "R/99-plan.R")

res_dir="applecorn-run"
if(!file_test("-d", res_dir)) {
   dir.create(res_dir, recursive = T)
}

setwd(res_dir)

report_fn <- "report.Rmd"
report_fn_full <- file.path(origin, report_fn)

file.copy(report_fn_full, report_fn)

libs2 <- file.path(origin, libs)

lapply(libs2, source)

# Now, your functions and workflow plan should be in your environment.
ls()

# Optionally plot the graph of your workflow.
# config <- drake_config(my_plan) # nolint
# vis_drake_graph(config)         # nolint

# Now it is time to actually run your project.
#make(applecorn,
#	 parallelism = "future",
#	 jobs = detectCores()) # Or make(my_plan, jobs = 2), etc.
make(applecorn)
## Now, if you make(whole_plan) again, no work will be done
## because the results are already up to date.
## But change the code and some targets will rebuild.
#
## Read the output report.md file
## for an overview of the project and the results.
end_time <- Sys.time()
msg <- paste0("Finished applecorn run: ", end_time)
cat(msg, sep = "\n")

wall_time <- end_time - start_time
wall_time
