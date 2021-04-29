#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Both MRI and demographics files are required to render the report..n", call.=FALSE)
}

Sys.setenv(RSTUDIO_PANDOC="/mnt/tigrlab/quarantine/rstudio/1.1.447/build/rstudio-1.1.447/bin/pandoc")
Sys.getenv("RSTUDIO_PANDOC")

library(renv)
renv::restore()

rmarkdown::render('OPTIMUM-Neuro_tracking.Rmd', 'html_document', 
			params=list(mri_file=args[1], demo_file=args[2]),
			output_file='index.html', output_dir='..')
