#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Both MRI and demographics files are required to render the report..n", call.=FALSE)
}

Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
Sys.getenv("RSTUDIO_PANDOC")
rmarkdown::render('OPTIMUM-Neuro_tracking.Rmd', 'html_document', params=list(mri_file=args[1], demo_file=args[2]))