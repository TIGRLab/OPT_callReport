#!/bin/bash -l 

module load lab-code R

python generate_report.py

Rscript --vanilla rmd2html.R OPT_report_mri.csv OPT_report_demo.csv

cd ..

git add index.html
git commit -m "Auto-updating recruitment doc"
git push origin master


