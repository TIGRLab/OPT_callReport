#!/bin/bash -l 

module load lab-code R
cwd=$(pwd)
base_dir=/archive/data/OPT/OPT_callReport

cd $base_dir/bin

python generate_report.py

Rscript --vanilla rmd2html.R OPT_report_mri.csv OPT_report_demo.csv

cd ..

# Start ssh agent
eval `ssh-agent`
ssh-add /home/clevis/.ssh/id_ed25519

git add index.html
git commit -m "Auto-updating recruitment doc"
git push origin master

cd $cwd

