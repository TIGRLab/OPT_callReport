#!/usr/bin/env python3
""" Generate the recruitment tracking html doc used for the OPT PI call. """

import csv
import os
import logging
import requests
from subprocess import Popen

import datman.config
from redcap import Project, RedcapError

logger = logging.getLogger(os.path.basename(__file__))


def read_token(token_file):
    if not os.path.isfile(token_file):
        logger.error('REDCap token file: {} not found'.format(token_file))
        raise IOError

    with open(token_file, 'r') as token_file:
        token = token_file.readline().strip()
        return token

def main():
    cfg = datman.config.config(study='OPT')

    api_url = cfg.get_key('RedcapUrl')
    redcap_url = api_url.replace('/api/', '/')

    dir_meta = cfg.get_path('meta')
    token_path = os.path.join(dir_meta, cfg.get_key('RedcapToken'))
    token = read_token(token_path)

    project = Project(api_url, token)
    
    print('Exporting reports from REDCap...')
    demo_data = project.export_reports(format='csv', report_id='45026')
    mri_data  = project.export_reports(format='csv', report_id='45089')
    
    mri_filepath = 'OPT_report_mri.csv'
    demo_filepath = 'OPT_report_demo.csv'
    
    print(f'Writing MRI data to {mri_filepath}...')
    with open(mri_filepath, 'w', newline='') as csvfile:
        csvfile.write(mri_data)

    print(f'Writing demographics data to {demo_filepath}...')
    with open(demo_filepath, 'w', newline='') as csvfile:
        csvfile.write(demo_data)
    
    print('Calling Rscript to knit html page...')
    Popen(['Rscript', '--vanilla', 'rmd2html.R',
        mri_filepath, demo_filepath], shell=False)
    
    print('Done!')

if __name__=='__main__':
    main()
