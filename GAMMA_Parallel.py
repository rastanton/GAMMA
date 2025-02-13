#!/usr/bin/env python3

import sys
import glob
import subprocess
from multiprocessing import Pool

##Written by Richard Stanton (njr5@cdc.gov/stanton.rich@gmail.com)
##Runs GAMMA in parallel for all *.fasta files in a directory against a multifasta gene database in the *.fasta format
##Requires Python/3+ and blat
##Usage: $ GAMMA_Parallel.py gene_db.fasta

def GAMMA_List(input_list):
    input_fasta = input_list[0]
    input_DB = input_list[1]
    DB_Name = input_DB.split('/')[-1]
    Name = input_fasta.split('/')[-1]
    subprocess.call('GAMMA.py ' + Name + ' ' + input_DB + ' ' + Name[0:-6] + '__' + DB_Name[0:-6], shell=True)

def GAMMA_Parallel(input_DB):
    List1 = glob.glob('*.fasta')
    Run_List = []
    for files in List1:
        Run_List.append([files, input_DB])
    with Pool() as pool:
        pool.map(GAMMA_List, Run_List)

GAMMA_Parallel(sys.argv[1])
