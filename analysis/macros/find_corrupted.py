#!/usr/bin/env python
import os
import sys

list_corrupted = []
for output in os.popen('grep \'rror\' logs/condor/run/err/*').read().split("\n"):
    err_file = output.split(":")[0]
    if '.stderr' not in err_file: continue
    for line in open(err_file).readlines():
        if '.root' not in line: continue
        rootfile = line.strip().split('root://')[1]
        if 'root://'+rootfile in list_corrupted: continue
        list_corrupted.append('root://'+rootfile)

if len(list_corrupted)<1:
    sys.exit('No corrupted files')

list = "data/corrupted_files.txt"
with open(list, "w") as fout:
    fout.writelines(list_corrupted)
