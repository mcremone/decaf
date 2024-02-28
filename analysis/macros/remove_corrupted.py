#!/usr/bin/env python
import json
from optparse import OptionParser
import os
import sys

parser = OptionParser()
parser.add_option('-m', '--metadata', help='metadata', dest='metadata', default='2018')
parser.add_option('-s', '--save', help='save', dest='save', default=None)
parser.add_option('-l', '--lists', help='lists', dest='lists', default=None)
(options, args) = parser.parse_args()

list_corrupted = []
for output in os.popen('grep \'rror\' logs/condor/run/err/*').read().split("\n"):
    err_file = output.split(":")[0]
    if '.stderr' not in err_file: continue
    for line in open(err_file).readlines():
        if '.root' not in line: continue
        rootfile = line.strip().split('root://')[1]
        if 'root://'+rootfile in list_corrupted: continue
        list_corrupted.append('root://'+rootfile)
if options.save:
    with open('/data/'+options.save+'_list_corrupted.pkl', 'wb') as fp:
        pickle.dump(list_corrupted, fp)    
if options.lists:
    with open ('/data/'+options.lists+'_list_corrupted.pkl', 'rb') as fp:
        list_corrupted = pickle.load(fp)
if len(list_corrupted)<1:
    sys.exit('No corrupted files')

dictionary={}
with open(options.metadata) as fin:
    dictionary.update(json.load(fin))

list_datasets = []
for key in dictionary:
    for rootfile in list_corrupted:
        if rootfile in dictionary[key]["files"]:
            if key not in list_datasets: list_datasets.append(key)
            dictionary[key]["files"].remove(rootfile)
print("Found corrupted file in", ','.join(list_datasets))


with open(options.metadata, "w") as fout:
    json.dump(dictionary, fout, indent=4)
