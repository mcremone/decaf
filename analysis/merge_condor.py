#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
import numexpr
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', '--folder', help='folder', dest='folder')
parser.add_option('-v', '--variable', help='variable', dest='variable', default='')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../cmssw_11_3_4.tgz '
              '--exclude=\'src/decaf/analysis/logs\' '
              '--exclude=\'src/decaf/analysis/plots\' '
              '--exclude=\'src/decaf/analysis/datacards\' '
              '--exclude=\'src/decaf/analysis/results\' '
              '--exclude=\'src/decaf/analysis/hists/*/*.futures\' '
              '--exclude=\'src/decaf/analysis/hists/*/*.merged\' '
              '../../../../CMSSW_11_3_4')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../pylocal_3_8.tgz -C ~/.local/lib/python3.8/ site-packages')

if options.cluster == 'kisti':
    if options.copy:
        os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/cmssw_11_3_4.tgz')
        print('cmssw removed')
        os.system('xrdcp -f ../../../../cmssw_11_3_4.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/cmssw_11_3_4.tgz')
        os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/pylocal_3_8.tgz') 
        print('pylocal removed')
        os.system('xrdcp -f ../../../../pylocal_3_8.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/pylocal_3_8.tgz')
    jdl = """universe = vanilla
Executable = merge.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = merge.sh, /tmp/x509up_u556950957
Output = logs/condor/merge/out/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/merge/err/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/merge/log/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE).merged=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE).merged"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(CLUSTER) $ENV(USER)
JobBatchName = $ENV(VARIABLE)
accounting_group=group_cms
request_cpus = 16
Queue 1"""

if options.cluster == 'lpc':
    if options.copy:
        os.system('xrdcp -f ../../../../cmssw_11_3_4.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/cmssw_11_3_4.tgz')
        os.system('xrdcp -f ../../../../pylocal_3_8.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal_3_8.tgz')
    jdl = """universe = vanilla
Executable = merge.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = merge.sh
Output = logs/condor/merge/out/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/merge/err/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/merge/log/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE).merged=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE).merged"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(CLUSTER) $ENV(USER)
request_cpus = 16
Queue 1"""

jdl_file = open("merge.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

tag=options.folder.split('/')[-1]

variables = []
for filename in os.listdir(options.folder):
          if '.reduced' not in filename: continue
          if filename.split('--')[0] not in variables: variables.append(filename.split('--')[0])

for variable in variables:
    if options.variable and options.variable not in variable: continue
    if options.variable:
        if not any(_variable==variable for _variable in options.variable.split(',')): continue
    os.system('mkdir -p logs/condor/merge/err/')
    os.system('rm -rf logs/condor/merge/err/*'+tag+'*'+variable+'*')
    os.system('mkdir -p logs/condor/merge/log/')
    os.system('rm -rf logs/condor/merge/log/*'+tag+'*'+variable+'*')
    os.system('mkdir -p logs/condor/merge/out/')
    os.system('rm -rf logs/condor/merge/out/*'+tag+'*'+variable+'*')
    os.environ['TAG'] = tag
    os.environ['FOLDER'] = options.folder
    os.environ['VARIABLE'] = variable
    os.environ['CLUSTER'] = options.cluster
    os.system('condor_submit merge.submit')
os.system('rm merge.submit')
