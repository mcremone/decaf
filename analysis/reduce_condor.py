#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
import os
from optparse import OptionParser
import numpy as np
from coffea.util import load

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset', default='')
parser.add_option('-e', '--exclude', help='exclude', dest='exclude', default='')
parser.add_option('-f', '--folder', help='folder', dest='folder')
parser.add_option('-v', '--variable', help='variable', dest='variable')
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
              '--exclude=\'src/decaf/analysis/hists/*/*.reduced\' '
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
        jdl = """universe = container
container_image = /cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el7:latest
+SingularityBind = "/cvmfs,/cms,/cms_scratch"
Executable = reduce.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = reduce.sh, /tmp/x509up_u556950957
Output = logs/condor/reduce/out/$ENV(TAG)_$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/reduce/err/$ENV(TAG)_$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/reduce/log/$ENV(TAG)_$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE)_$ENV(SAMPLE).reduced=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE)--$ENV(SAMPLE).reduced"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(SAMPLE) $ENV(CLUSTER) $ENV(USER)
JobBatchName = $ENV(VARIABLE)
accounting_group=group_cms
request_cpus = 16
request_disk = 10G
Queue 1"""

if options.cluster == 'lpc':
    if options.copy:
        os.system('xrdcp -f ../../../../cmssw_11_3_4.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/cmssw_11_3_4.tgz')
        os.system('xrdcp -f ../../../../pylocal_3_8.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal_3_8.tgz')
    jdl = """universe = vanilla
Executable = run.sh
+ApptainerImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7"
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = run.sh
Output = logs/condor/reduce/out/%TAG%_%SAMPLE%_%VARIABLE%_$(Cluster)_$(Process).stdout
Error = logs/condor/reduce/err/%TAG%_%SAMPLE%_%VARIABLE%_$(Cluster)_$(Process).stderr
Log = logs/condor/reduce/log/%TAG%_%SAMPLE%_%VARIABLE%_$(Cluster)_$(Process).log
TransferOutputRemaps = "%VARIABLE%_%SAMPLE%.reduced=$ENV(PWD)/%FOLDER%/%VARIABLE%--%SAMPLE%.reduced"
Arguments = %FOLDER% %VARIABLE% %SAMPLE% %CLUSTER% $ENV(USER)
request_cpus = 16
request_disk = 10G
request_memory = 6000
Queue 1"""

if options.cluster == 'lxplus':
    if options.copy:
        os.system('xrdcp -f ../../../../cmssw_11_3_4.tgz root://eosuser.cern.ch//eos/user/'+os.environ['USER'][0] +'/'+ os.environ['USER']+'/cmssw_11_3_4.tgz')
        os.system('xrdcp -f ../../../../pylocal_3_8.tgz root://eosuser.cern.ch//eos/user/'+os.environ['USER'][0] + '/'+os.environ['USER']+'/pylocal_3_8.tgz')
    jdl = """universe                = vanilla
executable              = run.sh
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = reduce.sh, /afs/cern.ch/user/m/mcremone/private/x509up
output                  = logs/condor/reduce/out/$ENV(TAG)_$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
error                   = logs/condor/reduce/err/$ENV(TAG)_$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
log                     = logs/condor/reduce/log/$ENV(TAG)_$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps    = "$ENV(VARIABLE)_$ENV(SAMPLE).reduced=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE)--$ENV(SAMPLE).reduced"
Arguments               = $ENV(FOLDER) $ENV(VARIABLE) $ENV(SAMPLE) $ENV(CLUSTER) $ENV(USER)
MY.SingularityImage     = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cat/cmssw-lxplus/cmssw-el7-lxplus:latest/"
request_cpus            = 16
request_memory          = 6000
JobBatchName            = $ENV(VARIABLE)
+JobFlavour             = "tomorrow"
Queue 1"""

jdl_file = open("reduce.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

pd = []
futurefile=''
for filename in os.listdir(options.folder):
    if '.futures' not in filename: continue
    futurefile=filename
    if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])

tag=options.folder.split('/')[-1]
variables=load(options.folder+'/'+futurefile).keys()
for pdi in pd:
    if options.dataset:
        if not any(_dataset in pdi for _dataset in options.dataset.split(',')): continue
    if options.exclude:
        if any(_dataset in pdi for _dataset in options.exclude.split(',')): continue
    for variable in variables:
        if options.variable and options.variable != variable: continue
        os.system('mkdir -p logs/condor/reduce/err/')
        os.system('rm -rf logs/condor/reduce/err/*'+tag+'*'+pdi+'*'+variable+'*')
        os.system('mkdir -p logs/condor/reduce/log/')
        os.system('rm -rf logs/condor/reduce/log/*'+tag+'*'+pdi+'*'+variable+'*')
        os.system('mkdir -p logs/condor/reduce/out/')
        os.system('rm -rf logs/condor/reduce/out/*'+tag+'*'+pdi+'*'+variable+'*')
        os.environ['TAG'] = tag
        os.environ['FOLDER'] = options.folder
        os.environ['SAMPLE'] = pdi
        os.environ['VARIABLE'] = variable
        os.environ['CLUSTER'] = options.cluster
        if options.cluster == 'lpc':
            jdl_dataset = jdl
            jdl_dataset = jdl_dataset.replace('%TAG%', tag)
            jdl_dataset = jdl_dataset.replace('%FOLDER%', options.folder)
            jdl_dataset = jdl_dataset.replace('%SAMPLE%', pdi)
            jdl_dataset = jdl_dataset.replace('%VARIABLE%', variable)
            jdl_dataset = jdl_dataset.replace('%CLUSTER%', options.cluster)
            jdl_file = open("run.submit", "w") 
            jdl_file.write(jdl_dataset) 
            jdl_file.close() 
        os.system('condor_submit reduce.submit')
os.system('rm reduce.submit')
