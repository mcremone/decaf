#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
#import numexpr
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset', default='')
parser.add_option('-e', '--exclude', help='exclude', dest='exclude', default='')
parser.add_option('-p', '--processor', help='processor', dest='processor', default='')
parser.add_option('-m', '--metadata', help='metadata', dest='metadata', default='')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

os.system("mkdir -p hists/"+options.processor)

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../cmssw_11_3_4.tgz '
              '--exclude=\'src/decaf/analysis/logs\' '
              '--exclude=\'src/decaf/analysis/plots\' '
              '--exclude=\'src/decaf/analysis/datacards\' '
              '--exclude=\'src/decaf/analysis/results\' '
              '--exclude=\'src/decaf/analysis/data/models\' '
              '--exclude=\'src/decaf/analysis/hists/*/*.futures\' '
              '--exclude=\'src/decaf/analysis/hists/*/*.merged\' '
              '--exclude=\'src/decaf/analysis/hists/*/*.reduced\' '
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
Executable = run.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = run.sh, run.py, /tmp/x509up_u556950957
Output = logs/condor/run/out/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/run/err/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/run/log/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(PROCESSOR)_$ENV(SAMPLE).futures=$ENV(PWD)/hists/$ENV(PROCESSOR)/$ENV(SAMPLE).futures"
Arguments = $ENV(METADATA) $ENV(SAMPLE) $ENV(PROCESSOR) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
JobBatchName = $ENV(BTCN)
request_cpus = 8
request_memory = 30000
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
Output = logs/condor/run/out/%PROCESSOR%_%SAMPLE%_$(Cluster)_$(Process).stdout
Error = logs/condor/run/err/%PROCESSOR%_%SAMPLE%_$(Cluster)_$(Process).stderr
Log = logs/condor/run/log/%PROCESSOR%_%SAMPLE%_$(Cluster)_$(Process).log
TransferOutputRemaps = "%PROCESSOR%_%SAMPLE%.futures=$ENV(PWD)/hists/%PROCESSOR%/%SAMPLE%.futures"
Arguments = %METADATA% %SAMPLE% %PROCESSOR% %CLUSTER% $ENV(USER) 
JobBatchName = %BTCN%
request_cpus = 8
request_memory = 30000
Queue 1"""

if options.cluster == 'lxplus':
    if options.copy:
        os.system('xrdcp -f ../../../../cmssw_11_3_4.tgz root://eosuser.cern.ch//eos/user/'+os.environ['USER'][0] +'/'+ os.environ['USER']+'/cmssw_11_3_4.tgz')
        os.system('xrdcp -f ../../../../pylocal_3_8.tgz root://eosuser.cern.ch//eos/user/'+os.environ['USER'][0] + '/'+os.environ['USER']+'/pylocal_3_8.tgz')
    jdl = """universe                = vanilla
executable              = run.sh
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = run.sh, /afs/cern.ch/user/$ENV(INITIAL)/$ENV(USER)/private/x509up
output                  = logs/condor/run/out/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
error                   = logs/condor/run/err/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
log                     = logs/condor/run/log/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps    = "$ENV(PROCESSOR)_$ENV(SAMPLE).futures=/eos/user/$ENV(INITIAL)/$ENV(USER)/CMSSW_11_3_4/src/decaf/analysis/hists/$ENV(PROCESSOR)/$ENV(SAMPLE).futures"
Arguments               = $ENV(METADATA) $ENV(SAMPLE) $ENV(PROCESSOR) $ENV(CLUSTER) $ENV(USER)
MY.SingularityImage     = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cat/cmssw-lxplus/cmssw-el7-lxplus:latest/"
output_destination      = root://eosuser.cern.ch//eos/user/$ENV(INITIAL)/$ENV(USER)/CMSSW_11_3_4/src/decaf/analysis/
request_cpus            = 8
request_memory          = 16000
JobBatchName            = $ENV(BTCN)
+JobFlavour             = "tomorrow"
Queue 1"""

jdl_file = open("run.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

with gzip.open('metadata/'+options.metadata+'.json.gz') as fin:
    datadef = json.load(fin)

for dataset, info in datadef.items():
    if options.dataset:
        if not any(_dataset in dataset for _dataset in options.dataset.split(',')): continue
    if options.exclude:
        if any(_dataset in dataset for _dataset in options.exclude.split(',')): continue
    os.system('mkdir -p logs/condor/run/err/')
    os.system('rm -rf logs/condor/run/err/*'+options.processor+'*'+dataset+'*')
    os.system('mkdir -p logs/condor/run/log/')
    os.system('rm -rf logs/condor/run/log/*'+options.processor+'*'+dataset+'*')
    os.system('mkdir -p logs/condor/run/out/')
    os.system('rm -rf logs/condor/run/out/*'+options.processor+'*'+dataset+'*')
    os.environ['SAMPLE'] = dataset
    os.environ['BTCN'] = dataset.split('____')[0]
    os.environ['PROCESSOR']   = options.processor
    os.environ['METADATA']   = options.metadata
    os.environ['CLUSTER'] = options.cluster
    if options.cluster == 'lxplus':
        os.environ['INITIAL'] = os.environ['USER'][0]
        os.system('mkdir -p /eos/user/'+os.environ['USER'][0] +'/'+ os.environ['USER']+'/CMSSW_11_3_4/src/decaf/analysis/logs/condor/run/err/')
        os.system('rm -rf /eos/user/'+os.environ['USER'][0] +'/'+ os.environ['USER']+'/CMSSW_11_3_4/src/decaf/analysis/logs/condor/run/err/*'+options.processor+'*'+dataset+'*')
        os.system('mkdir -p /eos/user/'+os.environ['USER'][0] +'/'+ os.environ['USER']+'/CMSSW_11_3_4/src/decaf/analysis/logs/condor/run/log/')
        os.system('rm -rf /eos/user/'+os.environ['USER'][0] +'/'+ os.environ['USER']+'/CMSSW_11_3_4/src/decaf/analysis/logs/condor/run/log/*'+options.processor+'*'+dataset+'*')
        os.system('mkdir -p /eos/user/'+os.environ['USER'][0] +'/'+ os.environ['USER']+'/CMSSW_11_3_4/src/decaf/analysis/logs/condor/run/out/')
        os.system('rm -rf /eos/user/'+os.environ['USER'][0] +'/'+ os.environ['USER']+'/CMSSW_11_3_4/src/decaf/analysis/logs/condor/run/out/*'+options.processor+'*'+dataset+'*')
    if options.cluster == 'lpc':
        jdl_dataset = jdl
        jdl_dataset = jdl_dataset.replace('%SAMPLE%', dataset)
        jdl_dataset = jdl_dataset.replace('%BTCN%', dataset.split('____')[0])
        jdl_dataset = jdl_dataset.replace('%PROCESSOR%', options.processor)
        jdl_dataset = jdl_dataset.replace('%METADATA%', options.metadata)
        jdl_dataset = jdl_dataset.replace('%CLUSTER%', options.cluster)
        jdl_file = open("run.submit", "w") 
        jdl_file.write(jdl_dataset) 
        jdl_file.close() 
    os.system('condor_submit run.submit')
os.system('rm run.submit')
