#!/usr/bin/env bash

export USER=${4}
echo "User is: ${4}"
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
echo $(hostname)
source /cvmfs/cms.cern.ch/cmsset_default.sh

if [ "${3}" == "kisti" ]; then
    env
    /usr/bin/voms-proxy-info -exists
    if [ $? -eq 0 ]; then
	echo "No need to copy"
	ls -l /tmp/x509up_u$(id -u)
	/usr/bin/voms-proxy-info -all
    else
	cp ./x509up_u* /tmp
	ls -l /tmp/x509up_u$(id -u)
	/usr/bin/voms-proxy-info -all
    fi
    xrdcp -s root://cms-xrdr.private.lo:2094//xrd/store/user/$USER/cmssw_11_3_4.tgz .
    echo "Decaf correctly copied"
    xrdcp -s root://cms-xrdr.private.lo:2094//xrd/store/user/$USER/pylocal_3_8.tgz .
    echo "Python correctly copied"
else
    xrdcp -s root://cmseos.fnal.gov//store/user/$USER/cmssw_11_3_4.tgz .
    echo "Decaf correctly copied"
    xrdcp -s root://cmseos.fnal.gov//store/user/$USER/pylocal_3_8.tgz .
    echo "Python correctly copied"
fi
tar -zxvf cmssw_11_3_4.tgz
tar -zxvf pylocal_3_8.tgz
rm cmssw_11_3_4.tgz
rm pylocal_3_8.tgz
export SCRAM_ARCH=slc7_amd64_gcc900
cd CMSSW_11_3_4/src
scramv1 b ProjectRename
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
export PYTHONPATH=${_CONDOR_SCRATCH_DIR}/site-packages:$PYTHONPATH
export PYTHONPATH=$(find ${_CONDOR_SCRATCH_DIR}/site-packages/ -name *.egg |tr '\n' ':')$PYTHONPATH
export PYTHONWARNINGS="ignore"
echo "Updated python path: " $PYTHONPATH
cd decaf/analysis
echo "python3 merge.py --folder ${1} --variable ${2}"
python3 merge.py --folder ${1} --variable ${2}
ls ${1}/${2}.merged
cp ${1}/${2}.merged ${_CONDOR_SCRATCH_DIR}/${2}.merged
