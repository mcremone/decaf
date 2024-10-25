cat $HOME/private/$USER.txt | voms-proxy-init -voms cms --valid 140:00 -pwstdin
source /cvmfs/cms.cern.ch/rucio/setup-py3.sh
export RUCIO_ACCOUNT=mcremone
export PYTHONPATH=~/.local/lib/python3.8/site-packages#:$PYTHONPATH
export PYTHONWARNINGS="ignore"
#export PATH=~/.local/bin:$PATH


