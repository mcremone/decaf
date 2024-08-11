Readme 




<img src="https://user-images.githubusercontent.com/10731328/193421563-cf992d8b-8e5e-4530-9179-7dbd507d2e02.png" width="350"/>

# **D**ark matter **E**xperience with the **C**offea **A**nalysis **F**ramework

## Initial Setup

(Note: remember to replace `<USERNAME>` with *your* username if copying commands below)

---

### LPC Setup

First, log into an LPC node:

```
kinit <USERNAME>@FNAL.GOV
ssh <USERNAME>@cmslpc-el9.fnal.gov
```

The CMSSW version used runs on slc7. You'll need to setup the correct OS environment using [singularity](https://cms-sw.github.io/singularity.html) (more LPC documentation [here](https://uscms.org/uscms_at_work/computing/setup/setup_software.shtml#apptainer)):

```
cmssw-el7 -p --bind `readlink $HOME` --bind `readlink -f ${HOME}/nobackup/` --bind /uscms_data --bind /cvmfs -- /bin/bash -l
```

Install `CMSSW_11_3_4` in your `nobackup` area:

```
cd ~/nobackup

source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
```

---

### KISTI Setup

First, log into KISTI:

```
ssh -p 4280 <USERNAME>@ui20.sdfarm.kr
```

The CMSSW version used runs on slc7. You'll need to setup the correct OS environment using [singularity](https://cms-sw.github.io/singularity.html). On KISTI, this can be done with:

```
setup_el7
```

Install CMSSW_11_3_4 on KISTI:

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc900 # this is set automatically on LPC but not on KISTI
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
```

---

### Installing Packages

The rest of the setup should be the same regardless of what cluster you are working on.

Install `combine` ([see detailed instructions](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#installation-instructions)):

```
cd $CMSSW_BASE/src
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.1.0 # current recommeneded tag (Jan 2024)
scramv1 b clean; scramv1 b # always make a clean build
```

Fork this repo on github and clone it into your `CMSSW_11_3_4/src` directory:

```
cd $CMSSW_BASE/src
git clone https://github.com/<USERNAME>/decaf.git -b UL
cd decaf
```

Then, setup the proper dependences:

```
source setup.sh
```

This script installs the necessary packages as user packages (Note: Pip gave errors when running `setup.sh` for the first time, but it seemed to install everything just fine. No errors showed up when running `setup.sh` a second time.). This is a one-time setup.

## Setup when returning

When you log in after doing all the one time installations, you only need to set up the environments. Consider using aliases, scripts, etc. to make your life easier.

Consider adding this line to your ~/.bashrc or ~/.bash_profile if you haven't done so already rather than sourcing every time:
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
```

Singularity on LPC:

```
cmssw-el7 -p --bind `readlink $HOME` --bind `readlink -f ${HOME}/nobackup/` --bind /uscms_data --bind /cvmfs -- /bin/bash -l
```

Singularity on KISTI:

```
setup_el7
```

Then, go to where you installed CMSSW and do:

```
cd CMSSW_11_3_4/src
cmsenv
cd decaf
source env.sh
```

By running this script you will also initialize your grid certificate (Note: `setup.sh` also runs `env.sh`). This requires you to save your grid certificate password in `$HOME/private/$USER.txt`. Alternatively, you can change this location, or you can comment this out and initialize it manually every time.


## Listing Input Files

The list of input files for the analyzer can be generated as a JSON file using the `macros/list.py` script. This script will run over the datasets listed in `data/process.py`, find the list of files for each dataset, “pack” them into small groups for condor jobs, and output the list of groups as a JSON file in `metadata/`.

The options for this script are:


1. **`-d` or `--dataset`**:
   - Select a specific dataset to pack. By default, it will run over all datasets in `process.py`.
   - **Usage**: `-d <dataset_name>`

2. **`-y` or `--year`**:
   - Data year. Options are `2016pre`, `2016post`, `2017`, and `2018`.
   - **Usage**: `-y <year>`

3. **`-m` or `--metadata`**:
   - Name of metadata output file. Output will be saved in `metadata/<NAME>.json`.
   - **Usage**: `-m <name>`

4. **`-p` or `--pack`**:
   - Size of file groups. The smaller the number, the more condor jobs will run. The larger the number, the longer each condor job will take. We tend to pick 32, but the decision is mostly arbitrary.
   - **Usage**: `-p <size>`

5. **`-s` or `--special`**:
   - Size of file groups for special datasets. For a specific dataset, use a different size with respect to the one established with `--pack`. The syntax is `-s <DATASET>:<NUMBER>`.
   - **Usage**: `-s <dataset>:<number>`

6. **`-c` or `--custom`**:
   - Boolean to decide to use public central NanoAODs (if `False`) or private custom NanoAODs (if `True`). Default is `False`.
   - **Usage**: `-c` (no argument needed)


As an example, to generate the JSON file for all 2017 data:

```
python3 macros/list.py -y 2017 -m 2017 -p 32
```

As a reminder, this script assumes that you are in the `decaf/analysis` directory when running. The output above will be saved in `metadata/2017.json`.

If using the `--custom` option, the script can take several hours to run, so it’s best to use a process manager such as `nohup` or `tmux` to avoid the program crashing in case of a lost connection. For example

```
nohup python3 macros/list.py -y 2017 -m 2017 -p 32 -c &
```

The `&` option at the end of the command lets it run in the background, and the std output and error is saved in `nohup.out`. 

The `nohup` command is useful and recommended for running most scripts, but you may also use tools like `tmux` or `screen`.


## Computing MC b-Tagging Efficiencies

MC b-tagging efficiencies are needed by most of the analyses to compute the b-tag event weight, once such efficiencies are corrected with the POG-provided b-tag SFs. To compute them, we first need to run the `common` module in `util`:

```
python utils/common.py
```

This will generate a series of auxiliary functions and information, like the AK4 b-tagging working points, and it will save such information in a `.coffea` file in the `data` folder. AK4 b-tagging working points are essential to measure the MC efficiencies and they are used by the `btag` processor in the `processors` folder. To generate the processor file: 

```
python3 processors/btageff.py -y 2018 -m 2018 -n 2018
```

The options for this script are:

1. **`-y` or `--year`**:
   - Data year. Options are `2016pre`, `2016post`, `2017`, and `2018`.
   - **Usage**: `-y <year>`

2. **`-m` or `--metadata`**:
   - Metadata file to be used as input.
   - **Usage**: `-m <metadata_file>`

3. **`-n` or `--name`**:
   - Name of the output processor file. In this case, it will generate a file called `btageff2018.processor` stored in the `data` folder.
   - **Usage**: `-n <output_name>`

To run the processor:

```
python3 run.py -p btageff2018 -m 2018 -d QCD
```

With this command you will run the `btag2018` processor over QCD MC datasets as defined by the `2018` metadata file. You will see a printout like:

Processing: QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8____4_
  Preprocessing 100% ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 32/32 [ 0:01:28 < 0:00:00 | ?   file/s ]
Merging (local) 100% ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 31/31 [ 0:00:23 < 0:00:00 | ? merges/s ]

This means an output file with histograms as defined in the btag processor file has been generated. In this case a folder called `btageff2018` inside the `hists` folder has been created. Inside this folder you can see a file called `QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8____4_.futures`, that stores the histograms. To take advantage of the parallelism offered by the HTCondor job scheduler, the `run_condor.py` script can be used:

```
python3 run_condor.py -p btag2018 -m 2018 -d QCD -c kisti -t -x
```

The options for this script are the same as for run.py, with the addition of:

1. **`-c` or `--cluster`**:
   - Specifies which cluster you are using. Currently supports `lpc` or `kisti`.
   - **Usage**: `-c <cluster_name>`

2. **`-t` or `--tar`**:
   - Tars the local python environment and the local CMSSW folder.
   - **Usage**: `-t` (no argument needed)

3. **`-x` or `--copy`**:
   - Copies these two tarballs to your EOS area. For example, to run the same setup but for a different year you won’t need to tar and copy again. You can simply do: `python run_condor.py -p btag2017 -m 2017 -d QCD -c kisti`
   - **Usage**: `-x` (no argument needed))

You can check the status of your HTCondor jobs by doing:

```
condor_q <YOUR_USERNAME>
```

After obtaining all the histograms, a first step of data reduction is needed. This step is achieved by running the `reduce.py` script:

```
python reduce.py -f hists/btag2018
```

The options of this script are:

1. **`-f` or `--folder`**:
   - Specifies the folder to be processed.
   - **Usage**: `-f <folder_name>`

2. **`-d` or `--dataset`**:
   - Specifies the dataset(s) to be processed. If not provided, defaults to `None`.
   - **Usage**: `-d <dataset_name>`

3. **`-e` or `--exclude`**:
   - Specifies the dataset(s) to be excluded from processing. Defaults to `None`.
   - **Usage**: `-e <dataset_name>`

4. **`-v` or `--variable`**:
   - Specifies the variable(s) to be processed. If not provided, defaults to `None`.
   - **Usage**: `-v <variable_name>`
  

All the different datasets produced at the previous step will be reduced. A different file for each variable for each reduced dataset will be produced. For example, the command above will produce the following reduced files:

```
hists/btageff2018/deepcsv--QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_120to170_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_15to30_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_170to300_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_300to470_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_30to50_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_470to600_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_50to80_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_600to800_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_800to1000_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_80to120_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_120to170_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_15to30_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_170to300_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_300to470_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_30to50_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_470to600_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_50to80_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_600to800_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_800to1000_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_80to120_TuneCP5_13TeV_pythia8.reduced
```

This step can be run in HTCondor by using the `reduce_condor.py` script. The `reduce_condor.py` script has the same options of `reduce.py`, with addition of the same `--cluster`, `--tar`, and `--copy` options descibed above when discussing `run_condor.py`.

A second step of data reduction is needed to merge all the `.reduced` files corresponding to a single variable. This is achieved by using the `merge.py` script:

```
python3 merge.py -f hists/btageff2018
```

The options of this script are:

1. **`-f` or `--folder`**:
   - Specifies the folder to be processed.
   - **Usage**: `-f <folder_name>`

2. **`-v` or `--variable`**:
   - Specifies the variable(s) to be processed. If not provided, defaults to `None`, meaning all variables are processed.
   - **Usage**: `-v <variable_name>`

3. **`-e` or `--exclude`**:
   - Specifies the variable(s) to be excluded from processing. Defaults to `None`.
   - **Usage**: `-e <variable_name>`

4. **`-p` or `--postprocess`**:
   - If specified, performs postprocessing on the merged files.
   - **Usage**: `-p` (no argument needed)

This command will produce the following files:

```
hists/btageff2018/deepcsv.merged  hists/btageff2018/deepflav.merged
```

The same script can be used to merge the the files corresponding to each single variable into a single file, using the `-p` or `—postprocess` option:

```
python3 merge.py -f hists/btageff2018 -p
```

Also this step can be run in HTCondor by using the `merge_condor.py` script. The `merge_condor.py` script has the same options of `merge.py`, with addition of the same `--cluster`, `--tar`, and `--copy` options descibed above when discussing `run_condor.py`.

---

This README is a work in progress
