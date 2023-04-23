# CMS Open Data for ML

This shows an exemplary way to read CMS Open Data MiniAOD, and end up with parquet files containing the data in the form of awkward arrays. This is a storage-friendly alternative to just saving in numpy / pandas etc. which would already fill empty or non-existing features (varying number of particles per event...) by padding with some placeholder value to create flat tuples.

Running on OpenData requires at least some CMS-software release to read and process the samples, packed in a virtual machine or from /cvmfs, if you have access. Then convert them to easier-to-process root-files, and make them usable by any deep learning framework.

## Overview
### CMS-part
The magic that needs to be done first:

```shell
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_10_6_30
cd CMSSW_10_6_30/src/
cmsenv
git cms-init
git cms-merge-topic 39040
git clone -b opendata git@github.com:AnnikaStein/OpenNano.git PhysicsTools/OpenNano
scram b -j 18
```

Using cmsRun on a "non-official" site though means a local site-configuration is necessary. Let's trick `cmsRun` to think we are an actual site, despite working locally (getting lucky is not a crime, there appears to be a config for RWTH-HPC already):  
```shell
mkdir -p /home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src/SITECONF/local/JobConfig
cp /cvmfs/cms.cern.ch/SITECONF/T2_DE_RWTH/RWTH-HPC/JobConfig/site-local-config.xml /home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src/SITECONF/local/JobConfig
export CMS_PATH=/home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src
scram b -j 18
```
(kudos to https://twiki.cern.ch/twiki/bin/view/Main/RobinGitlabCICMSSW which is a similar use case)

### Conda-Setup
No conda setup yet? Then do
```shell
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

We will need to use some conda environment containing relevant packages.
```shell
conda env create -f env.yml
```

With the help of a conda environment, create a proxy (e.g. `voms-proxy-init --voms cms --vomses .grid-security/vomses --valid=192:00`) and copy the proxy into some accessible directory (`cp /tmp/x509up_u40434 /home/um106329/BMBF_AISafety/OpenDataAISafety/CMS`)  
Problems with that? Do
```shell
mkdir ~/.grid-security
cp -r /cvmfs/grid.cern.ch/etc/grid-security/vomses ~/.grid-security
cp -r /cvmfs/grid.cern.ch/etc/grid-security/vomsdir ~/.grid-security/
```
Now assume there already is a MiniAOD file (after doing `xrdcp`, or if you used the https-option of the opendata-client).  
Example: `xrdcp root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TTToSemiLeptonic_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/00000/001CCEB6-4EC4-E511-B8DC-00259074AEAC.root /home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src/PhysicsTools/OpenNano/test/tt_miniaod.root`
Then one can finally do:
```shell
cmsDriver.py --python_filename custom_tt_cfg.py --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --fileout file:custom_tt_nanoaod.root --conditions 102X_mcRun2_asymptotic_v8 --step NANO \
  --filein file:tt_miniaod.root --era Run2_25ns,run2_nanoAOD_106X2015 --no_exec --mc -n 100 \
  --customise PhysicsTools/OpenNano/opennano_cff.Opennano_customizeMC_allPF_add_CustomTagger_and_Truth
```

and voilà: `cmsRun custom_tt_cfg.py` finally works, also inside a SLURM job.

## Processing OpenNano
The output of the previous step (a NanoAOD-like file) can be further processed, examples provided in `test/process_opennano.py`.

### Do all steps (after initial setup) via SLURM
It's convenient to have a set of files to run over, stored in a .txt file. Use something similar to that:
`xrdfs eospublic.cern.ch ls /eos/opendata/cms/mc/RunIIFall15MiniAODv2/TTToSemiLeptonic_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/00000 > samples_MiniAOD/ttsemi.txt`  
Some hints: the running number 00000 is not necessarily the only one, look for the next ones in line like 00001... and if you want, concatenate them together into a super-txt file containing more than 1K lines / filepaths.


Now one would like to do all those steps without manually pasting such commands for every file. That's what `sbatch test/opennano_rocky_multiFile.sh` is for. A set of filelists can be stored in `samples_MiniAOD` for that purpose. And because OpenData does not require authentication as CMS member, this will even work without the active proxy.