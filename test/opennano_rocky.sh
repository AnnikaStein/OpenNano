#!/usr/local_rwth/bin/zsh

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --time=00:05:00

# with/out gpu: 
# SBATCH --gres=gpu:1

#SBATCH --output=logs/output.%J.txt

# SBATCH --account=rwth1377

#SBATCH --mail-type=ALL
#SBATCH --mail-user=annika.stein@rwth-aachen.de

# %% Get MiniAOD Open Data file from remote storage (CONDA env, for xrdcp)
    source ~/miniconda3/bin/activate
    conda activate OpenDataTorch
    export XRD_RUNFORKHANDLER=1
    export X509_USER_PROXY=/home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/x509up_u40434
    export X509_CERT_DIR=/work/um106329/conda-storage/envs/OpenDataTorch/etc/grid-security/certificates
    voms-proxy-info
    xrdcp root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TTToSemiLeptonic_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/00000/001CCEB6-4EC4-E511-B8DC-00259074AEAC.root /home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src/PhysicsTools/OpenNano/test/tt_miniaod.root
    conda deactivate

# %% Run over Open Data in job (CMSSW env)
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src
    cmsenv
    export CMS_PATH=/home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src
    cd PhysicsTools/OpenNano/test
    cmsRun custom_tt_cfg.py

# %% Machine Learning, Data Processing (CONDA env)
    module load CUDA/11.8.0
    export PYTHONPATH=/work/um106329/conda-storage/envs/OpenDataTorch/lib/python3.10/site-packages:$PYTHONPATH
    export PATH=/work/um106329/conda-storage/envs/OpenDataTorch/bin:/home/um106329/miniconda3/condabin:$PATH
    source ~/miniconda3/bin/activate
    conda activate OpenDataTorch
   # export XRD_RUNFORKHANDLER=1
   # export X509_USER_PROXY=/home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/x509up_u40434
   # export X509_CERT_DIR=/work/um106329/conda-storage/envs/OpenDataTorch/etc/grid-security/certificates
    voms-proxy-info
    #python3 PhysicsTools/OpenNano/test/environment_test.py
    python3 process_opennano.py
