#!/usr/local_rwth/bin/zsh

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
# SBATCH --cpus-per-task=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1024M

# with/out gpu: 
# SBATCH --gres=gpu:1


# SBATCH --time=00:02:00
#SBATCH --time=00:50:00


#SBATCH --array=1-2

#SBATCH --output=logs/output-%A_%a.txt


# SBATCH --account=rwth1377


#SBATCH --mail-type=ALL
#SBATCH --mail-user=annika.stein@rwth-aachen.de



# %% Bookkeeping
    echo "Welcome to SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID and SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
    echo "# ================= Create directory in /hpcwork and get file ================= #"
    mkdir -p /hpcwork/rwth1377/OpenData/OpenNano/ttsemi_tests/custom_arrays_${SLURM_ARRAY_JOB_ID}
    INFILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src/PhysicsTools/OpenNano/samples_MiniAOD/ttsemi.txt)
    echo "Run over:" $INFILE
    echo "# ----------------- #"
    echo ""
    
# %% Get MiniAOD Open Data file from remote storage (CONDA env, for xrdcp)
    echo "# ================= Copy MiniAOD to /tmp ================= #"
    source ~/miniconda3/bin/activate
    conda activate OpenDataTorch
    export XRD_RUNFORKHANDLER=1
    export X509_USER_PROXY=/home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/x509up_u40434
    export X509_CERT_DIR=/work/um106329/conda-storage/envs/OpenDataTorch/etc/grid-security/certificates
    voms-proxy-info
    xrdcp root://eospublic.cern.ch/$INFILE $TMP/miniaod_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root
    if [[ $rc != 0 ]]
        then
            echo "xrdcp exit code: " $rc
            echo "Exiting."
            exit $rc
        fi
    conda deactivate
    echo "# ----------------- #"
    echo ""

# %% Run over Open Data in job (CMSSW env)
    echo "# ================= Run MiniAOD -> OpenNano ================= #"
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src
    cmsenv
    export CMS_PATH=/home/um106329/BMBF_AISafety/OpenDataAISafety/CMS/CMSSW_10_6_30/src
    cd PhysicsTools/OpenNano/test
    
    cmsDriver.py --python_filename $TMP/custom_cfg_job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py --eventcontent NANOAODSIM --datatier NANOAODSIM \
      --fileout file:$TMP/custom_nanoaod_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root --conditions 102X_mcRun2_asymptotic_v8 --step NANO \
      --filein file:$TMP/miniaod_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root --era Run2_25ns,run2_nanoAOD_106X2015 --mc -n -1 --nThreads 8 \
      --customise PhysicsTools/OpenNano/opennano_cff.Opennano_customizeMC_allPF_add_CustomTagger_and_Truth \
      --customise_commands "process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)));process.MessageLogger.cerr.FwkReport.reportEvery=1000;"
    if [[ $rc != 0 ]]
        then
            echo "python3 cmsDriver.py exit code: " $rc
            echo "Exiting."
            exit $rc
        fi
    echo "# ----------------- #"
    echo ""

# %% Machine Learning, Data Processing (CONDA env)
    echo "# ================= ROOT -> PARQUET ================= #"
    module load CUDA/11.8.0
    export PYTHONPATH=/work/um106329/conda-storage/envs/OpenDataTorch/lib/python3.10/site-packages:$PYTHONPATH
    export PATH=/work/um106329/conda-storage/envs/OpenDataTorch/bin:/home/um106329/miniconda3/condabin:$PATH
    source ~/miniconda3/bin/activate
    conda activate OpenDataTorch
    python3 process_opennano_parquet.py --fname_in=$TMP/custom_nanoaod_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root --fname_out_pattern=/hpcwork/rwth1377/OpenData/OpenNano/ttsemi_tests/custom_arrays_${SLURM_ARRAY_JOB_ID}/arrays_${SLURM_ARRAY_TASK_ID}
    if [[ $rc != 0 ]]
        then
            echo "python3 process_opennano_parquet.py exit code: " $rc
            echo "Exiting."
            exit $rc
        fi
    echo "# ----------------- #"
    echo ""
    echo "Success! Exiting."
