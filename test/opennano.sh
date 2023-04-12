#!/usr/local_rwth/bin/zsh

#SBATCH --ntasks=1

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=1

#SBATCH --output=output.%J.txt

#SBATCH --account=rwth1244

#SBATCH --time=00:01:00

# with gpu: 

#SBATCH --gres=gpu:1

#SBATCH --mail-type=ALL

#SBATCH --mail-user=annika.stein@rwth-aachen.de

cd /home/um106329/analysis/VHccAllIn/explore
module unload intelmpi; module switch intel gcc
module load cuda/11.6
module load cudnn
source ~/miniconda3/bin/activate
conda activate OpenDataTorch
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/home/um106329/analysis/VHccAllIn/explore/x509up_u40434
export X509_CERT_DIR=/work/um106329/conda-storage/envs/torch2Ana/etc/grid-security/certificates
voms-proxy-info
python3 datastructure.py