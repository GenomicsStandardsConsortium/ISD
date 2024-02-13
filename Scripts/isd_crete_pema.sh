#!/bin/bash -l

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=100G
#SBATCH --job-name="pema215_true_fastp"
#SBATCH --mail-user=s.paragkamian@hcmr.gr
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --requeue

start=`date +%s`

module purge # unloads all previous loads

module load singularity/3.7.1 #loads singularity

singularity run -B /home1/s.paragkamian/isd-crete/pema215_test2/:/mnt/analysis /home1/s.paragkamian/pema_v2.1.5.sif

### for checkpoints
#singularity exec \
#    -B /home1/s.paragkamian/isd-crete/pema_asv/:/mnt/analysis \
#    /mnt/big/containers/singularity/pema_v.2.1.4.sif \
#    /home/tools/BDS/.bds/bds -r /mnt/analysis/clustering.chp

module unload singularity/3.7.1 #unloads singularity

end=`date +%s`
runtime=$((end-start))
echo "Job ID: " $SLURM_JOB_ID
echo "Job name: " $SLURM_JOB_NAME
echo $runtime "in seconds" 
echo $((runtime/60)) "in minutes" 
