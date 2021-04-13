#!/bin/bash

#SBATCH --time=00:30:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=15360M   # memory per CPU core
#SBATCH -J "preproc1"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
release=release9
rootDir=/fslhome/mpeter55/fsl_groups/fslg_autism_asym/compute/HBN_Project
codeDir=${rootDir}/code/preproc
textDir=${rootDir}/code/preproc
dataDir=${rootDir}/data/reformat/${release}/t1
outDir=${rootDir}/preproc
tweight=T1

subj=$1

~/compute/research_bin/antsbin/bin/N4BiasFieldCorrection \
-v \
-d 3 \
-i ${dataDir}/${subj}/anat/${subj}_${tweight}.nii.gz \
-o ${outDir}/${release}/${subj}/anat/${subj}_${tweight}_n4.nii.gz \
-s 4 \
-b [200] \
-c [50x50x50x50,0.000001]

