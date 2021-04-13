#!/bin/bash

#SBATCH --time=1:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=15360M   # memory per CPU core
#SBATCH -J "preproc2"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HEre
release=release9

rootDir=/fslhome/mpeter55/fsl_groups/fslg_autism_asym/compute/HBN_Project
codeDir=${rootDir}/code/preproc
textDir=${rootDir}/code/preproc
dataDir=${rootDir}/preproc/${release}
outDir=${rootDir}/preproc/${release}
readyDir=${outDir}

tweight=T1

subj=$1

#mkdir -p ${readyDir}/${subj}/anat
~/compute/research_bin/c3d/bin/c3d \
${outDir}/${subj}/anat/${subj}_${tweight}_n4.nii.gz \
-resample-mm 1x1x1mm \
-o ${outDir}/${subj}/anat/${subj}_${tweight}_resampled.nii.gz

#cp ${outDir}/${time}/${subj}/anat/${subj}_${tweight}_resampled.nii.gz ${readyDir}/${time}/${subj}/anat/${subj}_${tweight}.nii.gz
			
