#!/bin/bash

release=release9
rootDir=/fslhome/mpeter55/fsl_groups/fslg_autism_asym/compute/HBN_Project
codeDir=${rootDir}/code/preproc
textDir=${rootDir}/code/preproc
dataDir=${rootDir}/data/reformat/${release}/t1
outDir=${rootDir}/preproc/${release}
var=`date +"%Y%m%d-%H%M%S"`

for subj in `cat ${textDir}/${release}_ids.txt`; do
		mkdir -p /fslhome/mpeter55/logfiles/${var}/${subj}
		mkdir -p ${outDir}/${subj}/anat
		sbatch \
		-o /fslhome/mpeter55/logfiles/${var}/${subj}/${subj}_output_log.txt \
		-e /fslhome/mpeter55/logfiles/${var}/${subj}/${subj}_error_log.txt \
		${codeDir}/preproc_parallel_batch.sh ${subj}
		
		sleep 1
done

