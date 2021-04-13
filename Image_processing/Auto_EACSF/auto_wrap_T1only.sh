#!/bin/bash

release=release1-7

rootDir=/fslhome/mpeter55/fsl_groups/fslg_csf_autism/compute/HBN_Project
scriptDir=${rootDir}/code/auto_code
subjidsDir=${rootDir}/code/subjids
var=`date +"%Y%m%d-%H%M%S"`

for i in `cat ${subjidsDir}/${release}_ids.txt`;do
	rm -r ${rootDir}/output/${release}/${i}
	mkdir -p /fslhome/mpeter55/logfiles/${var}
	sbatch \
	-o /fslhome/mpeter55/logfiles/$var/${i}_output_log.txt \
	-e /fslhome/mpeter55/logfiles/$var/${i}_error_log.txt \
	${scriptDir}/auto_job_T1only.sh ${i}

	sleep 1
done
