#!/bin/bash


for subj in `cat /fslhome/mpeter55/fsl_groups/fslg_csf_autism/compute/HBN_Project/code/subjids/release1-7_ids.txt`; do
	sbatch \
	-o ~/logfiles/${1}/output_${subj}.txt \
	-e ~/logfiles/${1}/error_${subj}.txt \
	~/fsl_groups/fslg_csf_autism/compute/HBN_Project/code/freesurfer_code/freesurfer_stats_table2.sh \
	${subj}
	sleep 1
done
