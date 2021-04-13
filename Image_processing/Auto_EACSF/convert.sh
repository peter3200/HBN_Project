#!/bin/bash

subjDir=/fslhome/mpeter55/fsl_groups/fslg_csf_autism/compute/HBN_Project/output/release1-7

for i in `cat /fslhome/mpeter55/fsl_groups/fslg_csf_autism/compute/HBN_Project/code/subjids/release1-7_ids.txt`; do
	c3d ${subjDir}/${i}/FinalMasking/${i}_T1_stx_stripped_EMS_withoutVent_MID02.nrrd -o ${subjDir}/${i}/${i}_MID02.nii.gz
done

