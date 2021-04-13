#!/bin/bash

GenDir=/fslhome/mpeter55/fsl_groups/fslg_csf_autism/compute/HBN_Project
ResultsDir=${GenDir}/output
time=release1-7

for i in `cat ${GenDir}/code/subjids/${time}_ids.txt`; do
	cmtk statistics ${ResultsDir}/${time}/${i}/FinalMasking/${i}*_MID02.nrrd > ${ResultsDir}/${time}/${i}/${i}_FinalStats.txt
	echo Calculating statistics for ${i}
done

