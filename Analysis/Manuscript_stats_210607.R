#----------------------------------------------------SETUP--------------------------------------------------
# Load libraries (install them first if necessary using install.packages("packagename") )
library(mosaic)
library(dplyr)
library(ggpubr)
library(broom)
library(tidyverse)
library(rstatix)
library(emmeans)
library(psych) #This is for the QC ICC analysis
library(readxl)
library(lsr) #This is for eta-squared
library(svglite) #This is for .svg figures
library(afex) #For raincloud plots


# Load files (requires Box Drive if using local machine) -- change paths
#Alternatively, first download files from Box and then point paths to your Downloads folder
ConsensusDx <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/phenotypic_data/CMI_relevant_phenotypic_data/9994_ConsensusDx_20200630.csv")
Basic_Demos <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/phenotypic_data/CMI_relevant_phenotypic_data/9994_Basic_Demos_20200630.csv")
release1.7_ids_T1 <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/code/subjids/release1-7_ids_T1.csv")
EACSF <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/final_results/EaCSF_results_1.7.7_T1only_release1-7_1611851892.csv")
#QC <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/QC/HBN_QC/HBN_QC.csv")
FS <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/freesurfer_results/aseg_stats_HBN.csv")
Sleep <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/phenotypic_data/assessment_data/9994_SDS_20200630.csv")
SRS <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/phenotypic_data/assessment_data/9994_SRS_20200630.csv")
WISC <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/phenotypic_data/assessment_data/9994_WISC_20200630.csv")
WAIS <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/phenotypic_data/assessment_data/9994_WAIS_Abbreviated_20200630.csv")
QC2 <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/QC/HBN_QC/resolve_discrepancies/HBN_QC_Resolve_Discrepancies2.csv")
CT1 <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/freesurfer_results/aparc_RH_stats_HBN.csv")
CT2 <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/freesurfer_results/aparc_LH_stats_HBN.csv")
available <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/phenotypic_data/available_pheno/all_available_pheno.csv")
#WASI <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/phenotypic_data/assessment_data/9994_WASI_20200630.csv")


#Data organization for each dataset
#Clean up available
available <- available[,c("EID", "Sex", "Age")]

#Clean up Basic_Demos
names(Basic_Demos)[1] <- "ID"
Basic_Demos <- Basic_Demos[,c("ID", "EID", "Study.Site")]
Basic_Demos <- subset(Basic_Demos, ID!="ID")

#Clean up ConsensusDx
names(ConsensusDx)[1] <- "ID"
ConsensusDx <- ConsensusDx[,c("ID", "EID", "DX_01", "DX_02", "DX_03", "DX_04", "DX_05", "DX_06", "DX_07", "DX_08", "DX_09", "DX_10", "Site")]
ConsensusDx <- subset(ConsensusDx, ID!="ID")

#Clean up release1.7_ids_T1
release1.7_ids_T1$EID <- gsub("^.{0,4}", "", release1.7_ids_T1$MRI_ID)
names(release1.7_ids_T1)[2] <- "T1"
release1.7_ids_T1 <- release1.7_ids_T1[,c("EID", "T1")]

#Clean up EACSF
EACSF$EID <- gsub("^.{0,4}", "", EACSF$EID)
EACSF$EID <- gsub(".{1}$", "", EACSF$EID)

#Clean up CT datasets
CT1 <- CT1[,c("EID", "rh_MeanThickness_thickness")]
CT2 <- CT2[,c("EID", "lh_MeanThickness_thickness")]
CT1$EID <- gsub("^.{0,4}", "", CT1$EID)
CT2$EID <- gsub("^.{0,4}", "", CT2$EID)

#Clean up FS
FS$EID <- gsub("^.{0,4}", "", FS$EID)
FS <- FS[,c("EID", "BrainSegVol", "EstimatedTotalIntraCranialVol")]

#Clean up Sleep
Sleep <- Sleep[,c("EID", "SDS_Total_Raw")]
Sleep<-subset(Sleep, EID!="SDSCMI_001")

#Clean up SRS
SRS <- subset(SRS, EID!="SRSHBN_001")
SRS$duplicate <- !duplicated(SRS$EID)
SRS <- subset(SRS, duplicate==TRUE)
SRS <- SRS[,c("EID", "SRS_Total")]

#Clean up WISC
WISC <- WISC[,c("EID", "WISC_FSIQ")]
WISC <- subset(WISC, EID!="WISC1_001")


#Merge datasets
Comb1 <- merge(ConsensusDx, Basic_Demos, by =c("EID"), all=TRUE)
Comb2 <- merge(Comb1, release1.7_ids_T1, by =c("EID"), all=TRUE)
Comb3 <- merge(Comb2, EACSF, by =c("EID"), all=TRUE)
Comb4 <- merge(Comb3, CT1, by =c("EID"), all=TRUE)
Comb5 <- merge(Comb4, CT2, by =c("EID"), all=TRUE)
Comb6 <- merge(Comb5, FS, by =c("EID"), all=TRUE)
Comb7 <- merge(Comb6, Sleep, by =c("EID"), all=TRUE)
Comb8 <- merge(Comb7, SRS, by =c("EID"), all=TRUE)
Comb9 <- merge(Comb8, WISC, by =c("EID"), all=TRUE)
Comb10 <- merge(Comb9, QC2, by =c("EID"), all=TRUE)
Comb11 <- merge(Comb10, available, by =c("EID"), all=TRUE)
Pheno_Comb2 <- Comb11




#Exclusion criteria

  #Drop duplicates
  Pheno_Comb2$duplicate <- !duplicated(Pheno_Comb2$EID)
  Pheno_Comb2 <- subset(Pheno_Comb2, duplicate==TRUE)

  #Drop subjects without a T1 image
  Pheno_Comb2 <- subset(Pheno_Comb2, T1!="NA")

  #Drop subjects with QC segmentation rating of 2 and higher
  Pheno_Comb2<-subset(Pheno_Comb2, Avg_Rating<"2")

  #Drop cases with a T1w rating of C3
    #Convert to numerical values (run trim lines first!!)
    Pheno_Comb2$T1w_QC_Anne <- gsub("^.{0,1}", "", Pheno_Comb2$T1w_QC_Anne)
    Pheno_Comb2$T1w_QC_Matt <- gsub("^.{0,1}", "", Pheno_Comb2$T1w_QC_Matt)
    Pheno_Comb2$T1w_QC_Christopher <- gsub("^.{0,1}", "", Pheno_Comb2$T1w_QC_Christopher)
    Pheno_Comb2$T1w_QC_Maddy <- gsub("^.{0,1}", "", Pheno_Comb2$T1w_QC_Maddy)

    #Avg across rows
    Pheno_Comb2$T1w_QC_Matt[ Pheno_Comb2$T1w_QC_Matt == ""] <- NA
    Pheno_Comb2$T1w_QC_Anne[ Pheno_Comb2$T1w_QC_Anne == ""] <- NA
    Pheno_Comb2$T1w_QC_Christopher[ Pheno_Comb2$T1w_QC_Christopher == ""] <- NA
    Pheno_Comb2$T1w_QC_Maddy[ Pheno_Comb2$T1w_QC_Maddy == ""] <- NA

    qc_t1 <- data.frame(Pheno_Comb2$T1w_QC_Anne, Pheno_Comb2$T1w_QC_Matt, Pheno_Comb2$T1w_QC_Christopher, Pheno_Comb2$T1w_QC_Maddy)
    qc_t1$Pheno_Comb2.T1w_QC_Anne <-as.numeric(qc_t1$Pheno_Comb2.T1w_QC_Anne)
    qc_t1$Pheno_Comb2.T1w_QC_Matt <-as.numeric(qc_t1$Pheno_Comb2.T1w_QC_Matt)
    qc_t1$Pheno_Comb2.T1w_QC_Christopher <-as.numeric(qc_t1$Pheno_Comb2.T1w_QC_Christopher)
    qc_t1$Pheno_Comb2.T1w_QC_Maddy <-as.numeric(qc_t1$Pheno_Comb2.T1w_QC_Maddy)

    Pheno_Comb2$T1w_Avg_Rating <- rowMeans(qc_t1, na.rm=TRUE)
    Pheno_Comb2<-subset(Pheno_Comb2, T1w_Avg_Rating!="3")

  #Drop subjects without age (N=1)
  Pheno_Comb2 <- subset(Pheno_Comb2, Age!="NA")

  #Drop subjects without sex (N=0)
  Pheno_Comb2 <- subset(Pheno_Comb2, Sex!="NA")

  #Drop subjects without eTIV (N=2, SUBJID = NDARRY268AF2,NDARWC737UJ0)
  missing <- subset(Pheno_Comb2, is.na(Pheno_Comb2$BrainSegVol))
  Pheno_Comb2 <- subset(Pheno_Comb2, BrainSegVol!="NA")

#--------------------------------DX NUMBERS-------------------------------------------

#We are including Dx groups with 50+ scans available in addition to ASD and NoDx

#Create dummy var for each group that had 50+ scans, then tally each column
Pheno_Comb2$ASD_dummy <- ifelse(Pheno_Comb2$DX_01 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_02 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_03 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_04 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_05 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_06 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_07 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_08 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_09 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_10 == "Autism Spectrum Disorder", "1", "0")
#Create var for no dx
Pheno_Comb2$NoDX_dummy <- ifelse(Pheno_Comb2$DX_01 == "No Diagnosis Given" | Pheno_Comb2$DX_02 == "No Diagnosis Given" | Pheno_Comb2$DX_03 == "No Diagnosis Given" | Pheno_Comb2$DX_04 == "No Diagnosis Given" | Pheno_Comb2$DX_05 == "No Diagnosis Given" | Pheno_Comb2$DX_06 == "No Diagnosis Given" | Pheno_Comb2$DX_07 == "No Diagnosis Given" | Pheno_Comb2$DX_08 == "No Diagnosis Given" | Pheno_Comb2$DX_09 == "No Diagnosis Given" | Pheno_Comb2$DX_10 == "No Diagnosis Given", "1", "0")
#Create var for other Dx
Pheno_Comb2$ADHDCombinedType_dummy <- ifelse(Pheno_Comb2$DX_01 == "ADHD-Combined Type" | Pheno_Comb2$DX_02 == "ADHD-Combined Type" | Pheno_Comb2$DX_03 == "ADHD-Combined Type" | Pheno_Comb2$DX_04 == "ADHD-Combined Type" | Pheno_Comb2$DX_05 == "ADHD-Combined Type" | Pheno_Comb2$DX_06 == "ADHD-Combined Type" | Pheno_Comb2$DX_07 == "ADHD-Combined Type" | Pheno_Comb2$DX_08 == "ADHD-Combined Type" | Pheno_Comb2$DX_09 == "ADHD-Combined Type" | Pheno_Comb2$DX_10 == "ADHD-Combined Type", "1", "0")
Pheno_Comb2$ADHDHyperactiveImpulsiveType_dummy <- ifelse(Pheno_Comb2$DX_01 == "ADHD-Hyperactive/Impulsive Type" | Pheno_Comb2$DX_02 == "ADHD-Hyperactive/Impulsive Type" | Pheno_Comb2$DX_03 == "ADHD-Hyperactive/Impulsive Type" | Pheno_Comb2$DX_04 == "ADHD-Hyperactive/Impulsive Type" | Pheno_Comb2$DX_05 == "ADHD-Hyperactive/Impulsive Type" | Pheno_Comb2$DX_06 == "ADHD-Hyperactive/Impulsive Type" | Pheno_Comb2$DX_07 == "ADHD-Hyperactive/Impulsive Type" | Pheno_Comb2$DX_08 == "ADHD-Hyperactive/Impulsive Type" | Pheno_Comb2$DX_09 == "ADHD-Hyperactive/Impulsive Type" | Pheno_Comb2$DX_10 == "ADHD-Hyperactive/Impulsive Type", "1", "0")
Pheno_Comb2$ADHDInattentiveType_dummy <- ifelse(Pheno_Comb2$DX_01 == "ADHD-Inattentive Type" | Pheno_Comb2$DX_02 == "ADHD-Inattentive Type" | Pheno_Comb2$DX_03 == "ADHD-Inattentive Type" | Pheno_Comb2$DX_04 == "ADHD-Inattentive Type" | Pheno_Comb2$DX_05 == "ADHD-Inattentive Type" | Pheno_Comb2$DX_06 == "ADHD-Inattentive Type" | Pheno_Comb2$DX_07 == "ADHD-Inattentive Type" | Pheno_Comb2$DX_08 == "ADHD-Inattentive Type" | Pheno_Comb2$DX_09 == "ADHD-Inattentive Type" | Pheno_Comb2$DX_10 == "ADHD-Inattentive Type", "1", "0")
Pheno_Comb2$AdjustmentDisorders_dummy <- ifelse(Pheno_Comb2$DX_01 == "Adjustment Disorders" | Pheno_Comb2$DX_02 == "Adjustment Disorders" | Pheno_Comb2$DX_03 == "Adjustment Disorders" | Pheno_Comb2$DX_04 == "Adjustment Disorders" | Pheno_Comb2$DX_05 == "Adjustment Disorders" | Pheno_Comb2$DX_06 == "Adjustment Disorders" | Pheno_Comb2$DX_07 == "Adjustment Disorders" | Pheno_Comb2$DX_08 == "Adjustment Disorders" | Pheno_Comb2$DX_09 == "Adjustment Disorders" | Pheno_Comb2$DX_10 == "Adjustment Disorders", "1", "0")
Pheno_Comb2$BipolarIIDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Bipolar II Disorder" | Pheno_Comb2$DX_02 == "Bipolar II Disorder" | Pheno_Comb2$DX_03 == "Bipolar II Disorder" | Pheno_Comb2$DX_04 == "Bipolar II Disorder" | Pheno_Comb2$DX_05 == "Bipolar II Disorder" | Pheno_Comb2$DX_06 == "Bipolar II Disorder" | Pheno_Comb2$DX_07 == "Bipolar II Disorder" | Pheno_Comb2$DX_08 == "Bipolar II Disorder" | Pheno_Comb2$DX_09 == "Bipolar II Disorder" | Pheno_Comb2$DX_10 == "Bipolar II Disorder", "1", "0")
Pheno_Comb2$BorderlineIntellectualFunctioning_dummy <- ifelse(Pheno_Comb2$DX_01 == "Borderline Intellectual Functioning" | Pheno_Comb2$DX_02 == "Borderline Intellectual Functioning" | Pheno_Comb2$DX_03 == "Borderline Intellectual Functioning" | Pheno_Comb2$DX_04 == "Borderline Intellectual Functioning" | Pheno_Comb2$DX_05 == "Borderline Intellectual Functioning" | Pheno_Comb2$DX_06 == "Borderline Intellectual Functioning" | Pheno_Comb2$DX_07 == "Borderline Intellectual Functioning" | Pheno_Comb2$DX_08 == "Borderline Intellectual Functioning" | Pheno_Comb2$DX_09 == "Borderline Intellectual Functioning" | Pheno_Comb2$DX_10 == "Borderline Intellectual Functioning", "1", "0")
Pheno_Comb2$ConductDisorderChildhoodonsettype_dummy <- ifelse(Pheno_Comb2$DX_01 == "Conduct Disorder-Childhood-onset type" | Pheno_Comb2$DX_02 == "Conduct Disorder-Childhood-onset type" | Pheno_Comb2$DX_03 == "Conduct Disorder-Childhood-onset type" | Pheno_Comb2$DX_04 == "Conduct Disorder-Childhood-onset type" | Pheno_Comb2$DX_05 == "Conduct Disorder-Childhood-onset type" | Pheno_Comb2$DX_06 == "Conduct Disorder-Childhood-onset type" | Pheno_Comb2$DX_07 == "Conduct Disorder-Childhood-onset type" | Pheno_Comb2$DX_08 == "Conduct Disorder-Childhood-onset type" | Pheno_Comb2$DX_09 == "Conduct Disorder-Childhood-onset type" | Pheno_Comb2$DX_10 == "Conduct Disorder-Childhood-onset type", "1", "0")
Pheno_Comb2$DisruptiveMoodDysregulationDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Disruptive Mood Dysregulation Disorder" | Pheno_Comb2$DX_02 == "Disruptive Mood Dysregulation Disorder" | Pheno_Comb2$DX_03 == "Disruptive Mood Dysregulation Disorder" | Pheno_Comb2$DX_04 == "Disruptive Mood Dysregulation Disorder" | Pheno_Comb2$DX_05 == "Disruptive Mood Dysregulation Disorder" | Pheno_Comb2$DX_06 == "Disruptive Mood Dysregulation Disorder" | Pheno_Comb2$DX_07 == "Disruptive Mood Dysregulation Disorder" | Pheno_Comb2$DX_08 == "Disruptive Mood Dysregulation Disorder" | Pheno_Comb2$DX_09 == "Disruptive Mood Dysregulation Disorder" | Pheno_Comb2$DX_10 == "Disruptive Mood Dysregulation Disorder", "1", "0")
Pheno_Comb2$Encopresis_dummy <- ifelse(Pheno_Comb2$DX_01 == "Encopresis" | Pheno_Comb2$DX_02 == "Encopresis" | Pheno_Comb2$DX_03 == "Encopresis" | Pheno_Comb2$DX_04 == "Encopresis" | Pheno_Comb2$DX_05 == "Encopresis" | Pheno_Comb2$DX_06 == "Encopresis" | Pheno_Comb2$DX_07 == "Encopresis" | Pheno_Comb2$DX_08 == "Encopresis" | Pheno_Comb2$DX_09 == "Encopresis" | Pheno_Comb2$DX_10 == "Encopresis", "1", "0")
Pheno_Comb2$Enuresis_dummy <- ifelse(Pheno_Comb2$DX_01 == "Enuresis" | Pheno_Comb2$DX_02 == "Enuresis" | Pheno_Comb2$DX_03 == "Enuresis" | Pheno_Comb2$DX_04 == "Enuresis" | Pheno_Comb2$DX_05 == "Enuresis" | Pheno_Comb2$DX_06 == "Enuresis" | Pheno_Comb2$DX_07 == "Enuresis" | Pheno_Comb2$DX_08 == "Enuresis" | Pheno_Comb2$DX_09 == "Enuresis" | Pheno_Comb2$DX_10 == "Enuresis", "1", "0")
Pheno_Comb2$ExcoriationDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Excoriation (Skin-Picking) Disorder" | Pheno_Comb2$DX_02 == "Excoriation (Skin-Picking) Disorder" | Pheno_Comb2$DX_03 == "Excoriation (Skin-Picking) Disorder" | Pheno_Comb2$DX_04 == "Excoriation (Skin-Picking) Disorder" | Pheno_Comb2$DX_05 == "Excoriation (Skin-Picking) Disorder" | Pheno_Comb2$DX_06 == "Excoriation (Skin-Picking) Disorder" | Pheno_Comb2$DX_07 == "Excoriation (Skin-Picking) Disorder" | Pheno_Comb2$DX_08 == "Excoriation (Skin-Picking) Disorder" | Pheno_Comb2$DX_09 == "Excoriation (Skin-Picking) Disorder" | Pheno_Comb2$DX_10 == "Excoriation (Skin-Picking) Disorder", "1", "0")
Pheno_Comb2$GeneralizedAnxietyDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Generalized Anxiety Disorder" | Pheno_Comb2$DX_02 == "Generalized Anxiety Disorder" | Pheno_Comb2$DX_03 == "Generalized Anxiety Disorder" | Pheno_Comb2$DX_04 == "Generalized Anxiety Disorder" | Pheno_Comb2$DX_05 == "Generalized Anxiety Disorder" | Pheno_Comb2$DX_06 == "Generalized Anxiety Disorder" | Pheno_Comb2$DX_07 == "Generalized Anxiety Disorder" | Pheno_Comb2$DX_08 == "Generalized Anxiety Disorder" | Pheno_Comb2$DX_09 == "Generalized Anxiety Disorder" | Pheno_Comb2$DX_10 == "Generalized Anxiety Disorder", "1", "0")
Pheno_Comb2$IntellectualDisabilityMild_dummy <- ifelse(Pheno_Comb2$DX_01 == "Intellectual Disability-Mild" | Pheno_Comb2$DX_02 == "Intellectual Disability-Mild" | Pheno_Comb2$DX_03 == "Intellectual Disability-Mild" | Pheno_Comb2$DX_04 == "Intellectual Disability-Mild" | Pheno_Comb2$DX_05 == "Intellectual Disability-Mild" | Pheno_Comb2$DX_06 == "Intellectual Disability-Mild" | Pheno_Comb2$DX_07 == "Intellectual Disability-Mild" | Pheno_Comb2$DX_08 == "Intellectual Disability-Mild" | Pheno_Comb2$DX_09 == "Intellectual Disability-Mild" | Pheno_Comb2$DX_10 == "Intellectual Disability-Mild", "1", "0")
Pheno_Comb2$IntellectualDisabilityModerate_dummy <- ifelse(Pheno_Comb2$DX_01 == "Intellectual Disability-Moderate" | Pheno_Comb2$DX_02 == "Intellectual Disability-Moderate" | Pheno_Comb2$DX_03 == "Intellectual Disability-Moderate" | Pheno_Comb2$DX_04 == "Intellectual Disability-Moderate" | Pheno_Comb2$DX_05 == "Intellectual Disability-Moderate" | Pheno_Comb2$DX_06 == "Intellectual Disability-Moderate" | Pheno_Comb2$DX_07 == "Intellectual Disability-Moderate" | Pheno_Comb2$DX_08 == "Intellectual Disability-Moderate" | Pheno_Comb2$DX_09 == "Intellectual Disability-Moderate" | Pheno_Comb2$DX_10 == "Intellectual Disability-Moderate", "1", "0")
Pheno_Comb2$LanguageDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Language Disorder" | Pheno_Comb2$DX_02 == "Language Disorder" | Pheno_Comb2$DX_03 == "Language Disorder" | Pheno_Comb2$DX_04 == "Language Disorder" | Pheno_Comb2$DX_05 == "Language Disorder" | Pheno_Comb2$DX_06 == "Language Disorder" | Pheno_Comb2$DX_07 == "Language Disorder" | Pheno_Comb2$DX_08 == "Language Disorder" | Pheno_Comb2$DX_09 == "Language Disorder" | Pheno_Comb2$DX_10 == "Language Disorder", "1", "0")
Pheno_Comb2$MajorDepressiveDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Major Depressive Disorder" | Pheno_Comb2$DX_02 == "Major Depressive Disorder" | Pheno_Comb2$DX_03 == "Major Depressive Disorder" | Pheno_Comb2$DX_04 == "Major Depressive Disorder" | Pheno_Comb2$DX_05 == "Major Depressive Disorder" | Pheno_Comb2$DX_06 == "Major Depressive Disorder" | Pheno_Comb2$DX_07 == "Major Depressive Disorder" | Pheno_Comb2$DX_08 == "Major Depressive Disorder" | Pheno_Comb2$DX_09 == "Major Depressive Disorder" | Pheno_Comb2$DX_10 == "Major Depressive Disorder", "1", "0")
Pheno_Comb2$ObsessiveCompulsiveDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Obsessive-Compulsive Disorder" | Pheno_Comb2$DX_02 == "Obsessive-Compulsive Disorder" | Pheno_Comb2$DX_03 == "Obsessive-Compulsive Disorder" | Pheno_Comb2$DX_04 == "Obsessive-Compulsive Disorder" | Pheno_Comb2$DX_05 == "Obsessive-Compulsive Disorder" | Pheno_Comb2$DX_06 == "Obsessive-Compulsive Disorder" | Pheno_Comb2$DX_07 == "Obsessive-Compulsive Disorder" | Pheno_Comb2$DX_08 == "Obsessive-Compulsive Disorder" | Pheno_Comb2$DX_09 == "Obsessive-Compulsive Disorder" | Pheno_Comb2$DX_10 == "Obsessive-Compulsive Disorder", "1", "0")
Pheno_Comb2$OppositionalDefiantDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Oppositional Defiant Disorder" | Pheno_Comb2$DX_02 == "Oppositional Defiant Disorder" | Pheno_Comb2$DX_03 == "Oppositional Defiant Disorder" | Pheno_Comb2$DX_04 == "Oppositional Defiant Disorder" | Pheno_Comb2$DX_05 == "Oppositional Defiant Disorder" | Pheno_Comb2$DX_06 == "Oppositional Defiant Disorder" | Pheno_Comb2$DX_07 == "Oppositional Defiant Disorder" | Pheno_Comb2$DX_08 == "Oppositional Defiant Disorder" | Pheno_Comb2$DX_09 == "Oppositional Defiant Disorder" | Pheno_Comb2$DX_10 == "Oppositional Defiant Disorder", "1", "0")
Pheno_Comb2$OtherSpecifiedAnxietyDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Other Specified Anxiety Disorder" | Pheno_Comb2$DX_02 == "Other Specified Anxiety Disorder" | Pheno_Comb2$DX_03 == "Other Specified Anxiety Disorder" | Pheno_Comb2$DX_04 == "Other Specified Anxiety Disorder" | Pheno_Comb2$DX_05 == "Other Specified Anxiety Disorder" | Pheno_Comb2$DX_06 == "Other Specified Anxiety Disorder" | Pheno_Comb2$DX_07 == "Other Specified Anxiety Disorder" | Pheno_Comb2$DX_08 == "Other Specified Anxiety Disorder" | Pheno_Comb2$DX_09 == "Other Specified Anxiety Disorder" | Pheno_Comb2$DX_10 == "Other Specified Anxiety Disorder", "1", "0")
Pheno_Comb2$OtherSpecifiedAttentionDeficitHyperactivityDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Other Specified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_02 == "Other Specified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_03 == "Other Specified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_04 == "Other Specified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_05 == "Other Specified Attention-Deficit/Hyperactivity Disorder"| Pheno_Comb2$DX_06 == "Other Specified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_07 == "Other Specified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_08 == "Other Specified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_09 == "Other Specified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_10 == "Other Specified Attention-Deficit/Hyperactivity Disorder", "1", "0")
Pheno_Comb2$OtherSpecifiedDepressiveDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Other Specified Depressive Disorder" | Pheno_Comb2$DX_02 == "Other Specified Depressive Disorder" | Pheno_Comb2$DX_03 == "Other Specified Depressive Disorder" | Pheno_Comb2$DX_04 == "Other Specified Depressive Disorder" | Pheno_Comb2$DX_05 == "Other Specified Depressive Disorder"| Pheno_Comb2$DX_06 == "Other Specified Depressive Disorder" | Pheno_Comb2$DX_07 == "Other Specified Depressive Disorder" | Pheno_Comb2$DX_08 == "Other Specified Depressive Disorder" | Pheno_Comb2$DX_09 == "Other Specified Depressive Disorder" | Pheno_Comb2$DX_10 == "Other Specified Depressive Disorder", "1", "0")
Pheno_Comb2$OtherSpecifiedDisruptive_dummy <- ifelse(Pheno_Comb2$DX_01 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder" | Pheno_Comb2$DX_02 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder" | Pheno_Comb2$DX_03 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder" | Pheno_Comb2$DX_04 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder" | Pheno_Comb2$DX_05 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder"| Pheno_Comb2$DX_06 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder" | Pheno_Comb2$DX_07 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder" | Pheno_Comb2$DX_08 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder" | Pheno_Comb2$DX_09 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder" | Pheno_Comb2$DX_10 == "Other Specified Disruptive, Impulse-Control, and Conduct Disorder", "1", "0")
Pheno_Comb2$OtherSpecifiedTicDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Other Specified Tic Disorder" | Pheno_Comb2$DX_02 == "Other Specified Tic Disorder" | Pheno_Comb2$DX_03 == "Other Specified Tic Disorder" | Pheno_Comb2$DX_04 == "Other Specified Tic Disorder" | Pheno_Comb2$DX_05 == "Other Specified Tic Disorder"| Pheno_Comb2$DX_06 == "Other Specified Tic Disorder" | Pheno_Comb2$DX_07 == "Other Specified Tic Disorder" | Pheno_Comb2$DX_08 == "Other Specified Tic Disorder" | Pheno_Comb2$DX_09 == "Other Specified Tic Disorder" | Pheno_Comb2$DX_10 == "Other Specified Tic Disorder", "1", "0")
Pheno_Comb2$OtherSpecifiedTrauma_dummy <- ifelse(Pheno_Comb2$DX_01 == "Other Specified Trauma- and Stressor-Related Disorder" | Pheno_Comb2$DX_02 == "Other Specified Trauma- and Stressor-Related Disorder" | Pheno_Comb2$DX_03 == "Other Specified Trauma- and Stressor-Related Disorder" | Pheno_Comb2$DX_04 == "Other Specified Trauma- and Stressor-Related Disorder" | Pheno_Comb2$DX_05 == "Other Specified Trauma- and Stressor-Related Disorder"| Pheno_Comb2$DX_06 == "Other Specified Trauma- and Stressor-Related Disorder" | Pheno_Comb2$DX_07 == "Other Specified Trauma- and Stressor-Related Disorder" | Pheno_Comb2$DX_08 == "Other Specified Trauma- and Stressor-Related Disorder" | Pheno_Comb2$DX_09 == "Other Specified Trauma- and Stressor-Related Disorder" | Pheno_Comb2$DX_10 == "Other Specified Trauma- and Stressor-Related Disorder", "1", "0")
Pheno_Comb2$PanicDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Panic Disorder" | Pheno_Comb2$DX_02 == "Panic Disorder" | Pheno_Comb2$DX_03 == "Panic Disorder" | Pheno_Comb2$DX_04 == "Panic Disorder" | Pheno_Comb2$DX_05 == "Panic Disorder"| Pheno_Comb2$DX_06 == "Panic Disorder" | Pheno_Comb2$DX_07 == "Panic Disorder" | Pheno_Comb2$DX_08 == "Panic Disorder" | Pheno_Comb2$DX_09 == "Panic Disorder" | Pheno_Comb2$DX_10 == "Panic Disorder", "1", "0")
Pheno_Comb2$ParentChildRelationalProblem_dummy <- ifelse(Pheno_Comb2$DX_01 == "Parent-Child Relational Problem" | Pheno_Comb2$DX_02 == "Parent-Child Relational Problem" | Pheno_Comb2$DX_03 == "Parent-Child Relational Problem" | Pheno_Comb2$DX_04 == "Parent-Child Relational Problem" | Pheno_Comb2$DX_05 == "Parent-Child Relational Problem"| Pheno_Comb2$DX_06 == "Parent-Child Relational Problem" | Pheno_Comb2$DX_07 == "Parent-Child Relational Problem" | Pheno_Comb2$DX_08 == "Parent-Child Relational Problem" | Pheno_Comb2$DX_09 == "Parent-Child Relational Problem" | Pheno_Comb2$DX_10 == "Parent-Child Relational Problem", "1", "0")
Pheno_Comb2$PersistentMotorDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Persistent (Chronic) Motor or Vocal Tic Disorder" | Pheno_Comb2$DX_02 == "Persistent (Chronic) Motor or Vocal Tic Disorder" | Pheno_Comb2$DX_03 == "Persistent (Chronic) Motor or Vocal Tic Disorder" | Pheno_Comb2$DX_04 == "Persistent (Chronic) Motor or Vocal Tic Disorder" | Pheno_Comb2$DX_05 == "Persistent (Chronic) Motor or Vocal Tic Disorder"| Pheno_Comb2$DX_06 == "Persistent (Chronic) Motor or Vocal Tic Disorder" | Pheno_Comb2$DX_07 == "Persistent (Chronic) Motor or Vocal Tic Disorder" | Pheno_Comb2$DX_08 == "Persistent (Chronic) Motor or Vocal Tic Disorder" | Pheno_Comb2$DX_09 == "Persistent (Chronic) Motor or Vocal Tic Disorder" | Pheno_Comb2$DX_10 == "Persistent (Chronic) Motor or Vocal Tic Disorder", "1", "0")
Pheno_Comb2$PersistentDepressiveDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Persistent Depressive Disorder (Dysthymia)" | Pheno_Comb2$DX_02 == "Persistent Depressive Disorder (Dysthymia)" | Pheno_Comb2$DX_03 == "Persistent Depressive Disorder (Dysthymia)" | Pheno_Comb2$DX_04 == "Persistent Depressive Disorder (Dysthymia)" | Pheno_Comb2$DX_05 == "Persistent Depressive Disorder (Dysthymia)" | Pheno_Comb2$DX_06 == "Persistent Depressive Disorder (Dysthymia)" | Pheno_Comb2$DX_07 == "Persistent Depressive Disorder (Dysthymia)" | Pheno_Comb2$DX_08 == "Persistent Depressive Disorder (Dysthymia)" | Pheno_Comb2$DX_09 == "Persistent Depressive Disorder (Dysthymia)" | Pheno_Comb2$DX_10 == "Persistent Depressive Disorder (Dysthymia)", "1", "0")
Pheno_Comb2$PosttraumaticStressDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Posttraumatic Stress Disorder" | Pheno_Comb2$DX_02 == "Posttraumatic Stress Disorder" | Pheno_Comb2$DX_03 == "Posttraumatic Stress Disorder" | Pheno_Comb2$DX_04 == "Posttraumatic Stress Disorder" | Pheno_Comb2$DX_05 == "Posttraumatic Stress Disorder" | Pheno_Comb2$DX_06 == "Posttraumatic Stress Disorder" | Pheno_Comb2$DX_07 == "Posttraumatic Stress Disorder" | Pheno_Comb2$DX_08 == "Posttraumatic Stress Disorder" | Pheno_Comb2$DX_09 == "Posttraumatic Stress Disorder" | Pheno_Comb2$DX_10 == "Posttraumatic Stress Disorder", "1", "0")
Pheno_Comb2$SelectiveMutism_dummy <- ifelse(Pheno_Comb2$DX_01 == "Selective Mutism" | Pheno_Comb2$DX_02 == "Selective Mutism" | Pheno_Comb2$DX_03 == "Selective Mutism" | Pheno_Comb2$DX_04 == "Selective Mutism" | Pheno_Comb2$DX_05 == "Selective Mutism" | Pheno_Comb2$DX_06 == "Selective Mutism" | Pheno_Comb2$DX_07 == "Selective Mutism" | Pheno_Comb2$DX_08 == "Selective Mutism" | Pheno_Comb2$DX_09 == "Selective Mutism" | Pheno_Comb2$DX_10 == "Selective Mutism", "1", "0")
Pheno_Comb2$SeparationAnxiety_dummy <- ifelse(Pheno_Comb2$DX_01 == "Separation Anxiety" | Pheno_Comb2$DX_02 == "Separation Anxiety" | Pheno_Comb2$DX_03 == "Separation Anxiety" | Pheno_Comb2$DX_04 == "Separation Anxiety" | Pheno_Comb2$DX_05 == "Separation Anxiety" | Pheno_Comb2$DX_06 == "Separation Anxiety" | Pheno_Comb2$DX_07 == "Separation Anxiety" | Pheno_Comb2$DX_08 == "Separation Anxiety" | Pheno_Comb2$DX_09 == "Separation Anxiety" | Pheno_Comb2$DX_10 == "Separation Anxiety", "1", "0")
Pheno_Comb2$SocialAnxiety_dummy <- ifelse(Pheno_Comb2$DX_01 == "Social Anxiety (Social Phobia)" | Pheno_Comb2$DX_02 == "Social Anxiety (Social Phobia)" | Pheno_Comb2$DX_03 == "Social Anxiety (Social Phobia)" | Pheno_Comb2$DX_04 == "Social Anxiety (Social Phobia)" | Pheno_Comb2$DX_05 == "Social Anxiety (Social Phobia)" | Pheno_Comb2$DX_06 == "Social Anxiety (Social Phobia)" | Pheno_Comb2$DX_07 == "Social Anxiety (Social Phobia)" | Pheno_Comb2$DX_08 == "Social Anxiety (Social Phobia)" | Pheno_Comb2$DX_09 == "Social Anxiety (Social Phobia)" | Pheno_Comb2$DX_10 == "Social Anxiety (Social Phobia)", "1", "0")
Pheno_Comb2$SLDIMath_dummy <- ifelse(Pheno_Comb2$DX_01 == "Specific Learning Disorder with Impairment in Mathematics" | Pheno_Comb2$DX_02 == "Specific Learning Disorder with Impairment in Mathematics" | Pheno_Comb2$DX_03 == "Specific Learning Disorder with Impairment in Mathematics" | Pheno_Comb2$DX_04 == "Specific Learning Disorder with Impairment in Mathematics" | Pheno_Comb2$DX_05 == "Specific Learning Disorder with Impairment in Mathematics" | Pheno_Comb2$DX_06 == "Specific Learning Disorder with Impairment in Mathematics" | Pheno_Comb2$DX_07 == "Specific Learning Disorder with Impairment in Mathematics" | Pheno_Comb2$DX_08 == "Specific Learning Disorder with Impairment in Mathematics" | Pheno_Comb2$DX_09 == "Specific Learning Disorder with Impairment in Mathematics" | Pheno_Comb2$DX_10 == "Specific Learning Disorder with Impairment in Mathematics", "1", "0")
Pheno_Comb2$SLDReading_dummy <- ifelse(Pheno_Comb2$DX_01 == "Specific Learning Disorder with Impairment in Reading" | Pheno_Comb2$DX_02 == "Specific Learning Disorder with Impairment in Reading" | Pheno_Comb2$DX_03 == "Specific Learning Disorder with Impairment in Reading" | Pheno_Comb2$DX_04 == "Specific Learning Disorder with Impairment in Reading" | Pheno_Comb2$DX_05 == "Specific Learning Disorder with Impairment in Reading" | Pheno_Comb2$DX_06 == "Specific Learning Disorder with Impairment in Reading" | Pheno_Comb2$DX_07 == "Specific Learning Disorder with Impairment in Reading" | Pheno_Comb2$DX_08 == "Specific Learning Disorder with Impairment in Reading" | Pheno_Comb2$DX_09 == "Specific Learning Disorder with Impairment in Reading" | Pheno_Comb2$DX_10 == "Specific Learning Disorder with Impairment in Reading", "1", "0")
Pheno_Comb2$SLDWrittenExpression_dummy <- ifelse(Pheno_Comb2$DX_01 == "Specific Learning Disorder with Impairment in Written Expression" | Pheno_Comb2$DX_02 == "Specific Learning Disorder with Impairment in Written Expression" | Pheno_Comb2$DX_03 == "Specific Learning Disorder with Impairment in Written Expression" | Pheno_Comb2$DX_04 == "Specific Learning Disorder with Impairment in Written Expression" | Pheno_Comb2$DX_05 == "Specific Learning Disorder with Impairment in Written Expression" | Pheno_Comb2$DX_06 == "Specific Learning Disorder with Impairment in Written Expression" | Pheno_Comb2$DX_07 == "Specific Learning Disorder with Impairment in Written Expression" | Pheno_Comb2$DX_08 == "Specific Learning Disorder with Impairment in Written Expression" | Pheno_Comb2$DX_09 == "Specific Learning Disorder with Impairment in Written Expression" | Pheno_Comb2$DX_10 == "Specific Learning Disorder with Impairment in Written Expression", "1", "0")
Pheno_Comb2$SpecificPhobia_dummy <- ifelse(Pheno_Comb2$DX_01 == "Specific Phobia" | Pheno_Comb2$DX_02 == "Specific Phobia" | Pheno_Comb2$DX_03 == "Specific Phobia" | Pheno_Comb2$DX_04 == "Specific Phobia" | Pheno_Comb2$DX_05 == "Specific Phobia" | Pheno_Comb2$DX_06 == "Specific Phobia" | Pheno_Comb2$DX_07 == "Specific Phobia" | Pheno_Comb2$DX_08 == "Specific Phobia" | Pheno_Comb2$DX_09 == "Specific Phobia" | Pheno_Comb2$DX_10 == "Specific Phobia", "1", "0")
Pheno_Comb2$SpeechSoundDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Speech Sound Disorder" | Pheno_Comb2$DX_02 == "Speech Sound Disorder" | Pheno_Comb2$DX_03 == "Speech Sound Disorder" | Pheno_Comb2$DX_04 == "Speech Sound Disorder" | Pheno_Comb2$DX_05 == "Speech Sound Disorder" | Pheno_Comb2$DX_06 == "Speech Sound Disorder" | Pheno_Comb2$DX_07 == "Speech Sound Disorder" | Pheno_Comb2$DX_08 == "Speech Sound Disorder" | Pheno_Comb2$DX_09 == "Speech Sound Disorder" | Pheno_Comb2$DX_10 == "Speech Sound Disorder", "1", "0")
Pheno_Comb2$TourettesDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Tourettes Disorder" | Pheno_Comb2$DX_02 == "Tourettes Disorder" | Pheno_Comb2$DX_03 == "Tourettes Disorder" | Pheno_Comb2$DX_04 == "Tourettes Disorder" | Pheno_Comb2$DX_05 == "Tourettes Disorder" | Pheno_Comb2$DX_06 == "Tourettes Disorder" | Pheno_Comb2$DX_07 == "Tourettes Disorder" | Pheno_Comb2$DX_08 == "Tourettes Disorder" | Pheno_Comb2$DX_09 == "Tourettes Disorder" | Pheno_Comb2$DX_10 == "Tourettes Disorder", "1", "0")
Pheno_Comb2$UnspecifiedAnxietyDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Unspecified Anxiety Disorder" | Pheno_Comb2$DX_02 == "Unspecified Anxiety Disorder" | Pheno_Comb2$DX_03 == "Unspecified Anxiety Disorder" | Pheno_Comb2$DX_04 == "Unspecified Anxiety Disorder" | Pheno_Comb2$DX_05 == "Unspecified Anxiety Disorder" | Pheno_Comb2$DX_06 == "Unspecified Anxiety Disorder" | Pheno_Comb2$DX_07 == "Unspecified Anxiety Disorder" | Pheno_Comb2$DX_08 == "Unspecified Anxiety Disorder" | Pheno_Comb2$DX_09 == "Unspecified Anxiety Disorder" | Pheno_Comb2$DX_10 == "Unspecified Anxiety Disorder", "1", "0")
Pheno_Comb2$UnspecifiedADHD_dummy <- ifelse(Pheno_Comb2$DX_01 == "Unspecified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_02 == "Unspecified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_03 == "Unspecified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_04 == "Unspecified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_05 == "Unspecified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_06 == "Unspecified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_07 == "Unspecified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_08 == "Unspecified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_09 == "Unspecified Attention-Deficit/Hyperactivity Disorder" | Pheno_Comb2$DX_10 == "Unspecified Attention-Deficit/Hyperactivity Disorder", "1", "0")
Pheno_Comb2$UnspecifiedNeurodevelopmentalDisorder_dummy <- ifelse(Pheno_Comb2$DX_01 == "Unspecified Neurodevelopmental Disorder" | Pheno_Comb2$DX_02 == "Unspecified Neurodevelopmental Disorder" | Pheno_Comb2$DX_03 == "Unspecified Neurodevelopmental Disorder" | Pheno_Comb2$DX_04 == "Unspecified Neurodevelopmental Disorder" | Pheno_Comb2$DX_05 == "Unspecified Neurodevelopmental Disorder" | Pheno_Comb2$DX_06 == "Unspecified Neurodevelopmental Disorder" | Pheno_Comb2$DX_07 == "Unspecified Neurodevelopmental Disorder" | Pheno_Comb2$DX_08 == "Unspecified Neurodevelopmental Disorder" | Pheno_Comb2$DX_09 == "Unspecified Neurodevelopmental Disorder" | Pheno_Comb2$DX_10 == "Unspecified Neurodevelopmental Disorder", "1", "0")

#Get subject counts
Pheno_Comb2$ASD_dummy <- as.numeric(Pheno_Comb2$ASD_dummy)  
sum(Pheno_Comb2$ASD_dummy, na.rm=TRUE) #85

Pheno_Comb2$NoDX_dummy <- as.numeric(Pheno_Comb2$NoDX_dummy)  
sum(Pheno_Comb2$NoDX_dummy, na.rm=TRUE) #60

Pheno_Comb2$ADHDCombinedType_dummy <- as.numeric(Pheno_Comb2$ADHDCombinedType_dummy)
sum(Pheno_Comb2$ADHDCombinedType_dummy, na.rm=TRUE) #152

Pheno_Comb2$ADHDHyperactiveImpulsiveType_dummy <- as.numeric(Pheno_Comb2$ADHDHyperactiveImpulsiveType_dummy)
sum(Pheno_Comb2$ADHDHyperactiveImpulsiveType_dummy, na.rm=TRUE) #14

Pheno_Comb2$ADHDInattentiveType_dummy <- as.numeric(Pheno_Comb2$ADHDInattentiveType_dummy)
sum(Pheno_Comb2$ADHDInattentiveType_dummy, na.rm=TRUE) #130

Pheno_Comb2$AdjustmentDisorders_dummy <- as.numeric(Pheno_Comb2$AdjustmentDisorders_dummy)
sum(Pheno_Comb2$AdjustmentDisorders_dummy, na.rm=TRUE) #14

Pheno_Comb2$BipolarIIDisorder_dummy <- as.numeric(Pheno_Comb2$BipolarIIDisorder_dummy)
sum(Pheno_Comb2$BipolarIIDisorder_dummy, na.rm=TRUE) #1

Pheno_Comb2$BorderlineIntellectualFunctioning_dummy <- as.numeric(Pheno_Comb2$BorderlineIntellectualFunctioning_dummy)
sum(Pheno_Comb2$BorderlineIntellectualFunctioning_dummy, na.rm=TRUE) #2

Pheno_Comb2$ConductDisorderChildhoodonsettype_dummy <- as.numeric(Pheno_Comb2$ConductDisorderChildhoodonsettype_dummy)
sum(Pheno_Comb2$ConductDisorderChildhoodonsettype_dummy, na.rm=TRUE) #1

Pheno_Comb2$DisruptiveMoodDysregulationDisorder_dummy <- as.numeric(Pheno_Comb2$DisruptiveMoodDysregulationDisorder_dummy)
sum(Pheno_Comb2$DisruptiveMoodDysregulationDisorder_dummy, na.rm=TRUE) #7

Pheno_Comb2$Encopresis_dummy <- as.numeric(Pheno_Comb2$Encopresis_dummy)
sum(Pheno_Comb2$Encopresis_dummy, na.rm=TRUE) #8

Pheno_Comb2$Enuresis_dummy <- as.numeric(Pheno_Comb2$Enuresis_dummy)
sum(Pheno_Comb2$Enuresis_dummy, na.rm=TRUE) #23

Pheno_Comb2$ExcoriationDisorder_dummy <- as.numeric(Pheno_Comb2$ExcoriationDisorder_dummy)
sum(Pheno_Comb2$ExcoriationDisorder_dummy, na.rm=TRUE) #10

Pheno_Comb2$GeneralizedAnxietyDisorder_dummy <- as.numeric(Pheno_Comb2$GeneralizedAnxietyDisorder_dummy)
sum(Pheno_Comb2$GeneralizedAnxietyDisorder_dummy, na.rm=TRUE) #63

Pheno_Comb2$IntellectualDisabilityMild_dummy <- as.numeric(Pheno_Comb2$IntellectualDisabilityMild_dummy)
sum(Pheno_Comb2$IntellectualDisabilityMild_dummy, na.rm=TRUE) #14

Pheno_Comb2$IntellectualDisabilityModerate_dummy <- as.numeric(Pheno_Comb2$IntellectualDisabilityModerate_dummy)
sum(Pheno_Comb2$IntellectualDisabilityModerate_dummy, na.rm=TRUE) #1

Pheno_Comb2$LanguageDisorder_dummy <- as.numeric(Pheno_Comb2$LanguageDisorder_dummy)
sum(Pheno_Comb2$LanguageDisorder_dummy, na.rm=TRUE) #39

Pheno_Comb2$MajorDepressiveDisorder_dummy <- as.numeric(Pheno_Comb2$MajorDepressiveDisorder_dummy)
sum(Pheno_Comb2$MajorDepressiveDisorder_dummy, na.rm=TRUE) #26

Pheno_Comb2$ObsessiveCompulsiveDisorder_dummy <- as.numeric(Pheno_Comb2$ObsessiveCompulsiveDisorder_dummy)
sum(Pheno_Comb2$ObsessiveCompulsiveDisorder_dummy, na.rm=TRUE) #18

Pheno_Comb2$OppositionalDefiantDisorder_dummy <- as.numeric(Pheno_Comb2$OppositionalDefiantDisorder_dummy)
sum(Pheno_Comb2$OppositionalDefiantDisorder_dummy, na.rm=TRUE) #33

Pheno_Comb2$OtherSpecifiedAnxietyDisorder_dummy <- as.numeric(Pheno_Comb2$OtherSpecifiedAnxietyDisorder_dummy)
sum(Pheno_Comb2$OtherSpecifiedAnxietyDisorder_dummy, na.rm=TRUE) #28

Pheno_Comb2$OtherSpecifiedAttentionDeficitHyperactivityDisorder_dummy <- as.numeric(Pheno_Comb2$OtherSpecifiedAttentionDeficitHyperactivityDisorder_dummy)
sum(Pheno_Comb2$OtherSpecifiedAttentionDeficitHyperactivityDisorder_dummy, na.rm=TRUE) #13

Pheno_Comb2$OtherSpecifiedDepressiveDisorder_dummy <- as.numeric(Pheno_Comb2$OtherSpecifiedDepressiveDisorder_dummy)
sum(Pheno_Comb2$OtherSpecifiedDepressiveDisorder_dummy, na.rm=TRUE) #7

Pheno_Comb2$OtherSpecifiedDisruptive_dummy <- as.numeric(Pheno_Comb2$OtherSpecifiedDisruptive_dummy)
sum(Pheno_Comb2$OtherSpecifiedDisruptive_dummy, na.rm=TRUE) #1

Pheno_Comb2$OtherSpecifiedTicDisorder_dummy <- as.numeric(Pheno_Comb2$OtherSpecifiedTicDisorder_dummy)
sum(Pheno_Comb2$OtherSpecifiedTicDisorder_dummy, na.rm=TRUE) #2

Pheno_Comb2$OtherSpecifiedTrauma_dummy <- as.numeric(Pheno_Comb2$OtherSpecifiedTrauma_dummy)
sum(Pheno_Comb2$OtherSpecifiedTrauma_dummy, na.rm=TRUE) #2

Pheno_Comb2$PanicDisorder_dummy <- as.numeric(Pheno_Comb2$PanicDisorder_dummy)
sum(Pheno_Comb2$PanicDisorder_dummy, na.rm=TRUE) #1

Pheno_Comb2$ParentChildRelationalProblem_dummy <- as.numeric(Pheno_Comb2$ParentChildRelationalProblem_dummy)
sum(Pheno_Comb2$ParentChildRelationalProblem_dummy, na.rm=TRUE) #0

Pheno_Comb2$PersistentMotorDisorder_dummy <- as.numeric(Pheno_Comb2$PersistentMotorDisorder_dummy)
sum(Pheno_Comb2$PersistentMotorDisorder_dummy, na.rm=TRUE) #14

Pheno_Comb2$PersistentDepressiveDisorder_dummy <- as.numeric(Pheno_Comb2$PersistentDepressiveDisorder_dummy)
sum(Pheno_Comb2$PersistentDepressiveDisorder_dummy, na.rm=TRUE) #5

Pheno_Comb2$PosttraumaticStressDisorder_dummy <- as.numeric(Pheno_Comb2$PosttraumaticStressDisorder_dummy)
sum(Pheno_Comb2$Pheno_Comb2$PosttraumaticStressDisorder_dummy, na.rm=TRUE) #0

Pheno_Comb2$SelectiveMutism_dummy <- as.numeric(Pheno_Comb2$SelectiveMutism_dummy)
sum(Pheno_Comb2$SelectiveMutism_dummy, na.rm=TRUE) #2

Pheno_Comb2$SeparationAnxiety_dummy <- as.numeric(Pheno_Comb2$SeparationAnxiety_dummy)
sum(Pheno_Comb2$SeparationAnxiety_dummy, na.rm=TRUE) #21

Pheno_Comb2$SocialAnxiety_dummy <- as.numeric(Pheno_Comb2$SocialAnxiety_dummy)
sum(Pheno_Comb2$SocialAnxiety_dummy, na.rm=TRUE) #51

Pheno_Comb2$SLDIMath_dummy <- as.numeric(Pheno_Comb2$SLDIMath_dummy)
sum(Pheno_Comb2$SLDIMath_dummy, na.rm=TRUE) #30

Pheno_Comb2$SLDReading_dummy <- as.numeric(Pheno_Comb2$SLDReading_dummy)
sum(Pheno_Comb2$SLDReading_dummy, na.rm=TRUE) #91

Pheno_Comb2$SLDWrittenExpression_dummy <- as.numeric(Pheno_Comb2$SLDWrittenExpression_dummy)
sum(Pheno_Comb2$SLDWrittenExpression_dummy, na.rm=TRUE) #20

Pheno_Comb2$SpecificPhobia_dummy <- as.numeric(Pheno_Comb2$SpecificPhobia_dummy)
sum(Pheno_Comb2$SpecificPhobia_dummy, na.rm=TRUE) #47

Pheno_Comb2$SpeechSoundDisorder_dummy <- as.numeric(Pheno_Comb2$SpeechSoundDisorder_dummy)
sum(Pheno_Comb2$SpeechSoundDisorder_dummy, na.rm=TRUE) #10

Pheno_Comb2$TourettesDisorder_dummy <- as.numeric(Pheno_Comb2$TourettesDisorder_dummy)
sum(Pheno_Comb2$TourettesDisorder_dummy, na.rm=TRUE) #13

Pheno_Comb2$UnspecifiedAnxietyDisorder_dummy <- as.numeric(Pheno_Comb2$UnspecifiedAnxietyDisorder_dummy)
sum(Pheno_Comb2$UnspecifiedAnxietyDisorder_dummy, na.rm=TRUE) #4

Pheno_Comb2$UnspecifiedADHD_dummy <- as.numeric(Pheno_Comb2$UnspecifiedADHD_dummy)
sum(Pheno_Comb2$UnspecifiedADHD_dummy, na.rm=TRUE) #1

Pheno_Comb2$UnspecifiedNeurodevelopmentalDisorder_dummy <- as.numeric(Pheno_Comb2$UnspecifiedNeurodevelopmentalDisorder_dummy)
sum(Pheno_Comb2$UnspecifiedNeurodevelopmentalDisorder_dummy, na.rm=TRUE) #1


#Dx that survive the 50 cutoff:
#ASD, NoDx, ADHD Combined, ADHD Inattentive, Generalized Anxiety Disorder, Social Anxiety, SLD Reading


#----------------------------------------------CREATE GROUPS----------------------------------------------------------
#Drop Dx groups with less than 50 (or preserve those that have 50+)
Pheno_Comb2$KeepDx <- ifelse(Pheno_Comb2$ASD_dummy == "1" | Pheno_Comb2$NoDX_dummy == "1" | Pheno_Comb2$ADHDCombinedType_dummy=="1" | Pheno_Comb2$ADHDInattentiveType_dummy =="1" | Pheno_Comb2$GeneralizedAnxietyDisorder_dummy =="1" | Pheno_Comb2$SocialAnxiety_dummy=="1" | Pheno_Comb2$SLDReading_dummy=="1", 1, 0)
Pheno_Comb2 <- subset(Pheno_Comb2, KeepDx!="0")


#Create var for any ASD diagnosis
Pheno_Comb2$ASDany <- ifelse(Pheno_Comb2$DX_01 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_02 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_03 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_04 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_05 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_06 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_07 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_08 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_09 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_10 == "Autism Spectrum Disorder", "ASD", NA)
#Create var for no dx
Pheno_Comb2$NoDX <- ifelse(Pheno_Comb2$DX_01 == "No Diagnosis Given" | Pheno_Comb2$DX_02 == "No Diagnosis Given" | Pheno_Comb2$DX_03 == "No Diagnosis Given" | Pheno_Comb2$DX_04 == "No Diagnosis Given" | Pheno_Comb2$DX_05 == "No Diagnosis Given" | Pheno_Comb2$DX_06 == "No Diagnosis Given" | Pheno_Comb2$DX_07 == "No Diagnosis Given" | Pheno_Comb2$DX_08 == "No Diagnosis Given" | Pheno_Comb2$DX_09 == "No Diagnosis Given" | Pheno_Comb2$DX_10 == "No Diagnosis Given", "NoDX", NA)
#Create var for other Dx
Pheno_Comb2$NotASDdx <- ifelse(is.na(Pheno_Comb2$ASDany) & is.na(Pheno_Comb2$NoDX), "OtherDx", " ")

#Create var for any ASD diagnosis
Pheno_Comb2$ASDany <- ifelse(Pheno_Comb2$DX_01 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_02 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_03 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_04 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_05 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_06 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_07 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_08 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_09 == "Autism Spectrum Disorder" | Pheno_Comb2$DX_10 == "Autism Spectrum Disorder", "ASD", " ")
#Create var for no dx
Pheno_Comb2$NoDX <- ifelse(Pheno_Comb2$DX_01 == "No Diagnosis Given" | Pheno_Comb2$DX_02 == "No Diagnosis Given" | Pheno_Comb2$DX_03 == "No Diagnosis Given" | Pheno_Comb2$DX_04 == "No Diagnosis Given" | Pheno_Comb2$DX_05 == "No Diagnosis Given" | Pheno_Comb2$DX_06 == "No Diagnosis Given" | Pheno_Comb2$DX_07 == "No Diagnosis Given" | Pheno_Comb2$DX_08 == "No Diagnosis Given" | Pheno_Comb2$DX_09 == "No Diagnosis Given" | Pheno_Comb2$DX_10 == "No Diagnosis Given", "NoDX", " ")

#Create overall group var
Pheno_Comb2$group <- (paste0(Pheno_Comb2$ASDany, Pheno_Comb2$NoDX, Pheno_Comb2$NotASDdx))
#Replace NANAOtherDx with OtherDx
Pheno_Comb2$group <- str_replace(Pheno_Comb2$group, "NANA", "")
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
Pheno_Comb2$group <- trim(Pheno_Comb2$group)


#Dx Count Table
table(Pheno_Comb2$group, Pheno_Comb2$T1)

#--------------------------------VARIABLES-------------------------------------

#Create variables

#Group
Pheno_Comb2$group <- as.factor(Pheno_Comb2$group)

#Sex
Pheno_Comb2$Sex <- as.factor(Pheno_Comb2$Sex)

#Site
Pheno_Comb2$Site <- as.factor(Pheno_Comb2$Site)

#Create CSF cm3 var
Pheno_Comb2$EACSF_cm <- Pheno_Comb2$rawEACSF/1000

#Create TBV cm3 var
Pheno_Comb2$BrainSegVolcm3 <- Pheno_Comb2$BrainSegVol/1000
Pheno_Comb2$BrainSegVolcm3_squared <- Pheno_Comb2$BrainSegVolcm3 * Pheno_Comb2$BrainSegVolcm3

#Create eTIV cm3 var
Pheno_Comb2$eTIV_cm <- Pheno_Comb2$EstimatedTotalIntraCranialVol/1000
Pheno_Comb2$eTIV_cm_squared <- Pheno_Comb2$eTIV_cm * Pheno_Comb2$eTIV_cm

#Create age squared
Pheno_Comb2$Age <- as.numeric(as.character(Pheno_Comb2$Age)) #converts age to "numeric" designation
Pheno_Comb2$Age_squared <- Pheno_Comb2$Age*Pheno_Comb2$Age

#Make age and age_squared mean-centered
meanage = mean(Pheno_Comb2$Age)
Pheno_Comb2$mc_age <- Pheno_Comb2$Age - meanage
Pheno_Comb2$mc_age_squared <- Pheno_Comb2$mc_age * Pheno_Comb2$mc_age
Pheno_Comb2$mc_age_squared <- as.numeric(Pheno_Comb2$mc_age_squared)

#Make etiv and etiv-squared mean-centered
etiv_mean = mean(Pheno_Comb2$eTIV_cm)
Pheno_Comb2$eTIV_mc <- Pheno_Comb2$eTIV_cm - etiv_mean
Pheno_Comb2$eTIV_squared_mc <- Pheno_Comb2$eTIV_mc * Pheno_Comb2$eTIV_mc
Pheno_Comb2$totalmean_CT <- (Pheno_Comb2$rh_MeanThickness_thickness + Pheno_Comb2$lh_MeanThickness_thickness)/2



#---------------------------------OUTLIERS-------------------------------------

#Outliers -- fence
#1. Determine percentiles of EACSF_cm *for every year* using favstats (we are taking variability due to age into account)

#Create age bins
Pheno_Comb2$Age <- as.numeric(Pheno_Comb2$Age)
Pheno_Comb2$Age_bins <- 9
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "5" & Pheno_Comb2$Age < "6", "5", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "6" & Pheno_Comb2$Age < "7", "6", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "7" & Pheno_Comb2$Age < "8", "7", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "8" & Pheno_Comb2$Age < "9", "8", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "9" & Pheno_Comb2$Age < "10", "9", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "10" & Pheno_Comb2$Age < "11", "10", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "11" & Pheno_Comb2$Age < "12", "11", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "12" & Pheno_Comb2$Age < "13", "12", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "13" & Pheno_Comb2$Age < "14", "13", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "14" & Pheno_Comb2$Age < "15", "14", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "15" & Pheno_Comb2$Age < "16", "15", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "16" & Pheno_Comb2$Age < "17", "16", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "17" & Pheno_Comb2$Age < "18", "17", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "18" & Pheno_Comb2$Age < "19", "18", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "19" & Pheno_Comb2$Age < "20", "19", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "20" & Pheno_Comb2$Age < "21", "20", Pheno_Comb2$Age_bins)
Pheno_Comb2$Age_bins <- ifelse(Pheno_Comb2$Age >= "21" & Pheno_Comb2$Age < "22", "21", Pheno_Comb2$Age_bins)

#Table to verify that Age_bins worked
table(Pheno_Comb2$Age_bins)

#Favstats
favstats(EACSF_cm ~ Age_bins, data=Pheno_Comb2)

#1.5 Subset dataset to each age_bin
age5 <- subset(Pheno_Comb2, Age_bins==5)
age6 <- subset(Pheno_Comb2, Age_bins==6)
age7 <- subset(Pheno_Comb2, Age_bins==7)
age8 <- subset(Pheno_Comb2, Age_bins==8)
age9 <- subset(Pheno_Comb2, Age_bins==9)
age10 <- subset(Pheno_Comb2, Age_bins==10)
age11 <- subset(Pheno_Comb2, Age_bins==11)
age12 <- subset(Pheno_Comb2, Age_bins==12)
age13 <- subset(Pheno_Comb2, Age_bins==13)
age14 <- subset(Pheno_Comb2, Age_bins==14)
age15 <- subset(Pheno_Comb2, Age_bins==15)
age16 <- subset(Pheno_Comb2, Age_bins==16)
age17 <- subset(Pheno_Comb2, Age_bins==17)
age18 <- subset(Pheno_Comb2, Age_bins==18)
age19 <- subset(Pheno_Comb2, Age_bins==19)
age20 <- subset(Pheno_Comb2, Age_bins==20)
age21 <- subset(Pheno_Comb2, Age_bins==21)


#2. Fence values below Q1- 1.5IQR. Replace the values with Q1 - (1.5*(Q3 - Q1))

#Create first quartile
q1_5 <- unname(quantile(age5$EACSF_cm, 0.25))
q1_6 <- unname(quantile(age6$EACSF_cm, 0.25))
q1_7 <- unname(quantile(age7$EACSF_cm, 0.25))
q1_8 <- unname(quantile(age8$EACSF_cm, 0.25))
q1_9 <- unname(quantile(age9$EACSF_cm, 0.25))
q1_10 <- unname(quantile(age10$EACSF_cm, 0.25))
q1_11 <- unname(quantile(age11$EACSF_cm, 0.25))
q1_12 <- unname(quantile(age12$EACSF_cm, 0.25))
q1_13 <- unname(quantile(age13$EACSF_cm, 0.25))
q1_14 <- unname(quantile(age14$EACSF_cm, 0.25))
q1_15 <- unname(quantile(age15$EACSF_cm, 0.25))
q1_16 <- unname(quantile(age16$EACSF_cm, 0.25))
q1_17 <- unname(quantile(age17$EACSF_cm, 0.25))
q1_18 <- unname(quantile(age18$EACSF_cm, 0.25))
q1_19 <- unname(quantile(age19$EACSF_cm, 0.25))
q1_20 <- unname(quantile(age20$EACSF_cm, 0.25))
q1_21 <- unname(quantile(age21$EACSF_cm, 0.25))


#Create 3rd quartile
q3_5 <- unname(quantile(age5$EACSF_cm, 0.75))
q3_6 <- unname(quantile(age6$EACSF_cm, 0.75))
q3_7 <- unname(quantile(age7$EACSF_cm, 0.75))
q3_8 <- unname(quantile(age8$EACSF_cm, 0.75))
q3_9 <- unname(quantile(age9$EACSF_cm, 0.75))
q3_10 <- unname(quantile(age10$EACSF_cm, 0.75))
q3_11 <- unname(quantile(age11$EACSF_cm, 0.75))
q3_12 <- unname(quantile(age12$EACSF_cm, 0.75))
q3_13 <- unname(quantile(age13$EACSF_cm, 0.75))
q3_14 <- unname(quantile(age14$EACSF_cm, 0.75))
q3_15 <- unname(quantile(age15$EACSF_cm, 0.75))
q3_16 <- unname(quantile(age16$EACSF_cm, 0.75))
q3_17 <- unname(quantile(age17$EACSF_cm, 0.75))
q3_18 <- unname(quantile(age18$EACSF_cm, 0.75))
q3_19 <- unname(quantile(age19$EACSF_cm, 0.75))
q3_20 <- unname(quantile(age20$EACSF_cm, 0.75))
q3_21 <- unname(quantile(age21$EACSF_cm, 0.75))


#Create lower fence
lowerfence5=q1_5 - (1.5*(q3_5 - q1_5))
lowerfence6=q1_6 - (1.5*(q3_6 - q1_6))
lowerfence7=q1_7 - (1.5*(q3_7 - q1_7))
lowerfence8=q1_8 - (1.5*(q3_8 - q1_8))
lowerfence9=q1_9 - (1.5*(q3_9 - q1_9))
lowerfence10=q1_10 - (1.5*(q3_10 - q1_10))
lowerfence11=q1_11 - (1.5*(q3_11 - q1_11))
lowerfence12=q1_12 - (1.5*(q3_12 - q1_12))
lowerfence13=q1_13 - (1.5*(q3_13 - q1_13))
lowerfence14=q1_14 - (1.5*(q3_14 - q1_14))
lowerfence15=q1_15 - (1.5*(q3_15 - q1_15))
lowerfence16=q1_16 - (1.5*(q3_16 - q1_16))
lowerfence17=q1_17 - (1.5*(q3_17 - q1_17))
lowerfence18=q1_18 - (1.5*(q3_18 - q1_18))
lowerfence19=q1_19 - (1.5*(q3_19 - q1_19))
lowerfence20=q1_20 - (1.5*(q3_20 - q1_20))
lowerfence21=q1_21 - (1.5*(q3_21 - q1_21))


#Replace bottom outliers with lower fence by age bin   
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "5" & Pheno_Comb2$EACSF_cm <= lowerfence5, lowerfence5, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "6" & Pheno_Comb2$EACSF_cm <= lowerfence6, lowerfence6, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "7" & Pheno_Comb2$EACSF_cm <= lowerfence7, lowerfence7, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "8" & Pheno_Comb2$EACSF_cm <= lowerfence8, lowerfence8, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "9" & Pheno_Comb2$EACSF_cm <= lowerfence9, lowerfence9, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "10" & Pheno_Comb2$EACSF_cm <= lowerfence10, lowerfence10, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "11" & Pheno_Comb2$EACSF_cm <= lowerfence11, lowerfence11, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "12" & Pheno_Comb2$EACSF_cm <= lowerfence12, lowerfence12, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "13" & Pheno_Comb2$EACSF_cm <= lowerfence13, lowerfence13, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "14" & Pheno_Comb2$EACSF_cm <= lowerfence14, lowerfence14, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "15" & Pheno_Comb2$EACSF_cm <= lowerfence15, lowerfence15, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "16" & Pheno_Comb2$EACSF_cm <= lowerfence16, lowerfence16, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "17" & Pheno_Comb2$EACSF_cm <= lowerfence17, lowerfence17, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "18" & Pheno_Comb2$EACSF_cm <= lowerfence18, lowerfence18, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "19" & Pheno_Comb2$EACSF_cm <= lowerfence19, lowerfence19, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "20" & Pheno_Comb2$EACSF_cm <= lowerfence20, lowerfence20, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins == "21" & Pheno_Comb2$EACSF_cm <= lowerfence21, lowerfence21, Pheno_Comb2$EACSF_cm)

#Check for number of changes made
table(Pheno_Comb2$EACSF_cm==lowerfence5)
table(Pheno_Comb2$EACSF_cm==lowerfence6)
table(Pheno_Comb2$EACSF_cm==lowerfence7)
table(Pheno_Comb2$EACSF_cm==lowerfence8)
table(Pheno_Comb2$EACSF_cm==lowerfence9)
table(Pheno_Comb2$EACSF_cm==lowerfence10)
table(Pheno_Comb2$EACSF_cm==lowerfence11)
table(Pheno_Comb2$EACSF_cm==lowerfence12)
table(Pheno_Comb2$EACSF_cm==lowerfence13)
table(Pheno_Comb2$EACSF_cm==lowerfence14)
table(Pheno_Comb2$EACSF_cm==lowerfence15)
table(Pheno_Comb2$EACSF_cm==lowerfence16)
table(Pheno_Comb2$EACSF_cm==lowerfence17)
table(Pheno_Comb2$EACSF_cm==lowerfence18)
table(Pheno_Comb2$EACSF_cm==lowerfence19)
table(Pheno_Comb2$EACSF_cm==lowerfence20)
table(Pheno_Comb2$EACSF_cm==lowerfence21)

#3. Fence values above Q3 + 1.5IQR. Replace the values with Q3 + (1.5*(Q3 - Q1))

#Create upper fence
upperfence5=q3_5 + (1.5*(q3_5 - q1_5))
upperfence6=q3_6 + (1.5*(q3_6 - q1_6))
upperfence7=q3_7 + (1.5*(q3_7 - q1_7))
upperfence8=q3_8 + (1.5*(q3_8 - q1_8))
upperfence9=q3_9 + (1.5*(q3_9 - q1_9))
upperfence10=q3_10 + (1.5*(q3_10 - q1_10))
upperfence11=q3_11 + (1.5*(q3_11 - q1_11))
upperfence12=q3_12 + (1.5*(q3_12 - q1_12))
upperfence13=q3_13 + (1.5*(q3_13 - q1_13))
upperfence14=q3_14 + (1.5*(q3_14 - q1_14))
upperfence15=q3_15 + (1.5*(q3_15 - q1_15))
upperfence16=q3_16 + (1.5*(q3_16 - q1_16))
upperfence17=q3_17 + (1.5*(q3_17 - q1_17))
upperfence18=q3_18 + (1.5*(q3_18 - q1_18))
upperfence19=q3_19 + (1.5*(q3_19 - q1_19))
upperfence20=q3_20 + (1.5*(q3_20 - q1_20))
upperfence21=q3_21 + (1.5*(q3_21 - q1_21))


#Replace upper outliers with upper fence for each age bin
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="5" & Pheno_Comb2$EACSF_cm >= upperfence5, upperfence5, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="6" & Pheno_Comb2$EACSF_cm >= upperfence6, upperfence6, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="7" & Pheno_Comb2$EACSF_cm >= upperfence7, upperfence7, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="8" & Pheno_Comb2$EACSF_cm >= upperfence8, upperfence8, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="9" & Pheno_Comb2$EACSF_cm >= upperfence9, upperfence9, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="10" & Pheno_Comb2$EACSF_cm >= upperfence10, upperfence10, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="11" & Pheno_Comb2$EACSF_cm >= upperfence11, upperfence11, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="12" & Pheno_Comb2$EACSF_cm >= upperfence12, upperfence12, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="13" & Pheno_Comb2$EACSF_cm >= upperfence13, upperfence13, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="14" & Pheno_Comb2$EACSF_cm >= upperfence14, upperfence14, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="15" & Pheno_Comb2$EACSF_cm >= upperfence15, upperfence15, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="16" & Pheno_Comb2$EACSF_cm >= upperfence16, upperfence16, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="17" & Pheno_Comb2$EACSF_cm >= upperfence17, upperfence17, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="18" & Pheno_Comb2$EACSF_cm >= upperfence18, upperfence18, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="19" & Pheno_Comb2$EACSF_cm >= upperfence19, upperfence19, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="20" & Pheno_Comb2$EACSF_cm >= upperfence20, upperfence20, Pheno_Comb2$EACSF_cm)
Pheno_Comb2$EACSF_cm <- ifelse(Pheno_Comb2$Age_bins=="21" & Pheno_Comb2$EACSF_cm >= upperfence21, upperfence21, Pheno_Comb2$EACSF_cm)


#Check for number of changes made for each age bin
table(Pheno_Comb2$EACSF_cm==upperfence5)
table(Pheno_Comb2$EACSF_cm==upperfence6)
table(Pheno_Comb2$EACSF_cm==upperfence7)
table(Pheno_Comb2$EACSF_cm==upperfence8)
table(Pheno_Comb2$EACSF_cm==upperfence9)
table(Pheno_Comb2$EACSF_cm==upperfence10)
table(Pheno_Comb2$EACSF_cm==upperfence11)
table(Pheno_Comb2$EACSF_cm==upperfence12)
table(Pheno_Comb2$EACSF_cm==upperfence13)
table(Pheno_Comb2$EACSF_cm==upperfence14)
table(Pheno_Comb2$EACSF_cm==upperfence15)
table(Pheno_Comb2$EACSF_cm==upperfence16)
table(Pheno_Comb2$EACSF_cm==upperfence17)
table(Pheno_Comb2$EACSF_cm==upperfence18)
table(Pheno_Comb2$EACSF_cm==upperfence19)
table(Pheno_Comb2$EACSF_cm==upperfence20)
table(Pheno_Comb2$EACSF_cm==upperfence21)

#--------------------------------------------------------TABLE 1-------------------------------------------------
#Number of participants per group
table(Pheno_Comb2$group)

#Sex by group (0 = male, 1 = female)
table(Pheno_Comb2$group, Pheno_Comb2$Sex)

#Age by groups
favstats(Age~group, data=Pheno_Comb2)
favstats(Pheno_Comb2$Age)

#LVCP vs. HVCP participants
Pheno_Comb2$WISC_FSIQ <- as.numeric(Pheno_Comb2$WISC_FSIQ) 
Pheno_Comb2$LVCP <- ifelse(Pheno_Comb2$WISC_FSIQ<=79, 1 ,0)
#1=LVCP, 0=HVCP
table(Pheno_Comb2$LVCP, Pheno_Comb2$group)

#FIQ
Pheno_Comb2$WISC_FSIQ <- as.numeric(Pheno_Comb2$WISC_FSIQ) 
favstats(WISC_FSIQ~group, data=Pheno_Comb2)

#SRS
Pheno_Comb2$SRS_Total <- as.numeric(Pheno_Comb2$SRS_Total) 
favstats(SRS_Total~group, data=Pheno_Comb2)

#SDS (sleep disturbance scale)
Pheno_Comb2$SDS_Total_Raw <- as.numeric(Pheno_Comb2$SDS_Total_Raw) 
favstats(SDS_Total_Raw~group, data=Pheno_Comb2)

#Table of Sites
table(Pheno_Comb2$Site)



#GROUP COMPARISONS 

  ## Chi-square test (biological sex)
  Pheno_CombX.data <- data.frame(Pheno_Comb2$Sex, Pheno_Comb2$group)
  # Create a table with the needed variables.
  Pheno_CombX.data = table(Pheno_Comb2$Sex, Pheno_Comb2$group) 
  # Perform the Chi-Square test.
  print(chisq.test(Pheno_CombX.data))
  #Multiple comparisons
  #chisq.multcomp(Pheno_CombX.data, p.method = "bonferroni")
    #ASD vs NoDX
    Pheno_CombASDND <- subset(Pheno_Comb2, group!="OtherDx") 
    Pheno_CombASDND$group <- factor(Pheno_CombASDND$group)
    Pheno_CombX.data <- data.frame(Pheno_CombASDND$Sex, Pheno_CombASDND$group)
    Pheno_CombX.data = table(Pheno_CombASDND$Sex, Pheno_CombASDND$group) 
    print(chisq.test(Pheno_CombX.data))
    
    #ASD vs OtherDx
    Pheno_CombASDOD <- subset(Pheno_Comb2, group!="NoDX")
    Pheno_CombASDOD$group <- factor(Pheno_CombASDOD$group)
    Pheno_CombX.data <- data.frame(Pheno_CombASDOD$Sex, Pheno_CombASDOD$group)
    Pheno_CombX.data = table(Pheno_CombASDOD$Sex, Pheno_CombASDOD$group) 
    print(chisq.test(Pheno_CombX.data))
    
    #NoDx vs OtherDx
    Pheno_CombNDOD <- subset(Pheno_Comb2, group!="ASD")
    Pheno_CombNDOD$group <- factor(Pheno_CombNDOD$group)
    Pheno_CombX.data <- data.frame(Pheno_CombNDOD$Sex, Pheno_CombNDOD$group)
    Pheno_CombX.data = table(Pheno_CombNDOD$Sex, Pheno_CombNDOD$group) 
    print(chisq.test(Pheno_CombX.data))
    
    
  ##ANOVA (Age)
  age_model <- aov(Age ~ group, data = Pheno_Comb2)
  summary(age_model)
  
  ##ANOVA (FSIQ)
  Pheno_Comb2$WISC_FSIQ <- as.numeric(Pheno_Comb2$WISC_FSIQ) 
  FSIQ_model <- aov(WISC_FSIQ ~ group, data = Pheno_Comb2)
  summary(FSIQ_model)
  pairwise.t.test(Pheno_Comb2$WISC_FSIQ, Pheno_Comb2$group, p.adjust.method="bonferroni")
  
  ##ANOVA (SRS)
  Pheno_Comb2$SRS_Total <- as.numeric(Pheno_Comb2$SRS_Total) 
  SRS_model <- aov(SRS_Total ~ group, data =Pheno_Comb2)
  summary(SRS_model)
  pairwise.t.test(Pheno_Comb2$SRS_Total, Pheno_Comb2$group, p.adjust.method="bonferroni")
  
  ##ANOVA (SDS)
  Pheno_Comb2$SDS_Total_Raw <- as.numeric(Pheno_Comb2$SDS_Total_Raw) 
  SDS_model <- aov(SDS_Total_Raw ~ group, data=Pheno_Comb2)
  summary(SDS_model)
  pairwise.t.test(Pheno_Comb2$SDS_Total_Raw, Pheno_Comb2$group, p.adjust.method="bonferroni")
  
  
#FIGURE 1 : Age violinplots by group
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")
ggplot(Pheno_Comb2, aes(x = group, y = Age, fill = group)) +
  labs(x = "", y = 'Age')+
  geom_violin() +
  geom_boxplot(width=0.1, color="black", position=position_dodge(.91)) +
  scale_colour_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
  scale_x_discrete(labels=c("Autism", "No Diagnosis", "Other Diagnosis")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.title=element_blank(), legend.text=element_blank(), 
        legend.background = element_blank(), axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#save the file
ggsave(filename = paste("Figure1_210520.png"), width = 4.4, height = 3.4,
       path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/manuscript", dpi = 300)



# FIGURE 1: Age Raincloud Plots by Group
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")
ggplot(Pheno_Comb2, aes(x = group, y = Age, fill=group)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = "NA") + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  coord_cartesian((xlim = c(1.2, NA)), (ylim=c(5, 22)), clip = "off")+
  labs(y="Age (Years)", x="Diagnostic Group")+
  scale_colour_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+  theme(axis.text=element_text(size = 9), axis.title = element_text(size = 9))+
  scale_x_discrete(labels=c("Autism", "No Diagnosis", "Other Diagnosis")) +
    theme(legend.position = "none")+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"))

#save the file
ggsave(filename = paste("Figure1_220802.pdf"), width = 3.4, height = 3.4,
       path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/manuscript", dpi = 300)





#Figure Legend (extract from this faux plot)
require(ggplot2)
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")
ggplot(Pheno_Comb2, aes(x=Age, y=Age, color=group))+
  geom_line(aes(group =group))+
  labs(x = "Age (Years)", y= "Age (Years)")+
  geom_point(aes(fill=group), colour="black", pch=21)+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
     scale_y_continuous(limits=c(5,22))+
     scale_x_continuous(limits=c(5,22))+
     theme_bw()+
     theme(panel.grid.major = element_line(color = "grey"))+
     theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
                     axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     legend.position = c(0.8, 0.2), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
                     legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
                     axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#save the file
ggsave(filename = paste("FigureLegend_220707.svg"), width = 3.4, height = 3.4,
       path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/manuscript", dpi = 300)




#Figure 3: Age by CSF (linear model)
require(ggplot2)
#Using all three groups
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")
ggplot(Pheno_Comb2, aes(x=Age, y=EACSF_cm, color=group))+
  labs(x = "Age (Years)", y = 'EA-CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=lm, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(20,150))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.95, .2), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#save the file
ggsave(filename = paste("Figure3_220707.svg"), width = 3.4, height = 3.4,
       path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/manuscript", dpi = 300)


#--------------------------------------------------VISUAL QC ANALYSIS--------------------------------------------------

#QC ICC analysis
#Compare ICC for MID02 (EA-CSF segmentation) ratings prior to resolving discrepancies
QC$Notes <- NULL
QC$EID <- NULL
QC$X <- NULL
QC$T1w_QC_Anne <- NULL
QC$T1w_QC_Christopher <- NULL
QC$T1w_QC_Maddy <- NULL
QC$T1w_QC_Matt <- NULL
QC$Avg_Rating <- NULL
#ICC overall is 0.95
ICC(QC)


#Dataset with all QC ratings + Pheno data
#Drop cases without a T1w image
Comb5 <- subset(Comb11, T1!="NA")
#Drop duplicates
Comb5$duplicate <- !duplicated(Comb5$EID)
Comb5 <- subset(Comb5, duplicate==TRUE)


#Create a group var and count EA-CSF seg fail by group

#Create var for any ASD diagnosis
Comb5$ASDany <- ifelse(Comb5$DX_01 == "Autism Spectrum Disorder" | Comb5$DX_02 == "Autism Spectrum Disorder" | Comb5$DX_03 == "Autism Spectrum Disorder" | Comb5$DX_04 == "Autism Spectrum Disorder" | Comb5$DX_05 == "Autism Spectrum Disorder" | Comb5$DX_06 == "Autism Spectrum Disorder" | Comb5$DX_07 == "Autism Spectrum Disorder" | Comb5$DX_08 == "Autism Spectrum Disorder" | Comb5$DX_09 == "Autism Spectrum Disorder" | Comb5$DX_10 == "Autism Spectrum Disorder", "ASD", NA)
#Create var for no dx
Comb5$NoDX <- ifelse(Comb5$DX_01 == "No Diagnosis Given" | Comb5$DX_02 == "No Diagnosis Given" | Comb5$DX_03 == "No Diagnosis Given" | Comb5$DX_04 == "No Diagnosis Given" | Comb5$DX_05 == "No Diagnosis Given" | Comb5$DX_06 == "No Diagnosis Given" | Comb5$DX_07 == "No Diagnosis Given" | Comb5$DX_08 == "No Diagnosis Given" | Comb5$DX_09 == "No Diagnosis Given" | Comb5$DX_10 == "No Diagnosis Given", "NoDX", NA)
#Create var for other Dx
Comb5$NotASDdx <- ifelse(is.na(Comb5$ASDany) & is.na(Comb5$NoDX), "OtherDx", " ")

#Create var for any ASD diagnosis
Comb5$ASDany <- ifelse(Comb5$DX_01 == "Autism Spectrum Disorder" | Comb5$DX_02 == "Autism Spectrum Disorder" | Comb5$DX_03 == "Autism Spectrum Disorder" | Comb5$DX_04 == "Autism Spectrum Disorder" | Comb5$DX_05 == "Autism Spectrum Disorder" | Comb5$DX_06 == "Autism Spectrum Disorder" | Comb5$DX_07 == "Autism Spectrum Disorder" | Comb5$DX_08 == "Autism Spectrum Disorder" | Comb5$DX_09 == "Autism Spectrum Disorder" | Comb5$DX_10 == "Autism Spectrum Disorder", "ASD", " ")
#Create var for no dx
Comb5$NoDX <- ifelse(Comb5$DX_01 == "No Diagnosis Given" | Comb5$DX_02 == "No Diagnosis Given" | Comb5$DX_03 == "No Diagnosis Given" | Comb5$DX_04 == "No Diagnosis Given" | Comb5$DX_05 == "No Diagnosis Given" | Comb5$DX_06 == "No Diagnosis Given" | Comb5$DX_07 == "No Diagnosis Given" | Comb5$DX_08 == "No Diagnosis Given" | Comb5$DX_09 == "No Diagnosis Given" | Comb5$DX_10 == "No Diagnosis Given", "NoDX", " ")

#Create overall group var
Comb5$group <- (paste0(Comb5$ASDany, Comb5$NoDX, Comb5$NotASDdx))
#Replace NANAOtherDx with OtherDx
Comb5$group <- str_replace(Comb5$group, "NANA", "")
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
Comb5$group <- trim(Comb5$group)

#number of scans in each group
table(Comb5$group)

#table of avg_Rating EA-CSF seg by group
table(Comb5$Avg_Rating<"2", Comb5$group)
table(Comb5$Avg_Rating<"2")

#T1w ratings by group
#Convert to numerical values (run trim lines first!!)
Comb5$T1w_QC_Anne <- gsub("^.{0,1}", "", Comb5$T1w_QC_Anne)
Comb5$T1w_QC_Matt <- gsub("^.{0,1}", "", Comb5$T1w_QC_Matt)
Comb5$T1w_QC_Christopher <- gsub("^.{0,1}", "", Comb5$T1w_QC_Christopher)
Comb5$T1w_QC_Maddy <- gsub("^.{0,1}", "", Comb5$T1w_QC_Maddy)

#Avg across rows
Comb5$T1w_QC_Matt[ Comb5$T1w_QC_Matt == ""] <- NA
Comb5$T1w_QC_Anne[ Comb5$T1w_QC_Anne == ""] <- NA
Comb5$T1w_QC_Christopher[ Comb5$T1w_QC_Christopher == ""] <- NA
Comb5$T1w_QC_Maddy[ Comb5$T1w_QC_Maddy == ""] <- NA

qc_t1 <- data.frame(Comb5$T1w_QC_Anne, Comb5$T1w_QC_Matt, Comb5$T1w_QC_Christopher, Comb5$T1w_QC_Maddy)
qc_t1$Comb5.T1w_QC_Anne <-as.numeric(qc_t1$Comb5.T1w_QC_Anne)
qc_t1$Comb5.T1w_QC_Matt <-as.numeric(qc_t1$Comb5.T1w_QC_Matt)
qc_t1$Comb5.T1w_QC_Christopher <-as.numeric(qc_t1$Comb5.T1w_QC_Christopher)
qc_t1$Comb5.T1w_QC_Maddy <-as.numeric(qc_t1$Comb5.T1w_QC_Maddy)

Comb5$T1w_Avg_Rating <- rowMeans(qc_t1, na.rm=TRUE)

#Table of T1w rating by group (T1w fail = 3)
table(Comb5$T1w_Avg_Rating, Comb5$group)


#Create overall fail var
Comb5$T1w_pass <- ifelse(Comb5$T1w_Avg_Rating<3, 1, 0)
table(Comb5$T1w_pass)

#Number of subjects that pass both T1w and EA-CSF seg QC (pass overlap = 727)
Comb5$qc_pass <- ifelse(Comb5$T1w_Avg_Rating<3 & Comb5$Avg_Rating<2, 1, 0) 
table(Comb5$qc_pass)
table(Comb5$qc_pass, Comb5$group)

#Number of subjects that fail both T1w and EA-CSF (fail overlap = 439)
Comb5$qc_overlap <- ifelse(Comb5$T1w_Avg_Rating==3 & Comb5$Avg_Rating>=2, 1, 0)
table(Comb5$qc_overlap)



#Is QC rating corr. with SRS?
#Drop subjects without SRS_Total
Comb5$SRS_Total <- as.numeric(Comb5$SRS_Total) 
CombSRS <- subset(Comb5, SRS_Total!="NA")
cor.test(CombSRS$Avg_Rating, CombSRS$SRS_Total)


#Is QC rating corr. with FSIQ?
Comb5$WISC_FSIQ <- as.numeric(Comb5$WISC_FSIQ) 
cor.test(Comb5$WISC_FSIQ, Comb5$Avg_Rating)


#Is there a group difference in EA-CSF ratings?
group.aov<-aov(Avg_Rating~group, data=Comb5) 
summary(group.aov)



#---------------------------------------------------ASSUMPTIONS TESTING------------------------------------------------

#Linearity assumption -- regression lines on graph 
ggscatter(
  Pheno_Comb2, x = "Age", y = "EACSF_cm",
  color = "group", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = group)
  )


#Homogeneity of Regression Slopes
##Tests that there is no sig. relationship between group and age
Pheno_Comb2 %>% anova_test(EACSF_cm ~ group*mc_age)
#group:Age term was not sig. so regression slopes look good!


#Assess normality of residuals - Shapiro-Wilk test
#1. Fit the model, the covariate goes first
model <- lm(EACSF_cm ~ mc_age + Site + Sex + eTIV_mc + eTIV_squared_mc + mc_age*group + mc_age_squared + mc_age_squared*group + group, data = Pheno_Comb2)
#2. Run model metrics
model.metrics <- augment(model)
head(model.metrics, 3) #this line just displays the metrics
#3. Run the test
shapiro_test(model.metrics$.resid)
#If p-value was significant, reject null hypothesis that residuals are normal


#Homogeneity of variance - Levene's test
model.metrics %>% levene_test(.resid ~ group)
#If p-value was significant, reject null hypothesis of homogeneity of variance for all groups

#--------------------------------MAIN ANCOVA ANALYSIS----------------------------------------------

#ANCOVA for overall group effects

#And model including eTIV_cm_squared covariate (since in the graph of eTIV by EA-CSF, it appears to have a nonlinear growth pattern)
eTIVmodel2 = lm(EACSF_cm ~ group + mc_age + Site + Sex + eTIV_mc + eTIV_squared_mc + mc_age*group + mc_age_squared + mc_age_squared*group, data = Pheno_Comb2)
#Run ANOVA on model
Anova(eTIVmodel2, type="III")  



#Eta-squared (for effect sizes)
etaSquared(eTIVmodel2, type=3)


#Parse out group differences if omnibus had sig. group effect

#Pairwise (Bonferroni) corrections using emmeans
##emmeans = estimated marginal means 
## = The estimated marginal means output gives the adjusted means (controlling for the covariates) for each group.
#1. Get marginal means
#marginal = emmeans(lmodel, ~ group + Age + Site + Sex + Age*group + Age_squared + Age_squared*group)
marginal = emmeans(eTIVmodel2, ~ group*mc_age + mc_age_squared*group)
#Display marginal means table
marginal
#2. Use Bonferroni corrections on marginal means, run pairwise corrections
pairs(marginal, adjust="bon") 






#ANCOVA Models with age re-centered 
#5
Pheno_Comb2$age5 <- Pheno_Comb2$Age - 5
Pheno_Comb2$age_squared5 <- Pheno_Comb2$age5 * Pheno_Comb2$age5
agemodel5 = lm(EACSF_cm ~ group + age5 + Site + Sex + age5*group + age_squared5 + age_squared5*group + eTIV_mc + eTIV_squared_mc, data = Pheno_Comb2)
Anova(agemodel5, type="III")

#10
Pheno_Comb2$age10 <- Pheno_Comb2$Age - 10
Pheno_Comb2$age_squared10 <- Pheno_Comb2$age10 * Pheno_Comb2$age10
agemodel10 = lm(EACSF_cm ~ group + age10 + Site + Sex + age10*group + age_squared10 + age_squared10*group + eTIV_mc + eTIV_squared_mc, data = Pheno_Comb2)
Anova(agemodel10, type="III")

#15
Pheno_Comb2$age15 <- Pheno_Comb2$Age - 15
Pheno_Comb2$age_squared15 <- Pheno_Comb2$age15 * Pheno_Comb2$age15
agemodel15 = lm(EACSF_cm ~ group + age15 + Site + Sex + age15*group + age_squared15 + age_squared15*group + eTIV_mc + eTIV_squared_mc, data = Pheno_Comb2)
Anova(agemodel15, type="III")

#20
Pheno_Comb2$age20 <- Pheno_Comb2$Age - 20
Pheno_Comb2$age_squared20 <- Pheno_Comb2$age20 * Pheno_Comb2$age20
agemodel20 = lm(EACSF_cm ~ group + age20 + Site + Sex + age20*group + age_squared20 + age_squared20*group + eTIV_mc + eTIV_squared_mc, data = Pheno_Comb2)
Anova(agemodel20, type="III")


#----------------------------ANCOVA: Dx GROUPS BY ASD---------------------------

#Compare groups that survived the cutoff against ASD: ADHD Combined, ADHD Inattentive, Generalized Anxiety Disorder, Social Anxiety, SLD Reading

#ADHD Combined Type vs ASD Analysis
  #Subset to ASD and ADHD Combined Type
  Pheno_Comb2$KeepDx1 <- ifelse(Pheno_Comb2$ASD_dummy == "1" | Pheno_Comb2$ADHDCombinedType_dummy=="1", 1, 0)
  #Filter out ASD subjects with comorbid ADHD_Comb diagnosis
  Pheno_Comb2$KeepDx2 <- ifelse(Pheno_Comb2$ASD_dummy == "1" & Pheno_Comb2$ADHDCombinedType_dummy=="1", 0, 1)
  ADHD_Comb <- subset(Pheno_Comb2, KeepDx1!="0")
  ADHD_Comb <- subset(ADHD_Comb, KeepDx2!="0")
  
  #Create group var
  #Create var for any ASD diagnosis
  ADHD_Comb$ASDany <- ifelse(ADHD_Comb$DX_01 == "Autism Spectrum Disorder" | ADHD_Comb$DX_02 == "Autism Spectrum Disorder" | ADHD_Comb$DX_03 == "Autism Spectrum Disorder" | ADHD_Comb$DX_04 == "Autism Spectrum Disorder" | ADHD_Comb$DX_05 == "Autism Spectrum Disorder" | ADHD_Comb$DX_06 == "Autism Spectrum Disorder" | ADHD_Comb$DX_07 == "Autism Spectrum Disorder" | ADHD_Comb$DX_08 == "Autism Spectrum Disorder" | ADHD_Comb$DX_09 == "Autism Spectrum Disorder" | ADHD_Comb$DX_10 == "Autism Spectrum Disorder", "ASD", "NA")
  #Create var for any ADHD Combined dx
  ADHD_Comb$ADHDComb <- ifelse(is.na(ADHD_Comb$ASDany), "ADHDComb", " ")
  #Create overall group var
    ADHD_Comb$group <- (paste0(ADHD_Comb$ASDany, ADHD_Comb$ADHDComb))
    #Replace NAADHDComb with ADHDComb
    ADHD_Comb$group <- str_replace(ADHD_Comb$group, "NA", "ADHDComb")
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    ADHD_Comb$group <- trim(ADHD_Comb$group)
    ADHD_Comb$group <- as.factor(ADHD_Comb$group)
  
  #Create new mean-age and eTIV vars for this group  
    meanage = mean(ADHD_Comb$Age)
    ADHD_Comb$mc_age <- ADHD_Comb$Age - meanage
    ADHD_Comb$mc_age_squared <- ADHD_Comb$mc_age * ADHD_Comb$mc_age
    ADHD_Comb$mc_age_squared <- as.numeric(ADHD_Comb$mc_age_squared)
    
  #Make etiv and etiv-squared mean-centered
    etiv_mean = mean(ADHD_Comb$eTIV_cm)
    ADHD_Comb$eTIV_mc <- ADHD_Comb$eTIV_cm - etiv_mean
    ADHD_Comb$eTIV_squared_mc <- ADHD_Comb$eTIV_mc * ADHD_Comb$eTIV_mc
    
  #Frequency Stats
    table(ADHD_Comb$group)
    favstats(data=ADHD_Comb, Age~group)
    
  #ANCOVA: ADHDComb and ASD
    adhdcomb_model = lm(EACSF_cm ~ group + mc_age + Site + Sex + mc_age*group + mc_age_squared + mc_age_squared*group + eTIV_mc + eTIV_squared_mc, data = ADHD_Comb)
    summary(adhdcomb_model)

#Plot for ADHDComb and ASD EA-CSF over Age
    require(ggplot2)
    #Plot
    Palette <- c("#E69F00", "#00BCB8", "#0072B2", "#0072B2")
    ggplot(ADHD_Comb, aes(x=Age, y=EACSF_cm, color=group))+
      labs(x = "Age (Years)", y = 'EA-CSF Volume ('~cm^3*')')+
      labs(fill = " ")+
      labs(color = " ")+
      scale_colour_manual(values=Palette, labels = c("ADHD Combined Type", "Autism"))+
      scale_fill_manual(values=Palette, labels = c("ADHD Combined Type", "Autism"))+
      geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
      geom_smooth(method=lm, aes(fill=group), se=TRUE) +
      scale_y_continuous(limits=c(20,150))+
      scale_x_continuous(limits=c(5,22), expand = c(0,0))+
      guides(color=guide_legend(override.aes=list(fill=NA)))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = c(.95, .2), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
            legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())
    
      ggsave(filename = paste("ADHDComb_220715.svg"), width = 3.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    
      #FIGURE: Age violinplots by group
      Palette <- c("#E69F00", "#00BCB8", "#0072B2", "#0072B2")
      ggplot(ADHD_Comb, aes(x = group, y = Age, fill = group)) +
        labs(x = "", y = 'Age')+
        geom_violin() +
        geom_boxplot(width=0.1, color="black", position=position_dodge(.91)) +
        scale_colour_manual(values=Palette, labels = c("ADHDComb", "Autism", "No Diagnosis"))+
        scale_fill_manual(values=Palette, labels = c("ADHDComb", "Autism", "No Diagnosis"))+
        scale_x_discrete(labels=c("ADHDComb", "Autism", "No Diagnosis")) +
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
              axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              legend.title=element_blank(), legend.text=element_blank(), 
              legend.background = element_blank(), axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
              axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())+
        theme(legend.position = "none")
      
      #save the file
      ggsave(filename = paste("ADHDComb_AgeViolin_220715.svg"), width = 4.4, height = 3.4,
             path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
      
    
    

#ADHD Inattentive Type vs ASD Analysis
    #Subset to ASD and ADHD Inattentive Type
    Pheno_Comb2$KeepDx1 <- ifelse(Pheno_Comb2$ASD_dummy == "1" | Pheno_Comb2$ADHDInattentiveType_dummy=="1", 1, 0)
    Pheno_Comb2$KeepDx2 <- ifelse(Pheno_Comb2$ASD_dummy == "1" & Pheno_Comb2$ADHDInattentiveType_dummy=="1", 0, 1)
    ADHD_IT <- subset(Pheno_Comb2, KeepDx1!="0")
    ADHD_IT <- subset(ADHD_IT, KeepDx2!="0")
    
    #Create group var
    #Create var for any ASD diagnosis
    ADHD_IT$ASDany <- ifelse(ADHD_IT$DX_01 == "Autism Spectrum Disorder" | ADHD_IT$DX_02 == "Autism Spectrum Disorder" | ADHD_IT$DX_03 == "Autism Spectrum Disorder" | ADHD_IT$DX_04 == "Autism Spectrum Disorder" | ADHD_IT$DX_05 == "Autism Spectrum Disorder" | ADHD_IT$DX_06 == "Autism Spectrum Disorder" | ADHD_IT$DX_07 == "Autism Spectrum Disorder" | ADHD_IT$DX_08 == "Autism Spectrum Disorder" | ADHD_IT$DX_09 == "Autism Spectrum Disorder" | ADHD_IT$DX_10 == "Autism Spectrum Disorder", "ASD", "NA")
    #Create var for any ADHD Inattentive dx
    ADHD_IT$ADHDIT <- ifelse(is.na(ADHD_IT$ASDany), "ADHDIT", " ")
    #Create overall group var
    ADHD_IT$group <- (paste0(ADHD_IT$ASDany, ADHD_IT$ADHDIT))
    #Replace NA with ADHDIT
    ADHD_IT$group <- str_replace(ADHD_IT$group, "NA", "ADHDIT")
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    ADHD_IT$group <- trim(ADHD_IT$group)
    ADHD_IT$group <- as.factor(ADHD_IT$group)
    
    #Create new mean-age and eTIV vars for this group  
    meanage = mean(ADHD_IT$Age)
    ADHD_IT$mc_age <- ADHD_IT$Age - meanage
    ADHD_IT$mc_age_squared <- ADHD_IT$mc_age * ADHD_IT$mc_age
    ADHD_IT$mc_age_squared <- as.numeric(ADHD_IT$mc_age_squared)
    
    #Make etiv and etiv-squared mean-centered
    etiv_mean = mean(ADHD_IT$eTIV_cm)
    ADHD_IT$eTIV_mc <- ADHD_IT$eTIV_cm - etiv_mean
    ADHD_IT$eTIV_squared_mc <- ADHD_IT$eTIV_mc * ADHD_IT$eTIV_mc
    
    #Frequency Stats
    table(ADHD_IT$group)
    favstats(data=ADHD_IT, Age~group)
    
    #ANCOVA: ADHD IT and ASD
    adhdit_model = lm(EACSF_cm ~ group + mc_age + Site + Sex + mc_age*group + mc_age_squared + mc_age_squared*group + eTIV_mc + eTIV_squared_mc, data = ADHD_IT)
    summary(adhdit_model)
    
    #Plot for ADHDIT and ASD EA-CSF over Age
    require(ggplot2)
    #Plot
    Palette <- c("#E69F00", "#00BCB8", "#0072B2", "#0072B2")
    ggplot(ADHD_IT, aes(x=Age, y=EACSF_cm, color=group))+
      labs(x = "Age (Years)", y = 'EA-CSF Volume ('~cm^3*')')+
      labs(fill = " ")+
      labs(color = " ")+
      scale_colour_manual(values=Palette, labels = c("ADHD Inattentive Type", "Autism", "No Diagnosis"))+
      scale_fill_manual(values=Palette, labels = c("ADHD Inattentive Type", "Autism", "No Diagnosis"))+
      geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
      geom_smooth(method=lm, aes(fill=group), se=TRUE) +
      scale_y_continuous(limits=c(20,150))+
      scale_x_continuous(limits=c(5,22), expand = c(0,0))+
      guides(color=guide_legend(override.aes=list(fill=NA)))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = c(.95, .2), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
            legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())
    
    ggsave(filename = paste("ADHDIT_220715.svg"), width = 3.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    #FIGURE: Age violinplots by group
    Palette <- c("#E69F00", "#00BCB8", "#0072B2", "#0072B2")
    ggplot(ADHD_IT, aes(x = group, y = Age, fill = group)) +
      labs(x = "", y = 'Age')+
      geom_violin() +
      geom_boxplot(width=0.1, color="black", position=position_dodge(.91)) +
      scale_colour_manual(values=Palette, labels = c("ADHD IT", "Autism", "No Diagnosis"))+
      scale_fill_manual(values=Palette, labels = c("ADHD IT", "Autism", "No Diagnosis"))+
      scale_x_discrete(labels=c("ADHD IT", "Autism", "No Diagnosis")) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.title=element_blank(), legend.text=element_blank(), 
            legend.background = element_blank(), axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())+
      theme(legend.position = "none")
    
    #save the file
    ggsave(filename = paste("ADHDIT_AgeViolin_220715.svg"), width = 4.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    

#Generalized Anxiety Disorder (GAD) vs ASD Analysis
    #Subset to ASD and GAD
    Pheno_Comb2$KeepDx1 <- ifelse(Pheno_Comb2$ASD_dummy == "1" | Pheno_Comb2$GeneralizedAnxietyDisorder_dummy=="1", 1, 0)
    Pheno_Comb2$KeepDx2 <- ifelse(Pheno_Comb2$ASD_dummy == "1" & Pheno_Comb2$GeneralizedAnxietyDisorder_dummy=="1", 0, 1)
    GAD <- subset(Pheno_Comb2, KeepDx1!="0")
    GAD <- subset(GAD, KeepDx2!="0")
    
    #Create group var
    #Create var for any ASD diagnosis
    GAD$ASDany <- ifelse(GAD$DX_01 == "Autism Spectrum Disorder" | GAD$DX_02 == "Autism Spectrum Disorder" | GAD$DX_03 == "Autism Spectrum Disorder" | GAD$DX_04 == "Autism Spectrum Disorder" | GAD$DX_05 == "Autism Spectrum Disorder" | GAD$DX_06 == "Autism Spectrum Disorder" | GAD$DX_07 == "Autism Spectrum Disorder" | GAD$DX_08 == "Autism Spectrum Disorder" | GAD$DX_09 == "Autism Spectrum Disorder" | GAD$DX_10 == "Autism Spectrum Disorder", "ASD", "NA")
    #Create var for any ADHD Combined dx
    GAD$GAD <- ifelse(is.na(GAD$ASDany), "GAD", " ")
    #Create overall group var
    GAD$group <- (paste0(GAD$ASDany, GAD$NoDX, GAD$GAD))
    #Replace NA with GAD
    GAD$group <- str_replace(GAD$group, "NA", "GAD")
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    GAD$group <- trim(GAD$group)
    GAD$group <- as.factor(GAD$group)
    
    #Create new mean-age and eTIV vars for this group  
    meanage = mean(GAD$Age)
    GAD$mc_age <- GAD$Age - meanage
    GAD$mc_age_squared <- GAD$mc_age * GAD$mc_age
    GAD$mc_age_squared <- as.numeric(GAD$mc_age_squared)
    
    #Make etiv and etiv-squared mean-centered
    etiv_mean = mean(GAD$eTIV_cm)
    GAD$eTIV_mc <- GAD$eTIV_cm - etiv_mean
    GAD$eTIV_squared_mc <- GAD$eTIV_mc * GAD$eTIV_mc
    
    #Frequency Stats
    table(GAD$group)
    favstats(data=GAD, Age~group)
    
    #ANCOVA: ADHDComb and ASD
    gad_model = lm(EACSF_cm ~ group + mc_age + Site + Sex + mc_age*group + mc_age_squared + mc_age_squared*group + eTIV_mc + eTIV_squared_mc, data = GAD)
    summary(gad_model)
   
    
    #Plot for GAD and ASD EA-CSF over Age
    require(ggplot2)
    #Plot
    Palette <- c("#00BCB8", "#E69F00", "#0072B2", "#0072B2")
    ggplot(GAD, aes(x=Age, y=EACSF_cm, color=group))+
      labs(x = "Age (Years)", y = 'EA-CSF Volume ('~cm^3*')')+
      labs(fill = " ")+
      labs(color = " ")+
      scale_colour_manual(values=Palette, labels = c("Autism", "GAD", "No Diagnosis"))+
      scale_fill_manual(values=Palette, labels = c("Autism","GAD", "No Diagnosis"))+
      geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
      geom_smooth(method=lm, aes(fill=group), se=TRUE) +
      scale_y_continuous(limits=c(20,150))+
      scale_x_continuous(limits=c(5,22), expand = c(0,0))+
      guides(color=guide_legend(override.aes=list(fill=NA)))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = c(.95, .2), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
            legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())
    
    ggsave(filename = paste("GAD_220715.svg"), width = 3.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    
    #FIGURE: Age violinplots by group
    Palette <- c("#00BCB8", "#E69F00", "#0072B2", "#0072B2")
    ggplot(GAD, aes(x = group, y = Age, fill = group)) +
      labs(x = "", y = 'Age')+
      geom_violin() +
      geom_boxplot(width=0.1, color="black", position=position_dodge(.91)) +
      scale_colour_manual(values=Palette, labels = c("Autism", "GAD", "No Diagnosis"))+
      scale_fill_manual(values=Palette, labels = c("Autism", "GAD", "No Diagnosis"))+
      scale_x_discrete(labels=c("Autism", "GAD", "No Diagnosis")) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.title=element_blank(), legend.text=element_blank(), 
            legend.background = element_blank(), axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())+
      theme(legend.position = "none")
    
    #save the file
    ggsave(filename = paste("GAD_AgeViolin_220715.svg"), width = 4.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    
    
#Social Anxiety (SA) vs ASD Analysis
    #Subset to ASD and ADHD Combined Type
    Pheno_Comb2$KeepDx1 <- ifelse(Pheno_Comb2$ASD_dummy == "1" | Pheno_Comb2$SocialAnxiety_dummy=="1", 1, 0)
    Pheno_Comb2$KeepDx2 <- ifelse(Pheno_Comb2$ASD_dummy == "1" & Pheno_Comb2$SocialAnxiety_dummy=="1", 0, 1)
    SA <- subset(Pheno_Comb2, KeepDx1!="0")
    SA <- subset(SA, KeepDx2!="0")
    
    #Create group var
    #Create var for any ASD diagnosis
    SA$ASDany <- ifelse(SA$DX_01 == "Autism Spectrum Disorder" | SA$DX_02 == "Autism Spectrum Disorder" | SA$DX_03 == "Autism Spectrum Disorder" | SA$DX_04 == "Autism Spectrum Disorder" | SA$DX_05 == "Autism Spectrum Disorder" | SA$DX_06 == "Autism Spectrum Disorder" | SA$DX_07 == "Autism Spectrum Disorder" | SA$DX_08 == "Autism Spectrum Disorder" | SA$DX_09 == "Autism Spectrum Disorder" | SA$DX_10 == "Autism Spectrum Disorder", "ASD", "NA")
    #Create var for any ADHD Combined dx
    SA$SA <- ifelse(is.na(SA$ASDany), "SA", " ")
    #Create overall group var
    SA$group <- (paste0(SA$ASDany, SA$SA))
    #Replace NA with SA
    SA$group <- str_replace(SA$group, "NA", "SA")
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    SA$group <- trim(SA$group)
    SA$group <- as.factor(SA$group)
    
    #Frequency Stats
    table(SA$group)
    favstats(data=SA, Age~group)
    
    #Create new mean-age and eTIV vars for this group  
    meanage = mean(SA$Age)
    SA$mc_age <- SA$Age - meanage
    SA$mc_age_squared <- SA$mc_age * SA$mc_age
    SA$mc_age_squared <- as.numeric(SA$mc_age_squared)
    
    #Make etiv and etiv-squared mean-centered
    etiv_mean = mean(SA$eTIV_cm)
    SA$eTIV_mc <- SA$eTIV_cm - etiv_mean
    SA$eTIV_squared_mc <- SA$eTIV_mc * SA$eTIV_mc
    
    #ANCOVA: SA and ASD
    sa_model = lm(EACSF_cm ~ group + mc_age + Site + Sex + mc_age*group + mc_age_squared + mc_age_squared*group + eTIV_mc + eTIV_squared_mc, data = SA)
    summary(sa_model)
    
    #Plot for SA and ASD EA-CSF over Age
    require(ggplot2)
    #Plot
    Palette <- c("#00BCB8", "#E69F00", "#0072B2")
    ggplot(SA, aes(x=Age, y=EACSF_cm, color=group))+
      labs(x = "Age (Years)", y = 'EA-CSF Volume ('~cm^3*')')+
      labs(fill = " ")+
      labs(color = " ")+
      scale_colour_manual(values=Palette, labels = c("Autism", "SA"))+
      scale_fill_manual(values=Palette, labels = c("Autism", "SA"))+
      geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
      geom_smooth(method=lm, aes(fill=group), se=TRUE) +
      scale_y_continuous(limits=c(20,150))+
      scale_x_continuous(limits=c(5,22), expand = c(0,0))+
      guides(color=guide_legend(override.aes=list(fill=NA)))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = c(.95, .2), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
            legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())
    
    ggsave(filename = paste("SA_220715.svg"), width = 3.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    
    #FIGURE: Age violinplots by group
    Palette <- c("#00BCB8", "#E69F00", "#0072B2")
    ggplot(SA, aes(x = group, y = Age, fill = group)) +
      labs(x = "", y = 'Age')+
      geom_violin() +
      geom_boxplot(width=0.1, color="black", position=position_dodge(.91)) +
      scale_colour_manual(values=Palette, labels = c("Autism", "Social Anxiety"))+
      scale_fill_manual(values=Palette, labels = c("Autism", "Social Anxiety"))+
      scale_x_discrete(labels=c("Autism", "Social Anxiety")) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.title=element_blank(), legend.text=element_blank(), 
            legend.background = element_blank(), axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())+
      theme(legend.position = "none")
    
    #save the file
    ggsave(filename = paste("SA_AgeViolin_220715.svg"), width = 4.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    
    
#Specific Learning Disorder-Reading (SLDR) vs ASD Analysis
    #Subset to ASD and SLDR
    Pheno_Comb2$KeepDx1 <- ifelse(Pheno_Comb2$ASD_dummy == "1" | Pheno_Comb2$SLDReading_dummy=="1", 1, 0)
    Pheno_Comb2$KeepDx2 <- ifelse(Pheno_Comb2$ASD_dummy == "1" & Pheno_Comb2$SLDReading_dummy=="1", 0, 1)
    SLDR <- subset(Pheno_Comb2, KeepDx1!="0")
    SLDR <- subset(SLDR, KeepDx2!="0")
    
    #Create group var
    #Create var for any ASD diagnosis
    SLDR$ASDany <- ifelse(SLDR$DX_01 == "Autism Spectrum Disorder" | SLDR$DX_02 == "Autism Spectrum Disorder" | SLDR$DX_03 == "Autism Spectrum Disorder" | SLDR$DX_04 == "Autism Spectrum Disorder" | SLDR$DX_05 == "Autism Spectrum Disorder" | SLDR$DX_06 == "Autism Spectrum Disorder" | SLDR$DX_07 == "Autism Spectrum Disorder" | SLDR$DX_08 == "Autism Spectrum Disorder" | SLDR$DX_09 == "Autism Spectrum Disorder" | SLDR$DX_10 == "Autism Spectrum Disorder", "ASD", "NA")
    #Create var for any SLDR dx
    SLDR$SLDR <- ifelse(is.na(SLDR$ASDany), "SLDR", " ")
    #Create overall group var
    SLDR$group <- (paste0(SLDR$ASDany, SLDR$SLDR))
    #Replace NA with SLDR
    SLDR$group <- str_replace(SLDR$group, "NA", "SLDR")
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    SLDR$group <- trim(SLDR$group)
    SLDR$group <- as.factor(SLDR$group)
    
    #Frequency Stats
    table(SLDR$group)
    favstats(data=SLDR, Age~group)
    
    #Create new mean-age and eTIV vars for this group  
    meanage = mean(SLDR$Age)
    SLDR$mc_age <- SLDR$Age - meanage
    SLDR$mc_age_squared <- SLDR$mc_age * SLDR$mc_age
    SLDR$mc_age_squared <- as.numeric(SLDR$mc_age_squared)
    
    #Make etiv and etiv-squared mean-centered
    etiv_mean = mean(SLDR$eTIV_cm)
    SLDR$eTIV_mc <- SLDR$eTIV_cm - etiv_mean
    SLDR$eTIV_squared_mc <- SLDR$eTIV_mc * SLDR$eTIV_mc
    
    #ANCOVA: SLDR and ASD
    sldr_model = lm(EACSF_cm ~ group + mc_age + Site + Sex + mc_age*group + mc_age_squared + mc_age_squared*group + eTIV_mc + eTIV_squared_mc, data = SLDR)
    summary(sldr_model)
    
   
    #Plot for SLDR and ASD EA-CSF over Age
    require(ggplot2)
    #Plot
    Palette <- c("#00BCB8", "#E69F00", "#0072B2")
    ggplot(SLDR, aes(x=Age, y=EACSF_cm, color=group))+
      labs(x = "Age (Years)", y = 'EA-CSF Volume ('~cm^3*')')+
      labs(fill = " ")+
      labs(color = " ")+
      scale_colour_manual(values=Palette, labels = c("Autism", "SLDR"))+
      scale_fill_manual(values=Palette, labels = c("Autism", "SLDR"))+
      geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
      geom_smooth(method=lm, aes(fill=group), se=TRUE) +
      scale_y_continuous(limits=c(20,150))+
      scale_x_continuous(limits=c(5,22), expand = c(0,0))+
      guides(color=guide_legend(override.aes=list(fill=NA)))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = c(.95, .2), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
            legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())
    
    ggsave(filename = paste("SLDR_220715.svg"), width = 3.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    
    #FIGURE: Age violinplots by group
    Palette <- c("#00BCB8", "#E69F00", "#0072B2")
    ggplot(SLDR, aes(x = group, y = Age, fill = group)) +
      labs(x = "", y = 'Age')+
      geom_violin() +
      geom_boxplot(width=0.1, color="black", position=position_dodge(.91)) +
      scale_colour_manual(values=Palette, labels = c("Autism", "SLD-Reading"))+
      scale_fill_manual(values=Palette, labels = c("Autism", "SLD-Reading"))+
      scale_x_discrete(labels=c("Autism", "SLD-Reading")) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.title=element_blank(), legend.text=element_blank(), 
            legend.background = element_blank(), axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())+
      theme(legend.position = "none")
    
    #save the file
    ggsave(filename = paste("SLDR_AgeViolin_220715.svg"), width = 4.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    


#NoDx vs ASD (No Comorbidities) Analysis
    #Subset to ASD and ADHD Combined Type
    Pheno_Comb2$KeepDx1 <- ifelse(Pheno_Comb2$ASD_dummy == "1" | Pheno_Comb2$NoDX_dummy=="1", 1, 0)
    #Filter out ASD subjects with comorbid diagnosis
    Pheno_Comb2$KeepDx2 <- ifelse(Pheno_Comb2$ADHDCombinedType_dummy=="1" | Pheno_Comb2$ADHDHyperactiveImpulsiveType_dummy=="1" | Pheno_Comb2$ADHDInattentiveType_dummy=="1" | Pheno_Comb2$AdjustmentDisorders_dummy=="1" | Pheno_Comb2$BipolarIIDisorder_dummy=="1"| Pheno_Comb2$BorderlineIntellectualFunctioning_dummy=="1" | Pheno_Comb2$ConductDisorderChildhoodonsettype_dummy=="1" | Pheno_Comb2$DisruptiveMoodDysregulationDisorder_dummy=="1" | Pheno_Comb2$Encopresis_dummy=="1" | Pheno_Comb2$Enuresis_dummy=="1" | Pheno_Comb2$ExcoriationDisorder_dummy=="1" | Pheno_Comb2$GeneralizedAnxietyDisorder_dummy=="1" | Pheno_Comb2$IntellectualDisabilityModerate_dummy=="1" | Pheno_Comb2$IntellectualDisabilityMild_dummy=="1" | Pheno_Comb2$LanguageDisorder_dummy=="1" | Pheno_Comb2$MajorDepressiveDisorder_dummy=="1" | Pheno_Comb2$ObsessiveCompulsiveDisorder_dummy=="1" | Pheno_Comb2$OppositionalDefiantDisorder_dummy=="1" | Pheno_Comb2$OtherSpecifiedAnxietyDisorder_dummy=="1" | Pheno_Comb2$OtherSpecifiedAttentionDeficitHyperactivityDisorder_dummy=="1" | Pheno_Comb2$OtherSpecifiedDepressiveDisorder_dummy=="1" | Pheno_Comb2$OtherSpecifiedDisruptive_dummy=="1" | Pheno_Comb2$OtherSpecifiedTicDisorder_dummy=="1" | Pheno_Comb2$OtherSpecifiedTrauma_dummy=="1" | Pheno_Comb2$PanicDisorder_dummy=="1" | Pheno_Comb2$ParentChildRelationalProblem_dummy=="1" | Pheno_Comb2$PersistentDepressiveDisorder_dummy=="1" | Pheno_Comb2$PersistentMotorDisorder_dummy=="1" | Pheno_Comb2$PosttraumaticStressDisorder_dummy=="1" | Pheno_Comb2$SelectiveMutism_dummy=="1" | Pheno_Comb2$SeparationAnxiety_dummy=="1" | Pheno_Comb2$SocialAnxiety_dummy=="1" | Pheno_Comb2$SLDIMath_dummy=="1" | Pheno_Comb2$SLDReading_dummy=="1" | Pheno_Comb2$SLDWrittenExpression_dummy=="1" | Pheno_Comb2$SpecificPhobia_dummy=="1" | Pheno_Comb2$SpeechSoundDisorder_dummy=="1" | Pheno_Comb2$TourettesDisorder_dummy=="1" | Pheno_Comb2$UnspecifiedADHD_dummy=="1" | Pheno_Comb2$UnspecifiedAnxietyDisorder_dummy=="1" | Pheno_Comb2$UnspecifiedNeurodevelopmentalDisorder_dummy=="1", 0, 1)
    NoDx_Comb <- subset(Pheno_Comb2, KeepDx1!="0")
    NoDx_Comb <- subset(NoDx_Comb, KeepDx2!="0")
    
    #Create group var
    #Create var for any ASD diagnosis
    NoDx_Comb$ASDany <- ifelse(NoDx_Comb$DX_01 == "Autism Spectrum Disorder" | NoDx_Comb$DX_02 == "Autism Spectrum Disorder" | NoDx_Comb$DX_03 == "Autism Spectrum Disorder" | NoDx_Comb$DX_04 == "Autism Spectrum Disorder" | NoDx_Comb$DX_05 == "Autism Spectrum Disorder" | NoDx_Comb$DX_06 == "Autism Spectrum Disorder" | NoDx_Comb$DX_07 == "Autism Spectrum Disorder" | NoDx_Comb$DX_08 == "Autism Spectrum Disorder" | NoDx_Comb$DX_09 == "Autism Spectrum Disorder" | NoDx_Comb$DX_10 == "Autism Spectrum Disorder", "ASD", "NA")
    #Create var for any ADHD Combined dx
    NoDx_Comb$NoDx <- ifelse(is.na(NoDx_Comb$ASDany), "NoDX", " ")
    #Create overall group var
    NoDx_Comb$group <- (paste0(NoDx_Comb$ASDany, NoDx_Comb$NoDx))
    #Replace NAADHDComb with ADHDComb
    NoDx_Comb$group <- str_replace(NoDx_Comb$group, "NA", "NoDX")
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    NoDx_Comb$group <- trim(NoDx_Comb$group)
    NoDx_Comb$group <- as.factor(NoDx_Comb$group)
    
    #Create new mean-age and eTIV vars for this group  
    meanage = mean(NoDx_Comb$Age)
    NoDx_Comb$mc_age <- NoDx_Comb$Age - meanage
    NoDx_Comb$mc_age_squared <- NoDx_Comb$mc_age * NoDx_Comb$mc_age
    NoDx_Comb$mc_age_squared <- as.numeric(NoDx_Comb$mc_age_squared)
    
    #Make etiv and etiv-squared mean-centered
    etiv_mean = mean(NoDx_Comb$eTIV_cm)
    NoDx_Comb$eTIV_mc <- NoDx_Comb$eTIV_cm - etiv_mean
    NoDx_Comb$eTIV_squared_mc <- NoDx_Comb$eTIV_mc * NoDx_Comb$eTIV_mc
    
    #Frequency Stats
    table(NoDx_Comb$group)
    favstats(data=NoDx_Comb, Age~group)
    
    #ANCOVA: ADHDComb and ASD
    nodx_model = lm(EACSF_cm ~ group + mc_age + Site + Sex + mc_age*group + mc_age_squared + mc_age_squared*group + eTIV_mc + eTIV_squared_mc, data = NoDx_Comb)
    summary(nodx_model)
    
    #Plot for ADHDComb and ASD EA-CSF over Age
    require(ggplot2)
    #Plot
    Palette <- c("#00BCB8", "#0072B2")
    ggplot(NoDx_Comb, aes(x=Age, y=EACSF_cm, color=group))+
      labs(x = "Age (Years)", y = 'EA-CSF Volume ('~cm^3*')')+
      labs(fill = " ")+
      labs(color = " ")+
      scale_colour_manual(values=Palette, labels = c("Autism", "No Diagnosis"))+
      scale_fill_manual(values=Palette, labels = c("Autism", "No Diagnosis"))+
      geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
      geom_smooth(method=lm, aes(fill=group), se=TRUE) +
      scale_y_continuous(limits=c(20,150))+
      scale_x_continuous(limits=c(5,22), expand = c(0,0))+
      guides(color=guide_legend(override.aes=list(fill=NA)))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = c(.95, .2), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
            legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())
    
    ggsave(filename = paste("NoDx_220804.svg"), width = 3.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
    
    #FIGURE: Age violinplots by group
    Palette <- c("#00BCB8", "#0072B2", "#E69F00")
    ggplot(NoDx_Comb, aes(x = group, y = Age, fill = group)) +
      labs(x = "", y = 'Age')+
      geom_violin() +
      geom_boxplot(width=0.1, color="black", position=position_dodge(.91)) +
      scale_colour_manual(values=Palette, labels = c("Autism", "No Diagnosis"))+
      scale_fill_manual(values=Palette, labels = c("Autism", "No Diagnosis"))+
      scale_x_discrete(labels=c("Autism", "No Diagnosis")) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.title=element_blank(), legend.text=element_blank(), 
            legend.background = element_blank(), axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())+
      theme(legend.position = "none")
    
    #save the file
    ggsave(filename = paste("NoDx_AgeViolin_220804.svg"), width = 4.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses/specificDx_plots", dpi = 300)
    
      

#------------------------------SLEEP ANALYSIS--------------------------------

## Sleep: Is there a relationship between sleep disturbances and EA-CSF vol generally?
Pheno_Comb2$SDS_Total_Raw <- as.numeric(Pheno_Comb2$SDS_Total_Raw)
cor.test(Pheno_Comb2$EACSF_cm, Pheno_Comb2$SDS_Total_Raw)


#Freq. stats on participants with SDS scores
Pheno_CombS <- subset(Pheno_Comb2, SDS_Total_Raw!="NA")
favstats(Pheno_CombS$Age)
table(Pheno_CombS$group)
table(Pheno_CombS$Sex)

#Create new mean-age and eTIV vars for this group  
meanage = mean(Pheno_CombS$Age)
Pheno_CombS$mc_age <- Pheno_CombS$Age - meanage
Pheno_CombS$mc_age_squared <- Pheno_CombS$mc_age * Pheno_CombS$mc_age
Pheno_CombS$mc_age_squared <- as.numeric(Pheno_CombS$mc_age_squared)

#Make etiv and etiv-squared mean-centered
etiv_mean = mean(Pheno_CombS$eTIV_cm)
Pheno_CombS$eTIV_mc <- Pheno_CombS$eTIV_cm - etiv_mean
Pheno_CombS$eTIV_squared_mc <- Pheno_CombS$eTIV_mc * Pheno_CombS$eTIV_mc


#Assumptions
#Assess normality of residuals - Shapiro-Wilk test
#1. Fit the model, the covariate goes first
modelS <- lm(EACSF_cm ~ SDS_Total_Raw + mc_age + group + Site + Sex + eTIV_mc + eTIV_squared_mc + mc_age_squared, data = Pheno_Comb2)
#2. Run model metrics
model.metrics <- augment(modelS)
head(model.metrics, 3) #this line just displays the metrics
#3. Run the test
shapiro_test(model.metrics$.resid)
#If p-value was significant, reject null hypothesis that residuals are normal


#Multiple regression to take covariates into consideration
Sleepmodel = lm(EACSF_cm ~ SDS_Total_Raw + mc_age + group + Site + Sex + eTIV_mc + eTIV_squared_mc + mc_age_squared, data = Pheno_Comb2)
summary(Sleepmodel)


#FIGURE 4: SDS and EACSF Using LINEAR MODEL
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")
Pheno_Comb2$SDS_Total_Raw <- as.numeric(Pheno_Comb2$SDS_Total_Raw)

ggplot(Pheno_Comb2, aes(x=EACSF_cm, y=SDS_Total_Raw, color=group))+
  labs(x = 'EA-CSF Volume ('~cm^3*')', y = "Sleep Disturbance Scale Raw Total Score")+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=lm, aes(fill=group), se=TRUE) +
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
  scale_y_continuous(limits=c(20,85))+
  scale_x_continuous(limits=c(35,165), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.major = element_line(color = "grey"))+
  theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(0.8, 0.9), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#save the file
ggsave(filename = paste("Figure4_220707.svg"), width = 3.4, height = 3.4,
       path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/manuscript", dpi = 300)


#------------------------------SHEN 2018 SLEEP ANALYSIS---------------------------
##Follow analysis in Shen et al. 2018 on relationship between sleep and EA-CSF vol
Pheno_Comb2$SDS_Total_Raw <- as.numeric(Pheno_Comb2$SDS_Total_Raw)

#1. Group diff in sleep problems, including subgroup of ASD?
  #Create subgroups of ASD (<1.5 SD EACSF volume and >1.5SD EACSF volume)
    #Grab SD for ASD
    sd(data=Pheno_CombS, EACSF_cm~group)
    SD = 15.88526*1.5
    mean(data=Pheno_CombS, EACSF_cm~group)
    M = 68.72992
    ASD_SD = SD + M
    #Create subgroup
    Pheno_CombS$ASD1.5 <- ifelse(Pheno_CombS$group=="ASD" & Pheno_CombS$EACSF_cm >= ASD_SD, "ASD1.5", 0)
    table(Pheno_CombS$ASD1.5)
      #Only 6 ASD1.5 participants = not enough for subgroup comparison
    
  #Group (ASD, NoDx, OtherDx) comparison on sleep scores
    sds_model = lm(SDS_Total_Raw ~ group + mc_age + Site + Sex, data = Pheno_CombS)
    #Run ANOVA on model
    Anova(sds_model, type="III")  
    
    #Eta-squared (for effect sizes)
    etaSquared(sds_model, type=3)
    
    #Parse out group differences if omnibus had sig. group effect
    
    #Pairwise (Bonferroni) corrections using emmeans
    ##emmeans = estimated marginal means 
    ## = The estimated marginal means output gives the adjusted means (controlling for the covariates) for each group.
    #1. Get marginal means
    #marginal = emmeans(lmodel, ~ group + Age + Site + Sex + Age*group + Age_squared + Age_squared*group)
    marginal = emmeans(sds_model, ~ group*mc_age)
    #Display marginal means table
    marginal
    #2. Use Bonferroni corrections on marginal means, run pairwise corrections
    pairs(marginal, adjust="bon") 
    

  #Raincloud plot of sleep disturbances scores by group
    Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")
    ggplot(Pheno_CombS, aes(x = group, y = SDS_Total_Raw, fill=group)) + 
      ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        .width = 0, 
        justification = -.3, 
        point_colour = "NA") + 
      geom_boxplot(
        width = .25, 
        outlier.shape = NA
      ) +
      geom_point(
        size = 1.3,
        alpha = .3,
        position = position_jitter(
          seed = 1, width = .1
        )
      ) + 
      coord_cartesian((xlim = c(1.2, NA)), (ylim=c(25, 90)), clip = "off")+
      labs(y="Age (Years)", x="Diagnostic Group")+
      scale_colour_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
      scale_fill_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+  theme(axis.text=element_text(size = 9), axis.title = element_text(size = 9))+
      scale_x_discrete(labels=c("Autism", "No Diagnosis", "Other Diagnosis")) +
      theme(legend.position = "none")+
      theme(panel.background = element_blank())+
      theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"))
    
    #save the file
    ggsave(filename = paste("SDSxGroup_220715.svg"), width = 3.4, height = 3.4,
           path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/manuscript", dpi = 300)
    
    

#2. Relationship between EACSF to BrainVol ratio and sleep problems
    
    #create var with ratio between EACSF and BrainVol
    Pheno_CombS$EACSFratio <- Pheno_CombS$EACSF_cm / Pheno_CombS$BrainSegVolcm3
    favstats(Pheno_CombS$EACSFratio)

    
    #regression model for ratio ~ SDS scores by group
    sleepratio_model = lm(EACSFratio ~ SDS_Total_Raw + mc_age  + eTIV_mc + eTIV_squared_mc + group + Site + Sex + mc_age_squared, data = Pheno_CombS)
    summary(sleepratio_model)
    
    
    #Plot ratio by SDS score
    Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")
    Pheno_CombS$SDS_Total_Raw <- as.numeric(Pheno_CombS$SDS_Total_Raw)
    
    ggplot(Pheno_CombS, aes(x=EACSFratio, y=SDS_Total_Raw, color=group))+
      labs(x = 'Ratio of EA-CSF Volume to Brain Segmentation Volume', y = "Sleep Disturbance Scale Raw Total Score")+
      geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
      geom_smooth(method=lm, aes(fill=group), se=TRUE) +
      labs(fill = " ")+
      labs(color = " ")+
      scale_colour_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
      scale_fill_manual(values=Palette, labels = c("Autism", "No Diagnosis", "Other Diagnosis"))+
      scale_y_continuous(limits=c(20,85))+
      scale_x_continuous(limits=c(0,.15), expand = c(0,0))+
      theme_bw()+
      theme(panel.grid.major = element_line(color = "grey"))+
      theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
            axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = c(0.8, 0.9), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
            legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
            axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())
    
    

#---------------------------------MEAN CT ANALYSIS-----------------------------

#Mean CT analysis--what is the relationship between mean CT and EA-CSF?
#Assumptions testing
#Homogeneity of Regression Slopes
##Tests that there is no sig. relationship between group and age
Pheno_Comb2 %>% anova_test(EACSF_cm ~ group*totalmean_CT)
#group:Age term was not sig. so regression slopes look good!

#Assess normality of residuals - Shapiro-Wilk test
#1. Fit the model, the covariate goes first
model <- lm(EACSF_cm ~ totalmean_CT + mc_age + Site + Sex + eTIV_mc + eTIV_squared_mc + mc_age*group + mc_age_squared + mc_age_squared*group + group, data = Pheno_Comb2)
#2. Run model metrics
model.metrics <- augment(model)
head(model.metrics, 3) #this line just displays the metrics
#3. Run the test
shapiro_test(model.metrics$.resid)
#If p-value was significant, reject null hypothesis that residuals are normal

#Homogeneity of variance - Levene's test
model.metrics %>% levene_test(.resid ~ group)
#If p-value was significant, reject null hypothesis of homogeneity of variance for all groups



#Multiple regression to take covariates into consideration
CTmodel = lm(EACSF_cm~totalmean_CT + group + mc_age + Site + Sex + eTIV_mc + eTIV_squared_mc + mc_age_squared + mc_age*group + mc_age_squared*group, data = Pheno_Comb2)
summary(CTmodel)


# Figure 5A of mean CT over age
#Using LM 
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")
ggplot(Pheno_Comb2, aes(x=Age, y=totalmean_CT, color=group))+
  labs(x = 'Age (Years)', y = 'Total Mean Cortical Thickness (mm)')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=lm, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(2.25,3))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none", legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())
#save the file
ggsave(filename = paste("Figure5A_220707.svg"), width = 3.4, height = 3.4,
       path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/manuscript", dpi = 300)


# Figure 5B of mean CT and EA-CSF
#Using LM
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")
ggplot(Pheno_Comb2, aes(x=totalmean_CT, y=EACSF_cm, color=group))+
  labs(x = 'Total Mean Cortical Thickness (mm)', y = 'EA-CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=lm, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(20,150))+
  scale_x_continuous(limits=c(2.25,3.03), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=9), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none", legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())
#save the file
ggsave(filename = paste("Figure5B_220707.svg"), width = 3.7, height = 3.4,
       path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/manuscript", dpi = 300)



#-------------------------------OTHER ANALYSES-------------------------------------       

#Why is age_squared but not age a significant effect?

#Start by looking at just age -- it is a significant predictor by itself
Agemodel= lm(EACSF_cm ~ Age, data=Pheno_Comb2)
Anova(Agemodel)
#Systematically add in predictors to see why age wasn't significant in the overall model
#Predictors:group + Age + Site + Sex + Age*group + Age_squared + Age_squared*group
Agemodel2= lm(EACSF_cm ~ group + Age, data=Pheno_Comb2)
Anova(Agemodel2)
Agemodel3= lm(EACSF_cm ~ group + Age + Site, data=Pheno_Comb2)
Anova(Agemodel3)
Agemodel4= lm(EACSF_cm ~ group + Age + Site + Sex, data=Pheno_Comb2)
Anova(Agemodel4)
Agemodel5= lm(EACSF_cm ~ group + Age + Site + Sex + Age_squared, data=Pheno_Comb2)
Anova(Agemodel5)

#Just age and age_squared
Agemodel1= lm(EACSF_cm ~ Age + Age_squared, data=Pheno_Comb2)
Anova(Agemodel1)



#Analysis with total brain volume covariate -- since sex is a significant predictor, could this be due to differences in brain size?
#Run linear model first (library: car)
Pheno_CombTBV <- subset(Pheno_Comb2, BrainSegVolcm3!="NA")
TBVmodel = lm(EACSF_cm ~ group + mc_age + Site + Sex + BrainSegVolcm3 + BrainSegVolcm3_squared + mc_age*group + mc_age_squared + mc_age_squared*group, data = Pheno_CombTBV)
#Run ANOVA on model
Anova(TBVmodel)


#Analysis with estimated total intracranial volume (eTIV) covariate -- since sex is a significant predictor, could this be due to differences in brain size?
#mean-center eTIV_cm
favstats(Pheno_Comb2$eTIV_cm)
etiv_mean= 1490.982
Pheno_Comb2$eTIV_mc <- etiv_mean - Pheno_Comb2$eTIV_cm
Pheno_Comb2$eTIV_squared_mc <- Pheno_Comb2$eTIV_mc * Pheno_Comb2$eTIV_mc

#Run linear model first with eTIV
Pheno_CombeTIV <- subset(Pheno_Comb2, eTIV_cm!="NA")
eTIVmodel = lm(EACSF_cm ~ group + mc_age + Site + Sex + eTIV_mc + mc_age*group + mc_age_squared + mc_age_squared*group, data = Pheno_CombeTIV)
#Run ANOVA on model
Anova(eTIVmodel)  

#And model including eTIV_cm_squared covariate (since in the graph of eTIV by EA-CSF, it appears to have a nonlinear growth pattern)
eTIVmodel2 = lm(EACSF_cm ~ group + mc_age + Site + Sex + eTIV_mc + eTIV_squared_mc + mc_age*group + mc_age_squared + mc_age_squared*group, data = Pheno_CombeTIV)
#Run ANOVA on model
Anova(eTIVmodel2)  




#QC ICC analysis
#Compare ICC for MID02 (EA-CSF segmentation) ratings
QC$Notes <- NULL
QC$EID <- NULL
QC$X <- NULL
QC$T1w_QC_Anne <- NULL
QC$T1w_QC_Christopher <- NULL
QC$T1w_QC_Maddy <- NULL
QC$T1w_QC_Matt <- NULL
QC$Avg_Rating <- NULL

#ICC overall is 0.91
ICC(QC)

#Correlations
#Maddy versus Christopher
QC_Maddy <- subset(QC, EACSF_QC_Maddy!="NA")
QC_Maddy <- subset(QC_Maddy, EACSF_QC_Christopher!="NA")
cor(QC_Maddy$EACSF_QC_Maddy, QC_Maddy$EACSF_QC_Christopher)
# 0.975
#Plot
plot(jitter(QC_Maddy$EACSF_QC_Maddy), jitter(QC_Maddy$EACSF_QC_Christopher), pch = 16, col = 'steelblue')
abline(a = 0, b= 1, col="orange")   
#List of SUBJIDS with discrepant ratings (N=15)
QC_Maddy$discrepant <- ifelse(QC_Maddy$EACSF_QC_Maddy != QC_Maddy$EACSF_QC_Christopher, 1, 0)
table(QC_Maddy$discrepant)
QC_Maddy_mc <- subset(QC_Maddy, discrepant=="1")
list_mc <- table(QC_Maddy_mc$EID)



#Maddy versus Matt
QC_Maddy <- subset(QC, EACSF_QC_Maddy!="NA")
QC_Maddy <- subset(QC_Maddy, EACSF_QC_Matt!="NA")
cor(QC_Maddy$EACSF_QC_Maddy, QC_Maddy$EACSF_QC_Matt)
# 0.970
#Plot
plot(jitter(QC_Maddy$EACSF_QC_Maddy), jitter(QC_Maddy$EACSF_QC_Matt), pch = 16, col = 'steelblue')
abline(a = 0, b= 1, col="orange")   
#List of SUBJIDS with discrepant ratings (N= 19)
QC_Maddy$discrepant <- ifelse(QC_Maddy$EACSF_QC_Maddy != QC_Maddy$EACSF_QC_Matt, 1, 0)
table(QC_Maddy$discrepant)
QC_Maddy_mm <- subset(QC_Maddy, discrepant=="1")
list_mm <- table(QC_Maddy_mm$EID)



#Matt versus Anne 
QC_Matt <- subset(QC, EACSF_QC_Matt!="NA")
QC_Matt <- subset(QC_Matt, EACSF_QC_Anne!="NA")
cor(QC_Matt$EACSF_QC_Matt, QC_Matt$EACSF_QC_Anne)
# 0.754
#Plot
plot(jitter(QC_Matt$EACSF_QC_Matt), jitter(QC_Matt$EACSF_QC_Anne), pch = 16, col = 'steelblue')
abline(a = 0, b= 1, col="orange")  
#List of SUBJIDS with discrepant ratings (N=108)
QC_Matt$discrepant <- ifelse(QC_Matt$EACSF_QC_Matt != QC_Matt$EACSF_QC_Anne, 1, 0)
table(QC_Matt$discrepant)
QC_Matt_ma <- subset(QC_Matt, discrepant=="1")
list_ma <- table(QC_Matt_ma$EID)     



#Anne versus Christopher
QC_Anne <- subset(QC, EACSF_QC_Anne!="NA")
QC_Anne <- subset(QC_Anne, EACSF_QC_Christopher!="NA")
cor(QC_Anne$EACSF_QC_Anne, QC_Anne$EACSF_QC_Christopher)
# 0.635
#Plot
plot(jitter(QC_Anne$EACSF_QC_Anne), jitter(QC_Anne$EACSF_QC_Christopher), pch = 16, col = 'steelblue')
abline(a = 0, b= 1, col="orange")  
#List of SUBJIDS with discrepant ratings (N=181)    
QC_Anne$discrepant <- ifelse(QC_Anne$EACSF_QC_Anne != QC_Anne$EACSF_QC_Christopher, 1, 0)
table(QC_Anne$discrepant)
QC_Anne_ac <- subset(QC_Anne, discrepant=="1")
list_ac <- table(QC_Anne_ac$EID)     




##QC ANALYSIS II: Discrepancies Resolved???
QC2 <- subset(QC2, select = c("EID","Avg_Rating","EACSF_QC_Maddy","EACSF_QC_Matt","EACSF_QC_Christopher","EACSF_QC_Anne"))

QC_Maddy <- subset(QC2, EACSF_QC_Maddy!="NA")
QC_Maddy <- subset(QC_Maddy, EACSF_QC_Christopher!="NA")
cor(QC_Maddy$EACSF_QC_Maddy, QC_Maddy$EACSF_QC_Christopher)
QC_Maddy$discrepant <- ifelse(QC_Maddy$EACSF_QC_Maddy != QC_Maddy$EACSF_QC_Matt, 1, 0)
table(QC_Maddy$discrepant)
QC_Maddy_mc <- subset(QC_Maddy, discrepant=="1")
list_mc <- table(QC_Maddy_mc$EID)
#Corr = 1
#0 discrepant still

QC_Maddy <- subset(QC2, EACSF_QC_Maddy!="NA")
QC_Maddy <- subset(QC_Maddy, EACSF_QC_Matt!="NA")
cor(QC_Maddy$EACSF_QC_Maddy, QC_Maddy$EACSF_QC_Matt)
QC_Maddy$discrepant <- ifelse(QC_Maddy$EACSF_QC_Maddy != QC_Maddy$EACSF_QC_Matt, 1, 0)
table(QC_Maddy$discrepant)
QC_Maddy_mm <- subset(QC_Maddy, discrepant=="1")
list_mm <- table(QC_Maddy_mm$EID)
#Corr = 0.99
# 3 discrepant still

QC_Matt <- subset(QC2, EACSF_QC_Matt!="NA")
QC_Matt <- subset(QC_Matt, EACSF_QC_Anne!="NA")
cor(QC_Matt$EACSF_QC_Matt, QC_Matt$EACSF_QC_Anne)
QC_Matt$discrepant <- ifelse(QC_Matt$EACSF_QC_Matt != QC_Matt$EACSF_QC_Anne, 1, 0)
table(QC_Matt$discrepant)
QC_Matt_ma <- subset(QC_Matt, discrepant=="1")
list_ma <- table(QC_Matt_ma$EID)  
#Corr = 0.95
# 16 discrepant still

QC_Anne <- subset(QC2, EACSF_QC_Anne!="NA")
QC_Anne <- subset(QC_Anne, EACSF_QC_Christopher!="NA")
cor(QC_Anne$EACSF_QC_Anne, QC_Anne$EACSF_QC_Christopher)
QC_Anne$discrepant <- ifelse(QC_Anne$EACSF_QC_Anne != QC_Anne$EACSF_QC_Christopher, 1, 0)
table(QC_Anne$discrepant)
QC_Anne_ac <- subset(QC_Anne, discrepant=="1")
list_ac <- table(QC_Anne_ac$EID)  
#Corr = 0.99
# 3 discrepant still



#Are there failed raw T1w QC ratings in the dataset? (after excluding the failed EACSF ratings)    
Pheno_Comb2$T1w_QC_Maddy <- trim(Pheno_Comb2$T1w_QC_Maddy)
Pheno_Comb2$T1w_QC_Matt <- trim(Pheno_Comb2$T1w_QC_Matt)
Pheno_Comb2$T1w_QC_Anne <- trim(Pheno_Comb2$T1w_QC_Anne)
Pheno_Comb2$T1w_QC_Christopher <- trim(Pheno_Comb2$T1w_QC_Christopher)

table(Pheno_Comb2$T1w_QC_Anne) # Check: 72, Fail: 15
table(Pheno_Comb2$T1w_QC_Matt) # Check: 102, Fail: 64
table(Pheno_Comb2$T1w_QC_Christopher) # Check 108, Fail: 41
table(Pheno_Comb2$T1w_QC_Maddy) # Check: 130, Fail: 66

#Create avg T1w QC variable -- see overall how many failed
#Convert to numerical values (run trim lines first!!)
Pheno_Comb2$T1w_QC_Anne <- gsub("^.{0,1}", "", Pheno_Comb2$T1w_QC_Anne)
Pheno_Comb2$T1w_QC_Matt <- gsub("^.{0,1}", "", Pheno_Comb2$T1w_QC_Matt)
Pheno_Comb2$T1w_QC_Christopher <- gsub("^.{0,1}", "", Pheno_Comb2$T1w_QC_Christopher)
Pheno_Comb2$T1w_QC_Maddy <- gsub("^.{0,1}", "", Pheno_Comb2$T1w_QC_Maddy)

#Avg across rows
Pheno_Comb2$T1w_QC_Matt[ Pheno_Comb2$T1w_QC_Matt == ""] <- NA
Pheno_Comb2$T1w_QC_Anne[ Pheno_Comb2$T1w_QC_Anne == ""] <- NA
Pheno_Comb2$T1w_QC_Christopher[ Pheno_Comb2$T1w_QC_Christopher == ""] <- NA
Pheno_Comb2$T1w_QC_Maddy[ Pheno_Comb2$T1w_QC_Maddy == ""] <- NA

qc_t1 <- data.frame(Pheno_Comb2$T1w_QC_Anne, Pheno_Comb2$T1w_QC_Matt, Pheno_Comb2$T1w_QC_Christopher, Pheno_Comb2$T1w_QC_Maddy)
qc_t1$Pheno_Comb2.T1w_QC_Anne <-as.numeric(qc_t1$Pheno_Comb2.T1w_QC_Anne)
qc_t1$Pheno_Comb2.T1w_QC_Matt <-as.numeric(qc_t1$Pheno_Comb2.T1w_QC_Matt)
qc_t1$Pheno_Comb2.T1w_QC_Christopher <-as.numeric(qc_t1$Pheno_Comb2.T1w_QC_Christopher)
qc_t1$Pheno_Comb2.T1w_QC_Maddy <-as.numeric(qc_t1$Pheno_Comb2.T1w_QC_Maddy)

Pheno_Comb2$T1w_Avg_Rating <- rowMeans(qc_t1, na.rm=TRUE)

# See how many had a failed average rating 
table(Pheno_Comb2$T1w_Avg_Rating)
#Check: 208, Fail: 77


#subset to failed T1w rating in order to get SUBJIDS to check out the situation
Pheno_CombFail <- subset(Pheno_Comb2, T1w_Avg_Rating == "3")
table(Pheno_CombFail$EID, Pheno_CombFail$T1)      

# Is there a correlation between a failed T1w and failed MID02 QC ratings??
#Get all ratings by running everything in the setup *except* the line where QC scores are removed
Pheno_Comb2$Fail_T1 <- ifelse(Pheno_Comb2$T1w_Avg_Rating == "3", 1, 0)
Pheno_Comb2$Fail_MID02 <- ifelse(Pheno_Comb2$Avg_Rating >= "2", 1, 0)
cor(Pheno_Comb2$Fail_T1, Pheno_Comb2$Fail_MID02) #Corr = 0.345


cor(Pheno_Comb2$Avg_Rating, Pheno_Comb2$T1w_Avg_Rating) #cor between T1 and MID02 ratings = 0.452
plot(Pheno_Comb2$Avg_Rating, Pheno_Comb2$T1w_Avg_Rating, main="Comparison of T1 and MID02 QC Ratings",
     xlab="MID02 (EACSF) Avg Rating ", ylab="T1w Avg Rating", pch=19) 
abline(lm(Pheno_Comb2$T1w_Avg_Rating~Pheno_Comb2$Avg_Rating), col="red")
abline(a = 0, b= 1, col="blue")         


# Re-centering age 5-21 versus mean
#Run linear model first-- age is mean-centered due to collinearity problems
mcAgemodel = lm(EACSF_cm ~ group + mc_age + Site + Sex + mc_age*group + mc_age_squared + mc_age_squared*group, data = Pheno_Comb2)
#Run ANOVA on model
Anova(mcAgemodel)

#5
Pheno_Comb2$age5 <- 5 - Pheno_Comb2$Age
Pheno_Comb2$age_squared5 <- Pheno_Comb2$age5 * Pheno_Comb2$age5
agemodel5 = lm(EACSF_cm ~ group + age5 + Site + Sex + age5*group + age_squared5 + age_squared5*group, data = Pheno_Comb2)
Anova(agemodel5)

#6
Pheno_Comb2$age6 <- 6 - Pheno_Comb2$Age
Pheno_Comb2$age_squared6 <- Pheno_Comb2$age6 * Pheno_Comb2$age6
agemodel6 = lm(EACSF_cm ~ group + age6 + Site + Sex + age6*group + age_squared6 + age_squared6*group, data = Pheno_Comb2)
Anova(agemodel6)

#7
Pheno_Comb2$age7 <- 7 - Pheno_Comb2$Age
Pheno_Comb2$age_squared7 <- Pheno_Comb2$age7 * Pheno_Comb2$age7
agemodel7 = lm(EACSF_cm ~ group + age7 + Site + Sex + age7*group + age_squared7 + age_squared7*group, data = Pheno_Comb2)
Anova(agemodel7)

#8
Pheno_Comb2$age8 <- 8 - Pheno_Comb2$Age
Pheno_Comb2$age_squared8 <- Pheno_Comb2$age8 * Pheno_Comb2$age8
agemodel8 = lm(EACSF_cm ~ group + age8 + Site + Sex + age8*group + age_squared8 + age_squared8*group, data = Pheno_Comb2)
Anova(agemodel8)

#9
Pheno_Comb2$age9 <- 9 - Pheno_Comb2$Age
Pheno_Comb2$age_squared9 <- Pheno_Comb2$age9 * Pheno_Comb2$age9
agemodel9 = lm(EACSF_cm ~ group + age9 + Site + Sex + age9*group + age_squared9 + age_squared9*group, data = Pheno_Comb2)
Anova(agemodel9)

#10
Pheno_Comb2$age10 <- 10 - Pheno_Comb2$Age
Pheno_Comb2$age_squared10 <- Pheno_Comb2$age10 * Pheno_Comb2$age10
agemodel10 = lm(EACSF_cm ~ group + age10 + Site + Sex + age10*group + age_squared10 + age_squared10*group, data = Pheno_Comb2)
Anova(agemodel10)

#11 (groupxage)
Pheno_Comb2$age11 <- 11 - Pheno_Comb2$Age
Pheno_Comb2$age_squared11 <- Pheno_Comb2$age11 * Pheno_Comb2$age11
agemodel11 = lm(EACSF_cm ~ group + age11 + Site + Sex + age11*group + age_squared11 + age_squared11*group, data = Pheno_Comb2)
Anova(agemodel11)

#12
Pheno_Comb2$age12 <- 12 - Pheno_Comb2$Age
Pheno_Comb2$age_squared12 <- Pheno_Comb2$age12 * Pheno_Comb2$age12
agemodel12 = lm(EACSF_cm ~ group + age12 + Site + Sex + age12*group + age_squared12 + age_squared12*group, data = Pheno_Comb2)
Anova(agemodel12)

#13
Pheno_Comb2$age13 <- 13 - Pheno_Comb2$Age
Pheno_Comb2$age_squared13 <- Pheno_Comb2$age13 * Pheno_Comb2$age13
agemodel13 = lm(EACSF_cm ~ group + age13 + Site + Sex + age13*group + age_squared13 + age_squared13*group, data = Pheno_Comb2)
Anova(agemodel13)

#14
Pheno_Comb2$age14 <- 14 - Pheno_Comb2$Age
Pheno_Comb2$age_squared14 <- Pheno_Comb2$age14 * Pheno_Comb2$age14
agemodel14 = lm(EACSF_cm ~ group + age14 + Site + Sex + age14*group + age_squared14 + age_squared14*group, data = Pheno_Comb2)
Anova(agemodel14)

#15
Pheno_Comb2$age15 <- 15 - Pheno_Comb2$Age
Pheno_Comb2$age_squared15 <- Pheno_Comb2$age15 * Pheno_Comb2$age15
agemodel15 = lm(EACSF_cm ~ group + age15 + Site + Sex + age15*group + age_squared15 + age_squared15*group, data = Pheno_Comb2)
Anova(agemodel15)

#16
Pheno_Comb2$age16 <- 16 - Pheno_Comb2$Age
Pheno_Comb2$age_squared16 <- Pheno_Comb2$age16 * Pheno_Comb2$age16
agemodel16 = lm(EACSF_cm ~ group + age16 + Site + Sex + age16*group + age_squared16 + age_squared16*group, data = Pheno_Comb2)
Anova(agemodel16)

#17
Pheno_Comb2$age17 <- 17 - Pheno_Comb2$Age
Pheno_Comb2$age_squared17 <- Pheno_Comb2$age17 * Pheno_Comb2$age17
agemodel17 = lm(EACSF_cm ~ group + age17 + Site + Sex + age17*group + age_squared17 + age_squared17*group, data = Pheno_Comb2)
Anova(agemodel17)

#18
Pheno_Comb2$age18 <- 18 - Pheno_Comb2$Age
Pheno_Comb2$age_squared18 <- Pheno_Comb2$age18 * Pheno_Comb2$age18
agemodel18 = lm(EACSF_cm ~ group + age18 + Site + Sex + age18*group + age_squared18 + age_squared18*group, data = Pheno_Comb2)
Anova(agemodel18)

#19
Pheno_Comb2$age19 <- 19 - Pheno_Comb2$Age
Pheno_Comb2$age_squared19 <- Pheno_Comb2$age19 * Pheno_Comb2$age19
agemodel19 = lm(EACSF_cm ~ group + age19 + Site + Sex + age19*group + age_squared19 + age_squared19*group, data = Pheno_Comb2)
Anova(agemodel19)

#20
Pheno_Comb2$age20 <- 20 - Pheno_Comb2$Age
Pheno_Comb2$age_squared20 <- Pheno_Comb2$age20 * Pheno_Comb2$age20
agemodel20 = lm(EACSF_cm ~ group + age20 + Site + Sex + age20*group + age_squared20 + age_squared20*group, data = Pheno_Comb2)
Anova(agemodel20)

#21
Pheno_Comb2$age21 <- 21 - Pheno_Comb2$Age
Pheno_Comb2$age_squared21 <- Pheno_Comb2$age21 * Pheno_Comb2$age21
agemodel21 = lm(EACSF_cm ~ group + age21 + Site + Sex + age21*group + age_squared21 + age_squared21*group, data = Pheno_Comb2)
Anova(agemodel21)




#How do QCistern and MID02 EACSF volumes compare??
MID02_EACSF <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/final_results/EaCSF_results_1.7.7_T1only_release1-7_1611851892.csv")
QCistern_EACSF <- read.csv("C:/Users/maddy/Box/Autism_CSF/HBN_Project/HBN_data/final_results/EaCSF_results_1.7.7_T1only_release1-7_QCistern1615997661.csv")
names(MID02_EACSF)[2] <- "MID02_EACSF"
names(QCistern_EACSF)[2] <- "QCistern_EACSF"
MID02_EACSF$EID <- gsub("^.{0,4}", "", MID02_EACSF$EID)
MID02_EACSF$EID <- gsub(".{1}$", "", MID02_EACSF$EID)
QCistern_EACSF$EID <- gsub("^.{0,4}", "", QCistern_EACSF$EID)
QCistern_EACSF$EID <- gsub(".{1}$", "", QCistern_EACSF$EID)
eacsfmerge <- merge(MID02_EACSF, QCistern_EACSF, by="EID", all=TRUE)
eacsfmerge$MID02_EACSF<- as.numeric(eacsfmerge$MID02_EACSF)
eacsfmerge$QCistern_EACSF<- as.numeric(eacsfmerge$QCistern_EACSF)
eacsfmerge <- subset(eacsfmerge, QCistern_EACSF!="NA")

#Correlation!
cor(eacsfmerge$MID02_EACSF, eacsfmerge$QCistern_EACSF) #0.992

#ICC= 0.99 (since EACSF volumes are coming from the same population (image))
eacsfmerge2 <- eacsfmerge
eacsfmerge2$EID <- NULL
ICC(eacsfmerge2)

#plot of QCistern by MID02 EACSF volumes
eacsfmerge$MID02_cm <- eacsfmerge$MID02_EACSF/1000
eacsfmerge$QCistern_cm <- eacsfmerge$QCistern_EACSF/1000

ggplot(eacsfmerge, aes(x=MID02_cm, y=QCistern_cm))+
  labs(x ='MID02 EACSF ('~cm^3*')', y ='QCistern EACSF ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  geom_abline(slope=1, intercept = 0, size=1, linetype= "dashed")+
  geom_point(colour="black", pch=21, cex=1)+
  geom_smooth(method=lm, se=TRUE) +
  scale_y_continuous(limits=c(0,170))+
  scale_x_continuous(limits=c(0,170), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=14), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.85, .2), legend.title=element_text(colour = "black", size = 14), legend.text=element_text(colour = "black", size = 18), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#On average, is there a difference? If so, is MID02 higher?
#reshape to long 
eacsfmerge_long <- reshape(eacsfmerge, direction='long', 
                           varying=c('MID02_cm', 'QCistern_cm'),
                           timevar='file',
                           times=c('MID02', 'QCistern'),
                           v.names=c('eacsf'),
                           idvar='EID')
t.test(eacsf~file, data=eacsfmerge_long, mu = 0, alternative = "greater", conf.level = 0.95)
#So yes, the MID02 mean is greater than the QCistern (69.57 versus 61.86)






#---------------------------------EXAMPLE FIGURES---------------------------------------------------------

#Plot marginal means 
plot(marginal)
plot(marginal, by = NULL, horizontal = FALSE, CIs= TRUE, xlab = "Estimated Marginal Mean", ylab = "Diagnostic Group", colors = "seagreen") + theme_bw()    


#Violin plot with marginal means
#install.packages("modelbased")
library(modelbased)
means_complex <- estimate_means(eTIVmodel2)


ggplot(Pheno_Comb2, aes(x = group, y = EACSF_cm, fill = group)) +
  labs(x = "Diagnostic Group", y = 'EA-CSF Volume ('~cm^3*')')+
  geom_violin() +
  geom_line(data = means_complex, aes(y = Mean, group = 1), size = 1) +
  geom_pointrange(data = means_complex,
                  aes(y = Mean, ymin = CI_low, ymax = CI_high), 
                  size = 1, color="white") +
  theme_bw()



#Example code for a histogram
Pheno_Comb2$Age_bins <- as.numeric(Pheno_Comb2$Age_bins)      
hist(Pheno_Comb2$Age_bins, col="skyblue", xlab="Age", main="Ages of Participants Included in Analysis")


#Example code for a boxplot of EACSF distribution by participant sex
boxplot(EACSF_cm~Sex,
        data=Pheno_Comb2,
        main="EACSF Volume by Participant Sex",
        xlab="Sex (0 = Male, 1 = Female)",
        ylab="EACSF (in cm3)",
        col="seagreen",
        border="green"
)



#Example code for a boxplot of EACSF distribution by group
boxplot(EACSF_cm~group,
        data=Pheno_Comb2,
        main="EACSF Volume by Group",
        xlab="Group",
        ylab="EACSF (in cm3)",
        col=Palette,
        border="black"
)

#Example code for a boxplot of EACSF distribution by group using GADgroup
boxplot(EACSF_cm~GADgroup,
        data=Pheno_CombG,
        main="EACSF Volume by Group",
        xlab="Group",
        ylab="EACSF (in cm3)",
        col=Palette,
        border="black"
)      


#Example code you can adapt to your purposes--this is sized for a poster and has high res.

#set palette, prep plot
require(ggplot2)

#Using all three groups
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")

#plot age by CSF (loess smoother)

ggplot(Pheno_Comb2, aes(x=Age, y=EACSF_cm, color=group))+
  labs(x = "Age (Years)", y = 'EA-CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(20,150))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=14), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.85, .2), legend.title=element_text(colour = "black", size = 14), legend.text=element_text(colour = "black", size = 14), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())


#Plot of TBV by EACSF
#Using all three groups
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")

#plot age by CSF (loess smoother)
Pheno_CombTBV <- subset(Pheno_Comb2, BrainSegVolcm3!="NA")

ggplot(Pheno_CombTBV, aes(x=EACSF_cm, y=BrainSegVolcm3, color=group))+
  labs(x ='EA-CSF Volume ('~cm^3*')', y ='Total Brain Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(700,1600))+
  scale_x_continuous(limits=c(25,150), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=14), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.85, .2), legend.title=element_text(colour = "black", size = 14), legend.text=element_text(colour = "black", size = 18), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())



#Plot of eTIV by EACSF
#Using all three groups
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")

#plot age by CSF (loess smoother)
Pheno_CombeTIV <- subset(Pheno_Comb2, eTIV_cm!="NA")

ggplot(Pheno_CombeTIV, aes(x=EACSF_cm, y=BrainSegVolcm3, color=group))+
  labs(x ='EA-CSF Volume ('~cm^3*')', y ='Estimated Total Intracranial Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(800,1600))+
  scale_x_continuous(limits=c(25,150), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.85, .2), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#Plot of eTIV by Age
#Using all three groups
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")

#plot age by CSF (loess smoother)
Pheno_CombeTIV <- subset(Pheno_Comb2, eTIV_cm!="NA")

ggplot(Pheno_CombeTIV, aes(x=Age, y=BrainSegVolcm3, color=group))+
  labs(x ='Age (years)', y ='Estimated Total Intracranial Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(800,1600))+
  scale_x_continuous(limits=c(5, 23), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(1.3, .12), legend.title=element_text(colour = "black", size = 9), legend.text=element_text(colour = "black", size = 9), 
        legend.background = element_rect(fill="white", size=0.4) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())



#save your figure to Box when you are done
#ggsave(filename = paste("QCistern_versus_MID02_210317.svg"), width = 13, height = 6,
#      path = "C:/Users/maddy/Box/Autism_CSF/HBN_Project/figures/other_analyses", dpi = 300)



#-----------------------------CSF COMPARTMENT FIGURES------------------------------------        
# A series of figures examining the relationship between various CSF compartments and EA-CSF volume

#Create variables in cm3
Pheno_Comb2$total_vent <- Pheno_Comb2$Left.Lateral.Ventricle + Pheno_Comb2$Right.Lateral.Ventricle + Pheno_Comb2$Left.Inf.Lat.Vent + Pheno_Comb2$Right.Inf.Lat.Vent + Pheno_Comb2$X3rd.Ventricle + Pheno_Comb2$X4th.Ventricle + Pheno_Comb2$X5th.Ventricle
Pheno_Comb2$total_vent_cm <- Pheno_Comb2$total_vent/1000  
Pheno_Comb2$Left.Lateral.Ventricle_cm <- Pheno_Comb2$Left.Lateral.Ventricle/1000
Pheno_Comb2$Right.Lateral.Ventricle_cm <- Pheno_Comb2$Right.Lateral.Ventricle/1000
Pheno_Comb2$Left.Inf.Lat.Vent_cm <- Pheno_Comb2$Left.Inf.Lat.Vent/1000
Pheno_Comb2$Right.Inf.Lat.Vent_cm <- Pheno_Comb2$Right.Inf.Lat.Vent/1000
Pheno_Comb2$X3rd.Ventricle_cm <- Pheno_Comb2$X3rd.Ventricle/1000
Pheno_Comb2$X4th.Ventricle_cm <- Pheno_Comb2$X4th.Ventricle/1000
Pheno_Comb2$X5th.Ventricle_cm <- Pheno_Comb2$X5th.Ventricle/1000

#Create the 8 plots        
Palette <- c("#00BCB8", "#0072B2", "#E69F00", "#0072B2")

#plot age by CSF (loess smoother)

#1. Total ventricular CSF by age
ggplot(Pheno_Comb2, aes(x=Age, y=total_vent_cm, color=group))+
  labs(x = "Age (Years)", y = 'Total Ventricular CSF ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(0,80))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.15, .9), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#2. Right Ventricle by Age
ggplot(Pheno_Comb2, aes(x=Age, y=Right.Lateral.Ventricle_cm, color=group))+
  labs(x = "Age (Years)", y = 'Right Lateral Ventricle CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(0, 30))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.15, .9), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#3. Left Ventricle by Age
ggplot(Pheno_Comb2, aes(x=Age, y=Left.Lateral.Ventricle_cm, color=group))+
  labs(x = "Age (Years)", y = 'Left Lateral Ventricle CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(0, 30))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.15, .9), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#4. Right Inferior Ventricle by Age
ggplot(Pheno_Comb2, aes(x=Age, y=Right.Inf.Lat.Vent_cm, color=group))+
  labs(x = "Age (Years)", y = 'Right Inf. Lateral Ventricle CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(0, 5))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.15, .9), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#5. Left Inferior Ventricle by Age
ggplot(Pheno_Comb2, aes(x=Age, y=Left.Inf.Lat.Vent_cm, color=group))+
  labs(x = "Age (Years)", y = 'Left Inf. Lateral Ventricle CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(0, 5))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.15, .9), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())


#6. 3rd Ventricle by Age
ggplot(Pheno_Comb2, aes(x=Age, y=X3rd.Ventricle_cm, color=group))+
  labs(x = "Age (Years)", y = 'Third Ventricle CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(0, 5))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.15, .9), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#7. 4th Ventricle by Age
ggplot(Pheno_Comb2, aes(x=Age, y=X4th.Ventricle_cm, color=group))+
  labs(x = "Age (Years)", y = 'Fourth Ventricle CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(0, 5))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.15, .9), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())

#8. 5th Ventricle by Age
ggplot(Pheno_Comb2, aes(x=Age, y=X5th.Ventricle_cm, color=group))+
  labs(x = "Age (Years)", y = 'Fifth Ventricle CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(0, 1))+
  scale_x_continuous(limits=c(5,22), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.15, .9), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())


#Create plot of EA-CSF by total ventricular CSF      

ggplot(Pheno_Comb2, aes(x=total_vent_cm, y=EACSF_cm, color=group))+
  labs(x = 'Total Ventricular CSF Volume ('~cm^3*')', y = 'EA-CSF Volume ('~cm^3*')')+
  labs(fill = " ")+
  labs(color = " ")+
  scale_colour_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  scale_fill_manual(values=Palette, labels = c("Autism", "Control", "OtherDx"))+
  geom_point(aes(fill=group), colour="black", pch=21, cex=1)+
  geom_smooth(method=loess, aes(fill=group), se=TRUE) +
  scale_y_continuous(limits=c(15,140))+
  scale_x_continuous(limits=c(5,60), expand = c(0,0))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title = element_text(colour = "black", size=12), axis.text.y =element_text(colour = "black", angle = 90, hjust = 0.6), 
        axis.text.x =element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.9, .18), legend.title=element_text(colour = "black", size = 12), legend.text=element_text(colour = "black", size = 12), 
        legend.background = element_rect(fill="white", size=0.5) , axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
        axis.ticks = element_line(colour = "black", size =1, linetype ="solid"), panel.border = element_blank(), panel.background = element_blank())







