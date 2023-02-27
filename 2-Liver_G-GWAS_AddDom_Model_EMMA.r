
##############################################################################################################
### ADDO_AddDom :
### 1. PLINK Input Format: (1) file.tped & file.tfam (2) file.phe (1st column name should be "id"; From 2nd column should start from covariates and the sex column should coded as 0=female and 1=male; The sex column is required!) (3) file.covs (1st column is phenotype name; 2nd column is corresponding covariates and all covariates should be separated by ",")
### 2. GenABEL Input Format: file.ABEL.dat (Just contain one GenABEL type variable named "dat") & file.covs (1st column is phenotype name; 2nd column is corresponding covariates and all covariates should be separated by ",")
### 3. Required Softwares: plink (v.1.90) & emmax-kin/gemma/gcta64 (just install needed one)
##############################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

# [1] Load library and function #

library(ADDO)

# [2] Specify directory and covariates types variable and focused phenotypes list (if want run all phenotypes, just set "PheList_Choose=F") # 

indir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/0_Data"
outdir = "/SAN/mottlab/heterosis/3_HSrats/8_sQTL_Thu/2_Liver_Genes"
Input_name = "HSrats_liver_Genes"
covariates_types = c("n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n","n");
names(covariates_types) = c("batch2","batch3","batch4","batch5","batch6","batch7","batch8","batch9","batch10","batch11","batch12","batch13","batch14","batch15","batch16","batch17","batch18","batch19","batch20","batch21","batch22","batch23","batch24","batch25","batch26","batch27","batch28","batch29","batch30","batch31","batch32","batch33","RIN","expr.pc1","pca1","pca2","pca3","pca4","pca5","pca6","pca7","pca8","pca9","pca10")

# [3] Use three functions one by one #

# ADDO_AddDom1_QC(indir=indir, outdir=outdir, Input_name=Input_name, Input_type="PLINK", Kinship_type="EMMA", PheList_Choose=F, PheList=PheList, Phe_ResDone = F, Phe_NormDone = F, Normal_method = "QUANTILE", covariates_sum=44, covariates_types=covariates_types, Phe_IndMinimum = 100, Phe_Extreme = 5, GT_maf = 0.05, GT_missing = 0.1, num_nodes=6)

ADDO_AddDom2_Pvalue(indir=indir, outdir=outdir, Input_name=Input_name, Kinship_type="EMMA", VarComponent_Method="EMMA_a", PheList_Choose=F, PheList=PheList, covariates_sum=44, Phe_IndMinimum=100, GT_IndMinimum=10, matrix_acceleration=T, num_nodes=10)

# ADDO_AddDom3_Plot(outdir=outdir, PheList_Choose=F, PheList=PheList, covariates_sum=44, RegionMan_chr_whole=F, RegionMan_chr_region = 2000000, chrs_sum=21, num_nodes=6)

# ADDO_AddDom4_IntePlot(outdir=outdir, PheList_Choose=F, PheList=PheList, covariates_sum=44, RegionMan_chr_whole=F, RegionMan_chr_region = 2000000, Plot_model="NvsAD", chrs_sum=21, num_nodes=6)


