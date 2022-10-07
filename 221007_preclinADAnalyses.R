################################################################################
# Statistical Analyses for Ferreiro et al., 'The gut microbiome as an early    #
# biomarker of preclinical Alzheimer disease'. (2022)                          #  
#                                                                              #
# Bioconductor version 3.14 (BiocManager 1.30.18), R 4.1.3 (2022-03-10)        #
#                                                                              #
# Contact: Aura Ferreiro (alferreiro@wustl.edu)                                #
################################################################################
#                                 10/07/2022                                   #
# IMPORTANT NOTE: THIS R SCRIPT IS BEING CLEANED UP FOR RELEASE WHILE THE      #
# PUBLICATION IS UNDER REVIEW, AND IS INCOMPLETE. FOR URGENT QUESTIONS ABOUT   #
# THE ANALYSIS, PLEASE CONTACT ALFERREIRO@WUSTL.EDU                            #
################################################################################


library(dplyr)    #v1.0.8   
library(tidyr)    #v1.2.0
library(phyloseq) #v1.38.0
library(ggplot2)  #v3.3.5
library(ape)      #v5.6.2
library(Maaslin2) #v1.10.0

##################### LOAD DATA ################################################
# meta: metadata for 164 participants (see Supplementary Data File 1 for Data Dictionary)
# metaphlan3.raw : taxonomic relative abundance table, unfiltered 
# humann3.raw : pathway relative abundance table, unfiltered
# m3.tree : MetaPhlAn3 phylogenetic tree (mpa_v30_CHOCOPhlAn_201901_species.tree.nwk)
# diet.intake : stool-matched intake of macro and essential nutrient categories (24 hour)

load('RData/ACS_AD/220322_STMRevision/Data/Github/221007_preclinADAnalyses.RData')



############### CREATE TAXONOMIC PHYLOSEQ OBJECT ###############################

## HOUSEKEEPING ##
rownames(meta) <- meta$Participant

#Convert appropriate variables to factor
vars <- c('Participant', 'Sex', 'Race', 'amyloid.positive.AF', 'APOE', 
          'APOE4.status', 'Cancer', 'Active.Depression', 'Hypercholesterolemia',
          'Hypertension', 'Diabetes', 'Liver.Disease', 'Alcohol.Abuse',
          'Tobacco.Use.PastOrPresent', 'Autoimmune.Disorder', 'Thyroid.Disease',
          'Cardiovascular.Disease')


meta[,vars] <- lapply(meta[,vars], factor)

