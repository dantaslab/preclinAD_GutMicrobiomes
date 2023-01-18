# AUTHOR INFORMATION ###########################################################
#                                                                              #
# Statistical Analyses for Ferreiro et al., 'The gut microbiome as an early    #
# biomarker of preclinical Alzheimer disease'. (2022)                          #  
#                                                                              #
# Bioconductor version 3.14 (BiocManager 1.30.18), R 4.1.3 (2022-03-10)        #
#                                                                              #
# Contact: Aura Ferreiro (alferreiro@wustl.edu)                                #
#                                                                              #   
#                                                                              #
# LOAD PACKAGES ################################################################

# If you are only interested in the Random Forest Analyses, proceed to that 
# section. There you will load workspace 230111_ps2filt.RData, and the packages
# specific to those analyses. Otherwise, continue! 

library(dplyr)        #v1.0.8   
library(tidyr)        #v1.2.0
library(phyloseq)     #v1.38.0
library(ggplot2)      #v3.3.5
library(ape)          #v5.6.2
library(Maaslin2)     #v1.10.0
library(tableone)     #v0.13.2
library(knitr)        #v1.39
library(RColorBrewer) #v1.1.3
library(forcats)      #v0.5.1
library(scales)       #v1.1.1
library(vegan)        #v2.5.7
library(ggExtra)      #v0.10.0
library(Hmisc)        #v4.7.0
library(corrplot)     #v0.92
library(corpcor)      #1.6.10
library(Boruta)       #v7.0.0
library(caret)        #v6.0.86
library(VIM)          #v6.1.0  
library(stringr)      #v1.4.0
library(egg)          #v0.4.5
library(purrr)        #v0.3.4
library(broom)        #v1.0.1





# LOAD DATA ####################################################################

# 230118_preclniicalAD_Ferreiro_etal_DATA.RData:
#
#  meta                  =  Metadata for 164 participants. See Supp. Data File 1 
#                           for Data Dictionary. Note: 'meta'contains some 
#                           variables that will be re-derived in this script 
#                           (e.g. PCoA coordinates) to illustrate analyses. 
#  metaphlan3.raw        =  Taxonomic relative abundance table, unfiltered. 
#  humann3.raw           =  Pathway relative abundance table, unfiltered.
#  m3.tree               =  MetaPhlAn3 phylogenetic tree: 
#                           (mpa_v30_CHOCOPhlAn_201901_species.tree.nwk)
#  genus_colors.df       =  Custom color palette for stacked taxonomic barplots.
#  diffspecies.all.df    =  Significant species from negative binomial 
#                           regression analyses (MaAslin2).
#  diffpath.all.df       =  Significant pathways from negative binomial
#                           regression analyses (MaAslin2).
#  diet.intake           =  Stool-matched 24hr intake of major nutrients 
#                           (Percent daily values) for each participant. 
#  total.cal             =  Stool-matched 24hr caloric intake by source.
#  cal.fromfat           =  Stool-matched 24hr caloric intake from fats. 


load('your/path/to/230118_preclinicalAD_Ferreiro_etal_DATA.RData')



# CREATE TAXONOMIC PHYLOSEQ OBJECT #############################################

## HOUSEKEEPING ##

## meta: convert appropriate variables to factor--------------------------------
rownames(meta) <- meta$Participant

vars <- c('Participant', 'Sex', 'Race', 'amyloid.positive.AF', 'APOE', 
          'APOE4.status', 'Cancer', 'Active.Depression', 'Hypercholesterolemia',
          'Hypertension', 'Diabetes', 'Liver.Disease', 'Alcohol.Abuse',
          'Tobacco.Use.PastOrPresent', 'Autoimmune.Disorder', 'Thyroid.Disease',
          'Cardiovascular.Disease')

meta[,vars] <- lapply(meta[,vars], factor)


## metaphlan3.raw: format for phyloseq------------------------------------------
rownames(metaphlan3.raw) <- metaphlan3.raw$clade_name
metaphlan3.raw[, c('clade_name', 'NCBI_tax_id')] <- list(NULL)

#take only rows that include s__ (species level annotation)
metaphlan3.species <- metaphlan3.raw[grepl('s__', rownames(metaphlan3.raw)), ]



## m3.tree: trim "GCA--" leading string from tip labels---_---------------------
m3.tree$tip.label <- gsub('GCA_[0-9]+\\|', '', m3.tree$tip.label)



## Create phyloseq taxa table---------------------------------------------------
species <- data.frame(Names = rownames(metaphlan3.species))
species <- data.frame(do.call('rbind', strsplit(as.character(species$Names),'|',fixed=TRUE)))
rownames(species) <- rownames(metaphlan3.species)
colnames(species) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
species <- as.matrix(species)


## Make phyloseq object---------------------------------------------------------
ps <- phyloseq(otu_table(metaphlan3.species, taxa_are_rows = TRUE),
               sample_data(meta),
               tax_table(species), 
               phy_tree(m3.tree))


## Filter low abundance taxa----------------------------------------------------

ps2.filt <- filter_taxa(ps, function(x) mean(x) > 0.1, TRUE)

# This filters out taxa with mean abundance <= 0.1% (Mphln3 relative abundance
# values are out of 100).
# ps2.filt should have 115 taxa across 164 samples, with 64 metadata variables.




# GENERATE TABLE ONE ###########################################################

## Table 1----------------------------------------------------------------------
tab1vars <- c('Age', 'Sex', 'Years.Education', 'Race', 'BMI', 'APOE4.status',
            'Active.Depression', 'Alcohol.Abuse', 'Autoimmune.Disorder', 
            'Cancer', 'Cardiovascular.Disease', 'Diabetes', 'Hypercholesterolemia',
            'Hypertension', 'Liver.Disease', 'Thyroid.Disease', 
            'Tobacco.Use.PastOrPresent','INT_PET', 'INT_CSF', 'INT_MRI', 'INT_CDR')

tab1 <- CreateTableOne(vars=tab1vars, strata='amyloid.positive.AF', data=meta)
tab1Mat <- print(tab1, showAllLevels = TRUE)

#Sex == 1: Male; Sex == 2: Female
#Race == 4: Black; Race == 5: White; Race == 6: Other
#INT_${Assessment}: Time interval in days between assessment and stool collection.
#See Supp. Data File 1 (Data Dictionary) for more information. 


## Sampling Interval Ranges-----------------------------------------------------

# max in years:
INT.max<- meta %>% 
  mutate(across(c(INT_PET, INT_CSF, INT_MRI, INT_CDR), as.numeric))%>%
  summarise(across(c(INT_PET,INT_CSF,INT_MRI,INT_CDR), ~ max(.x, na.rm = TRUE)/365))

# min in days:
INT.min<- meta %>% 
  mutate(across(c(INT_PET, INT_CSF, INT_MRI, INT_CDR), as.numeric))%>%
  summarise(across(c(INT_PET,INT_CSF,INT_MRI,INT_CDR), ~ min(.x, na.rm = TRUE)))




# TAXONOMIC ALPHA DIVERISTY COMPARISONS ########################################

# Alpha diversity metrics should be calculated using the unfiltered abundance 
# data, and require count data rather than fractional relative abundance: so,
# transform unfiltered relative abundances to integer counts, preserving 
# precision to 0.00001%. Function estimate_richness() will throw warning about 
# singletons that may be ignored (more relevant to ASV data).

ps.int <- transform_sample_counts(ps, function(x) trunc(x*100000))
adiv <- estimate_richness(ps.int, measures=c('Observed', 'Shannon'))


# These alpha div vars are already included in 'meta' dataframe, and therefore in 
# the phyloseq object. Going forward appending of re-derived variables to the 
# phyloseq objects will be commented out:

#sample_data(ps2.filt)$Richness <- adiv$Observed
#sample_data(ps2.filt)$Shannon <- adiv$Shannon


## Plot alpha diversity by AD status--------------------------------------------

# Access metadata and gather to long format by alpha diversity metric
ps2.filt.df <- data.frame(sample_data(ps2.filt))
ps2.filt.df.a <- gather(ps2.filt.df, Alpha_Measure, Value, Richness:Shannon, factor_key=TRUE)


p_alpha <- ggplot(ps2.filt.df.a, aes(x=amyloid.positive.AF, y=Value, color=amyloid.positive.AF))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2)+
  scale_color_brewer(palette='Dark2')+
  facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylab('Taxa (MetaPhlAn3)')+
  scale_y_continuous(limits=c(0, NA))+
  scale_x_discrete(labels=c('healthy', 'preclinical'))



## Test for differences in alpha diversity by AD status-------------------------
t.obs <- t.test(Richness ~ amyloid.positive.AF, data = ps2.filt.df)
t.sha <- t.test(Shannon ~ amyloid.positive.AF, data = ps2.filt.df)



# FIRMICUTES/BACTEROIDES RATIO #################################################

## Agglomerate at the phylum level and calculate ratio by individual------------
ps2.filt.phylum <- tax_glom(ps2.filt, taxrank = 'Phylum')
ps2.phylum.df <- psmelt(ps2.filt.phylum)

ps.phylum <- subset(ps2.phylum.df, 
                    select=c('Sample', 'Abundance', 'amyloid.positive.AF', 'Phylum'))

# Spread table wide by phylum and calculate ratio
ps.phylum.wide <- ps.phylum %>% spread(Phylum, Abundance)
ps.phylum.wide$FBratio <- ps.phylum.wide$p__Firmicutes/ps.phylum.wide$p__Bacteroidetes



## Summarize FBratio by group and test for differences--------------------------
# Remove 3 samples with Bacteroidetes == 0 or near 0 to avoid 'Inf' or nonsensical
# ratios. All three samples were from healthy (non preclinical AD) individuals.
ps.phylum.wide2 <- subset(ps.phylum.wide, !is.infinite(FBratio) & FBratio < 1000) 


ps.phylum.FBsumm <- ps.phylum.wide2 %>% group_by(amyloid.positive.AF) %>%
  summarise(ci = list(mean_cl_normal(FBratio) %>% rename(mean=y, lwr=ymin, upr=ymax))) %>%
  unnest(cols = c(ci))

t.FBratio <- t.test(FBratio ~ amyloid.positive.AF, ps.phylum.wide2, 
                    na.action = na.omit) 



# Retrieve FBratio vector:
FBratio <- subset(ps.phylum.wide, select=c('Sample', 'FBratio'))
rownames(FBratio) <- FBratio$Sample
FBratio <- FBratio %>% rename(Participant = Sample)

# Add FBratio to the phyloseq object:
#FBratio.ph <- sample_data(FBratio)
#ps2.filt <- merge_phyloseq(ps2.filt, FBratio.ph)




# STACKED TAXONOMIC BARPLOTS ###################################################

## Agglomerate at the genus level and color/order by phylogeny -----------------
ps.genus <- tax_glom(ps2.filt, taxrank = 'Genus')
ps.genus.df <- psmelt(ps.genus)

# prepare custom color palette for genera
genus_colors <- genus_colors.df$Color
names(genus_colors) <- genus_colors.df$Genus
genusColScale <- scale_fill_manual(values = genus_colors)

# reorder genus levels according to phylogeny (conveniently, this order has  
# been integrated into the custom color palette)
ps.genus.df$Genus <- factor(ps.genus.df$Genus, levels = genus_colors.df$Genus)


## Plot by AD status (abundances summarized by group)---------------------------
p_bar_genus.st <- ggplot(ps.genus.df, aes(x=amyloid.positive.AF, y=Abundance, fill=Genus))+
  genusColScale+
  geom_bar(stat='identity', position='fill')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))+
  scale_x_discrete(labels=c('healthy', 'preclinical'))


## Plot by individual participant-----------------------------------------------

# Order participants by increasing Firmicutes/Bacteroides ratio to improve
# aesthetics and interpretation of the participant-stratified barplots. 
FBratio <- FBratio %>%
  mutate(Participant = fct_reorder(Participant, FBratio, .desc = FALSE)) %>%
  arrange(Participant)

participant.order <- as.character(FBratio$Participant)

ps.genus.df <- ps.genus.df %>%
  mutate(Participant = fct_relevel(Participant, participant.order))


p_bar_genus.st.part <- ggplot(ps.genus.df, aes(x=Participant, y=Abundance, fill=Genus))+
  genusColScale+
  geom_bar(stat='identity', position='fill')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))+
  facet_grid(~amyloid.positive.AF, space='free', scales='free')




# TAXONOMIC ORDINATION ANALYSES ################################################

## Set up ordination wrapper function------------------------------------------- 

# Required packages: phyloseq, vegan, ggplot2, ggExtra
# Arguments: 
#  phyloseq.obj           = phyloseq object
#  method (string)        = ordination method supported by phyloseq (calls vegan)
#  distance (string)      = between sample distance metric supported by phyloseq 
#                           (see vegdist) 
#  adonisformula (string) = formula to pass to adonis2 (PERMANOVA);
#                           LHS should be 'ps.dist'
#  group (obj)            = variable by to which color sample points
#  seed (int)             = seed for adonis test
#  Rpalette (string)      = name of RColorBrewer palette for plotting
#  markershape (int)      = to recreate, enter 19 for metaphlan, 17 for humann
#  saveplot (logical)     = TRUE if figure should be saved to drive (will overwrite)
#                           NOTE! if TRUE, update filepath within function.
#
# Return:                 
#  list(ordination = ps.ordination, dist = ps.dist, plot = ordination.plot2, 
#       axis1test = axis1.test, axis2test = axis2.test, adonis = ps.adonis)



ordinate2.mphln.AF <- function(phyloseq.obj, method, distance, adonisformula, 
                               group, seed, Rpalette, markershape, saveplot) {
  #ARG HOUSEKEEPING
  metadata <- data.frame(phyloseq::sample_data(phyloseq.obj))
  
  group.name <- deparse(substitute(group))
  group.idx <- grep(group.name, colnames(metadata))
  
  
  #ORDINATE
  ps.ordination <- phyloseq::ordinate(phyloseq.obj, method=method, distance=distance)
  
  #PLOT
  ordination.plot <- phyloseq::plot_ordination(phyloseq.obj, ps.ordination,
                                               type="samples", color=group.name)+
    geom_point(size=3, shape=markershape)+
    stat_ellipse()+
    scale_color_brewer(palette=Rpalette)+ 
    theme_classic()+
    theme(axis.text.y = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          legend.text = element_text(size=20),
          legend.position = 'left',
          plot.title = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", size=3))
  
  ordination.plot2 <- ggMarginal(ordination.plot, type='boxplot', groupFill=TRUE, 
                                 size=10)
  
  #SAVE PLOT
  phyloseq.obj.name <- deparse(substitute(phyloseq.obj))
  filepath <- paste('path/for/your/Figure/', 
                    Sys.Date(),"_", 
                    phyloseq.obj.name,"_",
                    group.name,"_",
                    method,"_",
                    distance, ".pdf", sep="")
  
  if (saveplot == TRUE){
    ggsave(filepath, 
           plot = ordination.plot2, 
           device ='pdf', 
           width=20, height=13, units='cm')
  }
  
  #SAMPLE COORDINATES: MARGINAL t (or anova)-tests 
  coordinates <- data.frame(ps.ordination$vectors)
  print('Vector lengths: Coord axis 1, Coord axis 2, metadata$group')
  print(c(length(coordinates$Axis.1), length(coordinates$Axis.2), length(metadata[, group.idx])))
  print('Variable compared')
  print(group.name)
  
  df <- data.frame(Coord1 = coordinates$Axis.1,
                   Coord2 = coordinates$Axis.2,
                   Group = metadata[, group.idx])
  
  axis1.test <- aov(Coord1 ~ Group, data=df)
  axis2.test <- aov(Coord2 ~ Group, data=df)
  
  #ADONIS2 (PERMANOVA) TESTs
  ps.dist <- phyloseq::distance(phyloseq.obj, method=distance) 
  
  ps.adonis <- vegan::adonis2(formula(adonisformula),  
                              data=metadata,
                              na=na.omit,
                              permutations = 10000,
                              subset=complete.cases(ps.ordination))  
  
  #RETURN 
  list(ordination = ps.ordination, dist = ps.dist, plot = ordination.plot2, 
       axis1test = axis1.test, axis2test = axis2.test, 
       adonis = ps.adonis)
  
}



## Execute PCoA on taxonomic abundances, with group comparisons-----------------

adonis.formula <- 'ps.dist ~ Age + APOE4.status + Diabetes + BMI + Hypertension + Interval_Days_Amyloid + amyloid.positive.AF'

ordination3b <- ordinate2.mphln.AF(ps2.filt, method = 'PCoA', distance = 'unifrac', 
                                   adonisformula = adonis.formula,
                                   group = amyloid.positive.AF, seed = 13,
                                   Rpalette = 'Dark2', markershape = 19, 
                                   saveplot = FALSE)

# Add the PCOA1 and PCOA2 coordinates to the phyloseq obj metadata
#sample_data(ps2.filt)$Tax.PCoA1 <- ordination3b$ordination$vectors[,1]
#sample_data(ps2.filt)$Tax.PCoA2 <- ordination3b$ordination$vectors[,2]


## CAP analysis of taxonomic abundances-----------------------------------------

# CAP (see Anderson and Willis, Ecology 2003) has a number of implementations in
# R. The phyloseq implementation calls vegan::capscale, and is equivalent to 
# carrying out vegan::rda on PCoA coordinate vectors. The same 'distance' matrix
# that was calculated for the original PCoA is used, and in this case we fit the
# same model as in the PERMANOVA (adonis2) carried out in the above section.
# One idiosyncrasy is that rather than use the original PCoA coordinate vectors
# (generated in phyloseq by calling ape::pcoa), the phyloseq implementation of 
# CAP re-performs PCoA with vegan::cmdscale, and uses these as input. This has
# negligible impact on results, however.


# Access distance matrix from PCoA ordination and run CAP 
u.dist <- ordination3b$dist
ps2.filt.cap <- ordinate(ps2.filt,
                         method = 'CAP',
                         distance = u.dist,
                         formula = ~Age+APOE4.status+Diabetes+BMI+Hypertension+Interval_Days_Amyloid+amyloid.positive.AF,
                         na.action = na.omit)


# Biplot (out of curiosity)
cap.biplot <- ps2.filt.cap$CCA$biplot

# Plot. Note! phyloseq::plot_ordination tries to add values for % variance
# explained on each of the axes-- but these values are not relevant to CAP, 
# and indeed CAP plots are not shown with any '% variance explained' information
# on the axes. These auto-populated values should be ignored. 
p_ps2.filt.cap <- plot_ordination(ps2.filt, ps2.filt.cap, color = 'amyloid.positive.AF')+
  geom_point(size=3)+
  stat_ellipse()+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = 'left',
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))

p_ps2.filt.cap2 <- ggMarginal(p_ps2.filt.cap, type='boxplot', groupFill=TRUE, size=10)


# Test for the significance of the constraints (model terms)
cap.anova <- anova(ps2.filt.cap, by='term', permutations=10000)


# Access CAP axis coordinates and test for significant differences on CAP1 and CAP2 
ps2.filt.cap.df <- p_ps2.filt.cap$data

cap1.t <- t.test(CAP1 ~ amyloid.positive.AF, data=ps2.filt.cap.df) 
cap2.t <- t.test(CAP2 ~ amyloid.positive.AF, data=ps2.filt.cap.df)


# CREATE PATHWAYS PHYLOSEQ OBJECT ##############################################

## HOUSEKEEPING ##

## meta: source metadata from current phyloseq object (ps2.filt)----------------
meta.psh <- data.frame(sample_data(ps2.filt))

## humann3.raw: format for phyloseq---------------------------------------------
rownames(humann3.raw) <- humann3.raw$Pathway
humann3.raw$Pathway <- NULL

# Omit rows that break down the species associations (contain '|' delimiter)
humann3 <- humann3.raw[!grepl('\\|', rownames(humann3.raw)), ]

# Check pathway relative abundances sum to 1 (unlike taxa, which sum to 100)
humann3.sums <- colSums(humann3)

samples_humann <- colnames(humann3) 
#N = 169, creation of the phyloseq object will subset to those in meta.psh


## Create phyloseq 'taxa' table, adapted for pathways---------------------------
pathways <- data.frame(Pathway=rownames(humann3))
rownames(pathways) <- pathways$Pathway
pathways <- as.matrix(pathways)

## Create pathways phyloseq object----------------------------------------------
psh <- phyloseq(otu_table(humann3, taxa_are_rows = TRUE),
                sample_data(meta.psh),
                tax_table(pathways))

# 164 samples by 485 pathways

## Prune unintegrated/unmapped pathways and renormalize-------------------------
psh <- subset_taxa(psh, Pathway!="UNMAPPED"& Pathway!="UNINTEGRATED") #483 pathways
psh  <- transform_sample_counts(psh, function(x) (x / sum(x))*100 )



## Filter low abundance pathways------------------------------------------------
psh.filt <- filter_taxa(psh, function(x) mean(x) > 0.01, TRUE) 

# This filters out pathways with mean abundance <= 0.01% (Humann3 relative abundance
# values were re-normalized to 100 in the previous section).
# psh.filt should have 324 pathways across 164 samples, with 64 metadata variables.


# PATHWAY ALPHA DIVERSITY COMPARISONS ##########################################

# Alpha diversity metrics should be calculated using the unfiltered abundance 
# data, and require count data rather than fractional relative abundance: so,
# transform unfiltered relative abundances to integer counts, preserving 
# precision to 0.00001%. Function estimate_richness() will throw warning about 
# singletons that may be ignored (more relevant to ASV data).

psh.alpha.int <- transform_sample_counts(psh, function(x) trunc(x*100000))
h.adiv <- estimate_richness(psh.alpha.int, measures=c('Observed', 'Shannon'))

# These alpha div vars were already included in 'meta.psh' dataframe, and therefore in 
# the phyloseq object. Appending of re-derived variables to the phyloseq objects 
# is therefore commented out:
# 
# rownames(h.adiv) <- gsub("X", "", rownames(h.adiv))
# 
# # Add pathway alpha diversities to humann phyloseq object (no need to reorder)
# sample_data(psh.filt)$Fnl.Richness <- h.adiv$Observed
# sample_data(psh.filt)$Fnl.Shannon <- h.adiv$Shannon
# 
# # Add pathway alpha diversities to mphln phyloseq object 
# # Important! Sample order is inverted. This is handled by converting dataframe
# # to phyloseq 'sample_data' object before addition of new metadata. Actually,
# # this is the safe way to add any additional metadata to a phyloseq object.
# 
# h.adiv.ph <- sample_data(h.adiv)
# sample_data(ps2.filt)$Fnl.Richness <- h.adiv.ph$Fnl.Richness
# sample_data(ps2.filt)$Fnl.Shannon <- h.adiv.ph$Fnl.Shannon


## Plot pathway alpha diversity by AD status------------------------------------

# Access metadata and gather to long format by alpha diversity metric
psh.filt.df <- data.frame(sample_data(psh.filt))
psh.filt.df.a <- gather(psh.filt.df, Alpha_Measure, Value, Fnl.Richness:Fnl.Shannon, factor_key=TRUE)


p_alpha_h <- ggplot(psh.filt.df.a, aes(x=amyloid.positive.AF, y=Value, color=amyloid.positive.AF))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2, pch=17)+
  scale_color_brewer(palette='Dark2')+
  facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylab('Pathways (HUMAnN3)')+
  scale_y_continuous(limits=c(0, NA))+
  scale_x_discrete(labels=c('healthy', 'preclinical'))




## Test for differences in pathway alpha diversity by AD status-----------------
t.obs.h <- t.test(Fnl.Richness ~ amyloid.positive.AF, data = psh.filt.df)
t.sha.h <- t.test(Fnl.Shannon ~ amyloid.positive.AF, data = psh.filt.df) 



# PATHWAY ORDINATION ANALYSES ##################################################

## Set up ordination wrapper function-------------------------------------------
# Note: this function is identical to ordinate2.mphln.AF with the addition of the
# 'binary' argument to enable use of binary Bray-Curtis index.
# An intrepid scholar may one day combine these into one streamlined function. 
#
# Required packages: phyloseq, vegan, ggplot2, ggExtra
# Arguments: 
#  phyloseq.obj           = phyloseq object
#  method (string)        = ordination method supported by phyloseq (calls vegan)
#  distance (string)      = between sample distance metric supported by phyloseq 
#                           (see vegdist) 
#  adonisformula (string) = formula to pass to adonis2 (PERMANOVA);
#                           LHS should be 'ps.dist'
#  group (obj)            = variable by to which color sample points
#  seed (int)             = seed for adonis test
#  Rpalette (string)      = name of RColorBrewer palette for plotting
#  markershape (int)      = to recreate, enter 19 for metaphlan, 17 for humann
#  binary (logical)       = TRUE if the vegdist distance metric should be the 
#                           binary variant.
#  saveplot (logical)     = TRUE if figure should be saved to drive (will overwrite)
#                           NOTE! if TRUE, update filepath within function.
#
# Return:                 
#  list(ordination = ps.ordination, dist = ps.dist, plot = ordination.plot2, 
#       axis1test = axis1.test, axis2test = axis2.test, adonis = ps.adonis)


ordinate2.humann.AF <- function(phyloseq.obj, method, distance, adonisformula, 
                                group, seed, Rpalette, markershape, binary, 
                                saveplot) { 
  #ARG HOUSEKEEPING
  metadata <- data.frame(phyloseq::sample_data(phyloseq.obj))
  
  group.name <- deparse(substitute(group))
  group.idx <- grep(group.name, colnames(metadata))
  
  
  #ORDINATE
  ps.ordination <- phyloseq::ordinate(phyloseq.obj, method=method, 
                                      distance=distance, binary = binary)
  
  #PLOT
  ordination.plot <- phyloseq::plot_ordination(phyloseq.obj, ps.ordination, 
                                               type="samples", color=group.name)+
    geom_point(size=3, shape=markershape)+
    stat_ellipse()+
    scale_color_brewer(palette=Rpalette)+ 
    theme_classic()+
    theme(axis.text.y = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          legend.text = element_text(size=20),
          legend.position = 'left',
          plot.title = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", size=3))
  
  ordination.plot2 <- ggMarginal(ordination.plot, type='boxplot', groupFill=TRUE, 
                                 size=10)
  
  #SAVE PLOT
  phyloseq.obj.name <- deparse(substitute(phyloseq.obj))
  filepath <- paste('path/for/your/Figure/', 
                    Sys.Date(),"_", 
                    phyloseq.obj.name,"_",
                    group.name,"_",
                    method,"_",
                    distance, ".pdf", sep="")
  
  if (saveplot == TRUE){
    ggsave(filepath, 
           plot = ordination.plot2, 
           device ='pdf', 
           width=20, height=13, units='cm')
  }
  
  #SAMPLE COORDINATES: MARGINAL t (or anova)-tests (print to check)
  coordinates <- data.frame(ps.ordination$vectors)
  print('Vector lengths: Coord axis 1, Coord axis 2, metadata$group')
  print(c(length(coordinates$Axis.1), length(coordinates$Axis.2), length(metadata[, group.idx])))
  print('Variable compared')
  print(group.name)
  
  df <- data.frame(Coord1 = coordinates$Axis.1,
                   Coord2 = coordinates$Axis.2,
                   Group = metadata[, group.idx])
  
  axis1.test <- aov(Coord1 ~ Group, data=df)
  axis2.test <- aov(Coord2 ~ Group, data=df)
  
  #ADONIS TESTs 
  ps.dist <- phyloseq::distance(phyloseq.obj, method=distance, binary=binary)
  
  ps.adonis <- vegan::adonis2(formula(adonisformula),  
                              data=metadata,
                              na=na.omit,
                              permutations = 10000,
                              subset=complete.cases(ps.ordination))  
  
  #RETURN 
  list(ordination = ps.ordination, dist = ps.dist, plot = ordination.plot2, 
       axis1test = axis1.test, axis2test = axis2.test, 
       adonis = ps.adonis)
  
}



## Execute PCoA on pathway abundances, with group comparisons-------------------

# Using same model formula as before:
adonis.formula <- 'ps.dist ~ Age + APOE4.status + Diabetes + BMI + Hypertension + Interval_Days_Amyloid + amyloid.positive.AF'


h.ordination4b <- ordinate2.humann.AF(psh.filt, method = 'PCoA', distance = 'bray', binary = TRUE, 
                                      adonisformula = adonis.formula,
                                      group = amyloid.positive.AF, seed = 13,
                                      Rpalette = 'Dark2', markershape=17, saveplot = TRUE)


# # Add the pathway PCOA1 and PCOA2 coordinates to the phyloseq objects
# hmn.pcoa.vectors <- data.frame(h.ordination4b$ordination$vectors[,1:2])
# colnames(hmn.pcoa.vectors) <- c('Fnl.PCoA1', 'Fnl.PCoA2')
# hmn.pcoa.vectors.ph <- sample_data(hmn.pcoa.vectors)
# 
# ps2.filt <- merge_phyloseq(ps2.filt, hmn.pcoa.vectors.ph)
# psh.filt <- merge_phyloseq(psh.filt, hmn.pcoa.vectors.ph)



## CAP analysis of pathway abundances-------------------------------------------

# Access distance matrix from PCoA ordination and run CAP 
bb.dist <- h.ordination4b$dist
psh.filt.cap <- ordinate(psh.filt,
                         method = 'CAP',
                         distance = bb.dist,
                         formula = ~Age+APOE4.status+Diabetes+BMI+Hypertension+Interval_Days_Amyloid+amyloid.positive.AF,
                         na.action = na.omit)

# Biplot 
h.cap.biplot <- psh.filt.cap$CCA$biplot


# Plot. Note! phyloseq::plot_ordination tries to add values for % variance
# explained on each of the axes-- but these values are not relevant to CAP, 
# and indeed CAP plots are not shown with any '% variance explained' information
# on the axes. These auto-populated values should be ignored. 
p_psh.filt.cap <- plot_ordination(psh.filt, psh.filt.cap, color = 'amyloid.positive.AF')+
  geom_point(size=3, shape=17)+
  stat_ellipse()+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = 'left',
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))

p_psh.filt.cap2 <- ggMarginal(p_psh.filt.cap, type='boxplot', groupFill=TRUE, size=10)

# Test for the significance of the constraints (model terms) 
h.cap.anova <- anova(psh.filt.cap, by='term', permutations=10000)

# Access CAP axis coordinates and test for significant differences on CAP1 and CAP2 
psh.filt.cap.df <- p_psh.filt.cap$data

h.cap1.t <- t.test(CAP1 ~ amyloid.positive.AF, data=psh.filt.cap.df)
h.cap2.t <- t.test(CAP2 ~ amyloid.positive.AF, data=psh.filt.cap.df)




# PAIRWISE BIOMARKER CORRELATION ANALYSES ######################################

## Prepare biomarker data-------------------------------------------------------

# Access sample metadata, including GM metrics of alpha diversity and PCoA 
# coordinates 
metadata.corr <- data.frame(sample_data(ps2.filt))

# Subset to biomarkers and GM metrics of interest
corr.vars <- c('Tax.PCoA1', 'Tax.PCoA2', 'Fnl.PCoA1', 'Fnl.PCoA2',
               'combined.Centiloid', 'LUMIPULSE_CSF.ab42.ab40.ratio', 
               'Tauopathy', 'LUMIPULSE_CSF_pTau',
               'LUMIPULSE_CSF_tTau', 'CortSig_Thickness', 'MR_LV_RV_HIPPOCAMPUS', 
               'WMH_volume', 
               'Polygenic.Risk.Score', 'APOE4.status') 

metadata.corr2 <- metadata.corr[, corr.vars]

# Numerically encode APOE e4 carrier status.
metadata.corr2 <- metadata.corr2 %>% 
  mutate(APOE4.status = if_else(APOE4.status == 'e4-', 0, 1))


## Determine pairwise spearman correlations and plot correlogram----------------
corr.spearman <- rcorr(as.matrix(metadata.corr2, type='spearman'))

# Convert symmetric p-value matrix to non-redundant vector, adjust p-values, and 
# return to named symmetric matrix. Convert diagonal NAs to 0 for visualization. 
corr.spearman.p <- sm2vec(corr.spearman$P, diag = FALSE)
idx <- which(corr.spearman.p < 0.05)
corr.spearman.p[idx] <- p.adjust(corr.spearman.p[idx], method='BH')
corr.spearman.p.adj <- vec2sm(corr.spearman.p, diag = FALSE)
colnames(corr.spearman.p.adj) <- colnames(corr.spearman$P)
rownames(corr.spearman.p.adj) <- rownames(corr.spearman$P)
corr.spearman.p.adj[is.na(corr.spearman.p.adj)] <- 0


# Plot correlogram
p_corr.spearman <- corrplot(corr.spearman$r, order = 'original', p.mat=corr.spearman.p.adj,
                             sig.level=0.05, addCoef.col='black', insig='blank')


# BIOMARKER REGRESSION ANALYSES ################################################

## PET amyloid------------------------------------------------------------------
# Taxonomic PCoA1 v PET amyloid (combined.Centiloid)
p_pcoa1_amyloid <- ggplot(metadata.corr, aes(x=Tax.PCoA1, y=combined.Centiloid, 
                                             color=amyloid.positive.AF))+
  geom_smooth(method='lm', show.legend = FALSE)+
  geom_point(show.legend = FALSE)+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylim(0, NA)+
  ylab('PET amyloid')

# Fit null model (PET amyloid as explained by AD status only)
m.pcoa1_amyloid0b <- lm(combined.Centiloid ~ amyloid.positive.AF, 
                        metadata.corr, na.action=na.omit)
# Include interaction with Tax.PCoA1  
m.pcoa1_amyloid <- lm(combined.Centiloid ~ Tax.PCoA1*amyloid.positive.AF, 
                      metadata.corr, na.action=na.omit) 
# Test improvement in model fit (not significant)
aov.cent_pcoa1 <- anova(m.pcoa1_amyloid0b, m.pcoa1_amyloid)



# Taxonomic PCoA2 v PET amyloid (combined.Centiloid)
p_pcoa2_amyloid <- ggplot(metadata.corr, aes(x=Tax.PCoA2, y=combined.Centiloid,
                                             color=amyloid.positive.AF))+
  geom_smooth(method='lm', show.legend = FALSE)+
  geom_point(show.legend = FALSE)+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylim(0, NA)+
  ylab('PET amyloid')

# Fit null model (PET amyloid as explained by AD status only, same as above)
m.pcoa2_amyloid0b <- lm(combined.Centiloid ~ amyloid.positive.AF, 
                        metadata.corr, na.action=na.omit)
# Include interaction with Tax.PCoA2
m.pcoa2_amyloid <- lm(combined.Centiloid ~ Tax.PCoA2*amyloid.positive.AF, 
                      metadata.corr, na.action=na.omit) 
# Test improvement in model fit (p < 0.1)
aov.cent_pcoa2 <- anova(m.pcoa2_amyloid0b, m.pcoa2_amyloid)



# Functional PCoA1 v PET amyloid (combined.Centiloid) 
p_fnl.pcoa1_amyloid <- ggplot(metadata.corr, aes(x=Fnl.PCoA1, y=combined.Centiloid, 
                                                 color=amyloid.positive.AF))+
  geom_smooth(method='lm', show.legend = FALSE)+
  geom_point(show.legend = FALSE, shape=17)+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylim(0, NA)+
  ylab('PET amyloid')

# Fit null model (PET amyloid as explained by AD status only, same as above)
m.fnl.pcoa1_amyloid0b <- lm(combined.Centiloid ~ amyloid.positive.AF, 
                            metadata.corr, na.action=na.omit)
# Include interaction with Fnl.PCoA1
m.fnl.pcoa1_amyloid <- lm(combined.Centiloid ~ Fnl.PCoA1*amyloid.positive.AF, 
                          metadata.corr, na.action=na.omit) 
# Test improvement in model fit (p < 0.1)
aov.cent_fnl.pcoa1 <- anova(m.fnl.pcoa1_amyloid0b, m.fnl.pcoa1_amyloid)



## PET tau----------------------------------------------------------------------
# Taxonomic PCoA1 v PET tau
p_pcoa1_tau <- ggplot(metadata.corr, aes(x=Tax.PCoA1, y=Tauopathy, 
                                         color=amyloid.positive.AF))+
  geom_smooth(method='lm', show.legend = FALSE)+
  geom_point(show.legend = FALSE)+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylim(0.5, NA)+
  ylab('PET Tau')

# Fit null model (PET tau as explained by AD status only)
m.pcoa1_tau0b <- lm(Tauopathy ~ amyloid.positive.AF, 
                    metadata.corr, na.action=na.omit)
# Include interaction with Tax.PCoA1
m.pcoa1_tau <- lm(Tauopathy ~ Tax.PCoA1*amyloid.positive.AF, 
                  metadata.corr, na.action=na.omit)
# Test improvement in model fit (p < 0.1)
aov.tau_pcoa1 <- anova(m.pcoa1_tau0b, m.pcoa1_tau)



# Taxonomic PCoA2 v PET tau
p_pcoa2_tau <- ggplot(metadata.corr, aes(x=Tax.PCoA2, y=Tauopathy, 
                                         color=amyloid.positive.AF))+
  geom_smooth(method='lm', show.legend = FALSE)+
  geom_point(show.legend = FALSE)+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylim(0.5, NA)+
  ylab('PET Tau')

# Fit null model (PET tau as explained by AD status only, same as above)
m.pcoa2_tau0b <- lm(Tauopathy ~ amyloid.positive.AF, 
                    metadata.corr, na.action=na.omit)
# Include interaction with Tax.PCoA2
m.pcoa2_tau <- lm(Tauopathy ~ Tax.PCoA2*amyloid.positive.AF, 
                  metadata.corr, na.action=na.omit)
# Test improvement in model fit (p < 0.05)
aov.tau_pcoa2 <- anova(m.pcoa2_tau0b, m.pcoa2_tau)



# Functional PCoA1 v PET tau
p_fnl.pcoa1_tau <- ggplot(metadata.corr, aes(x=Fnl.PCoA1, y=Tauopathy,
                                             color=amyloid.positive.AF))+
  geom_smooth(method='lm', show.legend = FALSE)+
  geom_point(show.legend = FALSE, shape=17)+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylim(0.5, NA)+
  ylab('PET Tau')

# Fit null model (PET tau as explained by AD status only, same as above)
m.fnl.pcoa1_tau0b <- lm(Tauopathy ~ amyloid.positive.AF, 
                        metadata.corr, na.action=na.omit)
# Include interaction with Fnl.PCoA1
m.fnl.pcoa1_tau <- lm(Tauopathy ~ Fnl.PCoA1*amyloid.positive.AF, 
                      metadata.corr, na.action=na.omit)
# Test improvement in model fit (not significant)
aov.tau_fnl.pcoa1 <- anova(m.fnl.pcoa1_tau0b, m.fnl.pcoa1_tau)


## Neurodegeneration------------------------------------------------------------
# Taxonomic PCoA1 v Cortical Signature
p_pcoa1_cortsig <- ggplot(metadata.corr, aes(x=Tax.PCoA1, y=CortSig_Thickness,
                                             color=amyloid.positive.AF))+
  geom_smooth(method='lm', show.legend = FALSE)+
  geom_point(show.legend = FALSE)+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylim(2, NA)+
  ylab('Cortical Signature')

# Fit null model (Cortical Signature as explained by AD status only)
m.pcoa1_cort0b <- lm(CortSig_Thickness ~ amyloid.positive.AF, 
                     metadata.corr, na.action=na.omit)
# Include interaction with Tax.PCoA1
m.pcoa1_cort <- lm(CortSig_Thickness ~ Tax.PCoA1*amyloid.positive.AF, 
                   metadata.corr, na.action=na.omit)
# Test improvement in model fit (not significant)
aov.cort_pcoa1 <- anova(m.pcoa1_cort0b, m.pcoa1_cort)



# Taxonomic PCoA2 v Cortical Signature
p_pcoa2_cortsig <- ggplot(metadata.corr, aes(x=Tax.PCoA2, y=CortSig_Thickness, 
                                             color=amyloid.positive.AF))+
  geom_smooth(method='lm', show.legend = FALSE)+
  geom_point(show.legend = FALSE)+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylim(2, NA)+
  ylab('Cortical Signature')

# Fit null model (Cortical Signature as explained by AD status only, same as above)
m.pcoa2_cort0b <- lm(CortSig_Thickness ~ amyloid.positive.AF, 
                     metadata.corr, na.action=na.omit)
# Include interaction with Tax.PCoA2
m.pcoa2_cort <- lm(CortSig_Thickness ~ Tax.PCoA2*amyloid.positive.AF, 
                   metadata.corr, na.action=na.omit)
# Test improvement in model fit (not significant)
aov.cort_pcoa2 <- anova(m.pcoa2_cort0b, m.pcoa2_cort)


# Functional PCoA1 v Cortical Signature
p_fnl.pcoa1_cortsig <- ggplot(metadata.corr, aes(x=Fnl.PCoA1, y=CortSig_Thickness, 
                                                 color=amyloid.positive.AF))+
  geom_smooth(method='lm', show.legend = FALSE)+
  geom_point(show.legend = FALSE, shape=17)+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=20),
        panel.border = element_rect(fill=NA, colour="black", size=2))+
  ylim(2, NA)+
  ylab('Cortical Signature')

# Fit null model (Cortical Signature as explained by AD status only, same as above)
m.fnl.pcoa1_cort0b <- lm(CortSig_Thickness ~ amyloid.positive.AF, 
                         metadata.corr, na.action=na.omit)
# Include interaction with Fnl.PCoA1
m.fnl.pcoa1_cort <- lm(CortSig_Thickness ~ Fnl.PCoA1*amyloid.positive.AF,
                       metadata.corr, na.action=na.omit)
# Test improvement in model fit (not significant)
aov.cort_fnl.pcoa1 <- anova(m.fnl.pcoa1_cort0b, m.fnl.pcoa1_cort)



## Summarize models for reporting-----------------------------------------------
# Taxonomic PCoA1 
reg_models.pcoa1 <- list(m.pcoa1_amyloid,
                         m.pcoa1_tau,
                         m.pcoa1_cort)

modelsummary(reg_models.pcoa1, 
             fmt=3,
             estimate = "{estimate} ({std.error}){stars}",
             statistic = 'conf.int',
             conf_level= .95,
             output='path/to/your/file/RegModelSummaries_pcoa1.docx')

# Taxonomic PCoA2
reg_models.pcoa2 <- list(m.pcoa2_amyloid,
                         m.pcoa2_tau,
                         m.pcoa2_cort)

modelsummary(reg_models.pcoa2, 
             fmt=3,
             estimate = "{estimate} ({std.error}){stars}",
             statistic = 'conf.int',
             conf_level= .95,
             output='path/to/your/file/RegModelSummaries_pcoa2.docx')

# Functional (pathway) PCoA1
reg_models.fnl.pcoa1 <- list(m.fnl.pcoa1_amyloid,
                             m.fnl.pcoa1_tau,
                             m.fnl.pcoa1_cort)

modelsummary(reg_models.fnl.pcoa1, 
             fmt=3,
             estimate = "{estimate} ({std.error}){stars}",
             statistic = 'conf.int',
             conf_level= .95,
             output='path/to/your/file/RegModelSummaries_fnl-pcoa1.docx')





# SPECIES ASSOCIATIONS WITH AD STATUS (MAASLIN2) ###############################
# Useful to specify a separate output folder for Maaslin analyses: 
# setwd('path/to/your/output/folder')

## Prepare data and fit negative binomial model----------------------------------

# Use taxa abundances from phyloseq object ps2.filt, which has already been 
# filtered to omit taxa with < 0.1% mean relative abundance across all samples. 
# To fit a negative binomial model, count data is required. Here we transform
# relative abundances to integer count data, preserving accuracy to 0.00001%.  
ps2.filt.int <- transform_sample_counts(ps2.filt, function(x) trunc(x*100000))

# Access count-transformed species abundances, and sample metadata 
species.ps2.int <- data.frame(otu_table(ps2.filt.int), check.names = FALSE) 
metadata.ps2.int <- data.frame(sample_data(ps2.filt.int), check.names = FALSE)


#A run Maaslin at species level, with pre-clinical status --
fit_Maaslin2 <- Maaslin2(
  input_data = species.ps2.int,
  input_metadata = metadata.ps2.int,
  normalization = 'NONE',      # Already normalized (count-transformed relabund)
  standardize = 'TRUE',        # We'd like to z-score continuous metadata
  min_prevalence = 0,          # No filtering on prevalence
  min_abundance = 0,           # Data already filtered on abundance
  transform = 'NONE',          # Further transform not needed for NEGBIN model
  analysis_method = 'NEGBIN',
  max_significance = 0.05,     # For q-values (BH-adjusted p-values). Default is 0.25.
  output = "YYMMDD_Mphln3_NEGBIN",
  fixed_effects = c("Age", "APOE4", "Diabetes", "BMI", 
                    "Hypertension", "Interval_Days", 
                    "amyloid.positive.AF")
)


## Filter significant results---------------------------------------------------

# Maaslin will generate a file 'significant_results.tsv' in the indicated output
# folder. Edit externally to replace '.' in taxa names with '|' and save as csv. 
# This is to maintain nomenclature consistency with the original phyloseq object.
# (Maaslin automatically replaces non-alphanumerics). The resulting file can be 
# read back in as a data frame. 

# diffspecies.all.df <- read.csv('path/to/edited/significant_results.csv')

# FOR CONVENIENCE, diffspecies.all.df HAS ALREADY BEEN INCLUDED IN THE 
# ACCOMPANYING WORKSPACE (preclinADAnalyses.RData). 

# Subset to features (taxa) significantly associated with AD status, detected in
# at least 25 samples, and with magnitude of the coefficient > 0.15.
# Rearrange features in order of descending coefficient. 
diffspecies.a.df <- subset(diffspecies.all.df, metadata == 'amyloid.positive.AF' & 
                             N.not.0 > 25 & abs(coef) > 0.15)

diffspecies.a.df <- diffspecies.a.df %>% arrange(desc(coef))



## Generate coefficient plot----------------------------------------------------

# Separate out taxa names at species level and reorder taxa by coefficient value.
diffspecies.a.df2 <- diffspecies.a.df %>% separate(feature, c(NA, 'species'), '\\|s__')
diffspecies.a.df2 <- diffspecies.a.df2 %>% mutate(species = fct_reorder(species, coef))


# Calculate 95% CI for coefficients
diffspecies.a.df2$coefmax <- diffspecies.a.df2$coef + diffspecies.a.df2$stderr
diffspecies.a.df2$coefmin <- diffspecies.a.df2$coef - diffspecies.a.df2$stderr

diffspecies.a.df2$CI95upper <- diffspecies.a.df2$coef + 1.96*diffspecies.a.df2$stderr
diffspecies.a.df2$CI95lower <- diffspecies.a.df2$coef - 1.96*diffspecies.a.df2$stderr

# Plot coefficients
p_speciescoeff_A <- ggplot(diffspecies.a.df2, aes(x=coef, y=species, 
                                                  color=as.factor(coef<0)))+
  geom_point()+
  geom_errorbarh(aes(xmax = CI95upper, xmin = CI95lower, height =0.1))+
  geom_vline(xintercept = 0, linetype='dotted', color='grey', size=1.5)+
  scale_color_manual(values=c("#d95f02", "#1b9e77"))+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))+
  xlim(-1.0, 1.0)+
  xlab('coefficient')




## Access abundance and prevalence data for significant taxa--------------------

# Melt the phyloseq object into a data frame and subset to significant taxa
ps2.filt.int.df <- psmelt(ps2.filt.int) 
ps2.filt.int.a.df <- subset(ps2.filt.int.df, OTU %in% diffspecies.a.df$feature)


# Re-order taxa according to coefficient value.
# This order is preserved in diffspecies.a.df, which was previously ordered: 
diffspecies.a <- unique(diffspecies.a.df$feature)
ps2.filt.int.a.df$OTU <- as.factor(ps2.filt.int.a.df$OTU)
ps2.filt.int.a.df <- ps2.filt.int.a.df %>% mutate(OTU = fct_relevel(OTU, diffspecies.a))



# Make named facet labels vector (human-readable taxa names for abundance plots) 
# That is: split out just the species names, but name the list entries with the 
# full taxon names.
diffspecies.a.OTU <- data.frame(OTU = unique(ps2.filt.int.a.df$OTU))
diffspecies.a.labels.df <- diffspecies.a.OTU %>% separate(OTU, c(NA, 'species'), '\\|s__')
diffspecies.a.labels <- diffspecies.a.labels.df$species
names(diffspecies.a.labels) <- diffspecies.a.OTU$OTU




## Generate taxa prevalence plot------------------------------------------------

# For each significant taxon in each sample, set AbunNotZero == 1 if its 
# abundance is greater than 0. 
ps2.filt.int.a.df$AbunNotZero <- 1*(ps2.filt.int.a.df$Abundance > 0)

# Summarize species prevalence by preclinical AD status group.
ps2.filt.int.a.n0 <- ps2.filt.int.a.df %>% 
  group_by(OTU, amyloid.positive.AF) %>%
  summarise(N = n(), SumNotZero = sum(AbunNotZero)) %>%
  mutate(PercentNotZero = 100*(SumNotZero / N)) %>%
  ungroup()

# Separate out species name from full taxon name. 
ps2.filt.int.a.n0 <- ps2.filt.int.a.n0 %>% separate(OTU, c(NA, 'species'), '\\|s__')
ps2.filt.int.a.n0$species <- factor(ps2.filt.int.a.n0$species, 
                                    levels=unique(ps2.filt.int.a.n0$species))

# Plot prevalence by preclinical AD status group.
p_speciesPnot0_A_tile <- ggplot(ps2.filt.int.a.n0, aes(x=amyloid.positive.AF,
                                                       y=species, fill=PercentNotZero))+
  geom_tile()+
  scale_fill_distiller(palette ='RdBu', type='seq', breaks = c(0,50,100), 'Prevalence (%)')+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))+
  scale_y_discrete(limits=rev)



## Generate taxa abundance plots------------------------------------------------

# For brevity, we are only interested in plotting the top five preclinical AD 
# associated species, and the top five healthy status associated species:
# From diffspecies.a.df2 these are:  

top5species <- c('s__Bacteroides_intestinalis', 's__Eubacterium_ventriosum', 
                 's__Dorea_formicigenerans', 's__Blautia_obeum', 
                 's__Methanobrevibacter_smithii',
                 's__Bacteroides_caccae', 's__Bifidobacterium_longum', 
                 's__Bacteroides_massiliensis', 's__Bacteroides_salyersiae', 
                 's__Bacteroides_faecis')

ps2.filt.int.a.top5 <- subset(ps2.filt.int.a.df, Species %in% top5species)

# Plot abundances. Note transformation of abundances back to percent scale, with
# addition of a small pseudo-count to enable plotting on a log scale:  
p_diffspecies.A5 <- ggplot(ps2.filt.int.a.top5, 
                           aes(x=amyloid.positive.AF, y=Abundance/100000+0.0001))+
  geom_violin(aes(color=amyloid.positive.AF), draw_quantiles = c(0.25,0.75))+
  geom_jitter(aes(color=amyloid.positive.AF),width=0.2, size=1, alpha=0.8, pch=21)+
  scale_color_brewer(palette ='Dark2')+
  theme_classic()+
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = "none",
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=1))+
  facet_wrap(~OTU, labeller = labeller(OTU = diffspecies.a.labels), nrow=2)+
  scale_y_continuous(trans="log10", labels=comma) +
  ylab('Abundance (%)')




# PATHWAY ASSOCIATIONS WITH AD STATUS (MAASLIN2) ###############################

## Prepare data and fit negative binomial model---------------------------------

# Use pathway abundances from phyloseq object psh.filt, which has already been 
# filtered to omit paths with < 0.01% mean relative abundance across all samples. 
# To fit a negative binomial model, count data is required. Here we transform
# relative abundances to integer count data, preserving accuracy to 0.00001%.  
psh.filt.int <- transform_sample_counts(psh.filt, function(x) trunc(x*100000))

# Access count-transformed pathway abundances, and sample metadata 
path.psh.int <- data.frame(otu_table(psh.filt.int), check.names = FALSE) 
metadata.psh.int <- data.frame(sample_data(psh.filt.int), check.names = FALSE)

#A run Maaslin at species level, with pre-clinical status --
h.fit_Maaslin2 <- Maaslin2(
  input_data = path.psh.int,
  input_metadata = metadata.psh.int,
  normalization = 'NONE',    # Already normalized (count-transformed relabund)
  standardize = 'TRUE',      # We'd like to z-score continuous metadata
  min_prevalence = 0,        # No filtering on prevalence
  min_abundance = 0,         # Data already filtered on abundance
  transform = 'NONE',        # Further transform not needed for NEGBIN model
  analysis_method = 'NEGBIN',
  max_significance = 0.05,   # For q-values (BH-adjusted p-values). Default is 0.25
  output = "YYMMDD_Hmn3_NEGBIN",
  fixed_effects = c("Age", "APOE4", "Diabetes", "BMI", 
                    "Hypertension", "Interval_Days", 
                    "amyloid.positive.AF")
)



## Filter significant results---------------------------------------------------

# Maaslin will generate a file 'significant_results.tsv' in the indicated output
# folder. Edit externally to replace '..' separating MetaCyc path IDs and 
# corresponding descriptors with ':' and save as csv. Again, this is to help  
# maintain nomenclature consistency with the original phyloseq object. (Maaslin  
# automatically replaces non-alphanumerics). The resulting file can be read 
# back in as a data frame. 

# diffpath.all.df <- read.csv('path/to/edited/significant_results.csv')

# FOR CONVENIENCE, diffpath.all.df HAS ALREADY BEEN INCLUDED IN THE  
# ACCOMPANYING WORKSPACE (preclinADAnalyses.RData). 

# Subset to features (pathways) significantly associated with AD status, detected 
# in at least 25 samples, and with magnitude of the coefficient > 0.15.
# Rearrange features in order of descending coefficient. 
diffpath.a.df <- subset(diffpath.all.df, metadata == 'amyloid.positive.AF' & 
                          N.not.0 > 25 & abs(coef) > 0.15)

diffpath.a.df <- diffpath.a.df %>% arrange(desc(coef))


## Pre-process pathway names and data, and prepare label vectors----------------
# Pathway names have more idiosyncrasies and require a little extra parsing 
# compared to taxa names, we'll do that here, relying on unique pathway IDs. 

# Duplicate pathway names and separate out the MetaCyc Path ID.
diffpath.a.df$feature2 <- diffpath.a.df$feature
diffpath.a.df <- diffpath.a.df %>% 
  separate(feature2, c("PathID", "PATHNAME.err"), sep=":", extra='merge')

# Melt the phyloseq object into a data frame and similarly separate out the
# MetaCyc Path ID.
psh.filt.int.df <- psmelt(psh.filt.int) 
psh.filt.int.df$Paths <- psh.filt.int.df$OTU
psh.filt.int.df <- psh.filt.int.df %>% 
  separate(Paths, c("PathID", "PATHNAME"), sep=":", extra='merge')

# Subset phyloseq object to significant pathways of interest given filtering 
# criteria. 
psh.filt.int.a.df <- subset(psh.filt.int.df, PathID %in% diffpath.a.df$PathID)

# Check to make sure psh.filt.int.a.df and diffpath.a.df have the same number
# of unique pathways (N = 46). If not, this may indicate issue with name parsing. 
#unique(psh.filt.int.a.df$OTU)

# Add original, correct pathway names ($OTU) to diffpath.a.df 
pathnames <- unique(subset(psh.filt.int.a.df, select = c("PathID", "OTU")))
diffpath.a.df <- merge(diffpath.a.df, pathnames, by='PathID')
diffpath.a.df$PATHNAME.err <- NULL

# Order pathways in melted phyloseq object (which contains pathway abundances) 
# according to descending coefficient values, using order preserved in 
# diffpath.a.df. This is for plotting pathway prevalences/abundances.
diffpath.otu <- unique(diffpath.a.df$OTU)
psh.filt.int.a.df$OTU <- as.factor(psh.filt.int.a.df$OTU)
psh.filt.int.a.df$OTU <- fct_relevel(psh.filt.int.a.df$OTU, levels(diffpath.otu))

# Generate pathway label strings with new lines to enable text wrapping in
# faceted figure panels. 
psh.filt.int.a.df$Pathwrap <- unlist(lapply(strwrap(psh.filt.int.a.df$OTU, 
                                                    width=25, simplify=FALSE), 
                                            paste, collapse="\n"))

# Make named facet labels vector (human-readable pathway names with text 
# wrapping). The list entries are named with the full pathway ($OTU) names.
diffpath.a.OTU <- unique(psh.filt.int.a.df$OTU)
diffpath.a.labels <- unique(psh.filt.int.a.df$Pathwrap)
names(diffpath.a.labels) <- diffpath.a.OTU



## Generate coefficient plot----------------------------------------------------

# Reorder pathways (OTU) by coefficient values. 
diffpath.a.df <- diffpath.a.df %>% mutate(OTU= fct_reorder(OTU, coef))

# Calculate 95% CI for coefficients
diffpath.a.df$coefmax <- diffpath.a.df$coef + diffpath.a.df$stderr
diffpath.a.df$coefmin <- diffpath.a.df$coef - diffpath.a.df$stderr

diffpath.a.df$CI95upper <- diffpath.a.df$coef + 1.96*diffpath.a.df$stderr
diffpath.a.df$CI95lower <- diffpath.a.df$coef - 1.96*diffpath.a.df$stderr

# Plot coefficients
p_pathcoeff_A <- ggplot(diffpath.a.df, aes(x=coef, y=OTU, color=as.factor(coef<0)))+
  geom_point()+
  geom_errorbarh(aes(xmax = CI95upper, xmin = CI95lower, height =0.1))+
  geom_vline(xintercept = 0, linetype='dotted', color='grey', size=1.5)+
  scale_color_manual(values=c("#d95f02", "#1b9e77"))+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))+
  xlab('coefficient')



## Generate pathway prevalence plot---------------------------------------------

# For each significant pathway in each sample, set AbunNotZero == 1 if its 
# abundance is greater than 0. 
psh.filt.int.a.df$AbunNotZero <- 1*(psh.filt.int.a.df$Abundance > 0)

# Summarize pathway prevalence by preclinical AD status group.
psh.filt.int.a.n0 <- psh.filt.int.a.df %>% 
  group_by(OTU, amyloid.positive.AF) %>%
  summarise(N = n(), SumNotZero = sum(AbunNotZero)) %>%
  mutate(PercentNotZero = 100*(SumNotZero / N)) %>%
  ungroup()

# Plot prevalence by preclinical AD status group.
p_pathPnot0_A_tile <- ggplot(psh.filt.int.a.n0, 
                             aes(x=amyloid.positive.AF, y=OTU, fill=PercentNotZero))+
  geom_tile()+
  scale_fill_distiller(palette ='RdBu', type='seq', breaks = c(0,50,100), 'Prevalence (%)')+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))



## Generate pathway abundance plots---------------------------------------------

# For brevity, we are only interested in plotting the top five preclinical AD 
# associated pathways, and the top five healthy status associated pathways:
# From diffpath.otu these are:  

top5path <- c(rev(tail(levels(diffpath.otu),5)), head(levels(diffpath.otu), 5))

psh.filt.int.a.top5 <- subset(psh.filt.int.a.df, OTU %in% top5path)
psh.filt.int.a.top5$OTU <- fct_relevel(psh.filt.int.a.top5$OTU, top5path)

# Plot abundances. Note transformation of abundances back to percent scale, with
# addition of a small pseudo-count to enable plotting on a log scale:  
p_diffpath.A5 <- ggplot(psh.filt.int.a.top5, 
                        aes(x=amyloid.positive.AF, y=Abundance/100000+0.0001))+
  geom_violin(aes(color=amyloid.positive.AF), draw_quantiles = c(0.25,0.75))+
  geom_jitter(aes(color=amyloid.positive.AF),width=0.2, size=1, alpha=0.8, pch=21)+
  scale_color_brewer(palette ='Dark2')+
  theme_classic()+
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = "none",
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=1))+
  facet_wrap(~OTU, labeller = labeller(OTU = diffpath.a.labels), nrow=2)+
  scale_y_continuous(trans="log10", labels=comma) +
  ylab('Abundance (%)')




# DIET - NUTRITIONAL INTAKE COMPARISONS ########################################

## Compare nutritional intake by group (Student's t)----------------------------

# Nutrients: 
# Using variable 'diet.intake', which was loaded in at the top with the main 
# workspace. Subset to only look at percent daily value intake.
diet.intake$Nutrients <- as.factor(diet.intake$Nutrients)

per.intake.df <- subset(diet.intake, select=c('Nutrients', 'PercentDV', 'amyloid.positive.AF'))

# For each nutrient category, carry out a t-test to compare preclinical AD and
# healthy cohort intakes. 
per.t.tests <- per.intake.df %>% group_by(Nutrients, amyloid.positive.AF) %>%
  nest() %>%
  spread(key = amyloid.positive.AF, value = data) %>%
  mutate(t_test = map2(healthy, preclinical, ~{t.test(.x$PercentDV, .y$PercentDV) %>% tidy()}),
         A = map(healthy, nrow),
         B = map(preclinical, nrow)) %>%
  unnest(cols = c(Nutrients, t_test, A, B))

# Adjust for multiple hypothesis testing.
per.t.tests$p.value.BH <- p.adjust(per.t.tests$p.value, method='BH')
per.t.tests.sig <- subset(per.t.tests, p.value.BH < 0.05) #empty because there are none




# Total Calories by source:
# Using variable 'total.cal', which was loaded in at the top with the main 
# workspace. 

total.cal2 <- subset(total.cal, select=c('Source.of.Total.Calories', 'Percent.Total.Calories', 'amyloid.positive.AF'))

total.cal.t.tests <- total.cal2 %>% group_by(Source.of.Total.Calories, amyloid.positive.AF) %>%
  nest() %>%
  spread(key = amyloid.positive.AF, value = data) %>%
  mutate(t_test = map2(healthy, preclinical, ~{t.test(.x$Percent.Total.Calories, .y$Percent.Total.Calories) %>% tidy()}),
         A = map(healthy, nrow),
         B = map(preclinical, nrow)) %>%
  unnest(cols = c(Source.of.Total.Calories, t_test, A, B))


total.cal.t.tests$p.value.BH <- p.adjust(total.cal.t.tests$p.value, method='BH')
total.cal.t.tests.sig <- subset(total.cal.t.tests, p.value.BH < 0.05) #empty because there are none




# Calories by Fat:
# Using variable 'cal.fromfat', which was loaded in at the top with the main
# workspace. 

cff.df2 <- subset(cal.fromfat, select=c('Source.of.Calories.from.Fat', 'Percent.Calories.From.Fat', 'amyloid.positive.AF'))

cff.t.tests <- cff.df2 %>% group_by(Source.of.Calories.from.Fat, amyloid.positive.AF) %>%
  nest() %>%
  spread(key = amyloid.positive.AF, value = data) %>%
  mutate(t_test = map2(healthy, preclinical, ~{t.test(.x$Percent.Calories.From.Fat, .y$Percent.Calories.From.Fat) %>% tidy()}),
         A = map(healthy, nrow),
         B = map(preclinical, nrow)) %>%
  unnest(cols = c(Source.of.Calories.from.Fat, t_test, A, B))


cff.t.tests$p.value.BH <- p.adjust(cff.t.tests$p.value, method='BH')
cff.t.tests.sig <- subset(cff.t.tests, p.value.BH < 0.05) #empty because there are none.



## Plot nutritional intake by group---------------------------------------------

# FUNCTION: makeDietBoxplot will plot boxplots of intake by nutrient,
# optionally subsetting to provided nutritional category, and comparing 
# between preclinical AD and healthy cohorts.
# Required packages: ggplot2
# Arguments:
#  data (dataframe)                  = nutritional intake by participant
#  subset (logical)                  = subset data to the indicated nutrient
#                                      category?
#  subset.category (str)             = if subset == TRUE, use this category 
#  wrap.var (str)                    = column name to wrap figure facets by
#  y.val (str)                       = column name to plot on y axis
#  nrows (int)                       = num rows for wrapping 
# Return:
#  Boxplot by nutrient, for the indicated nutritional category, comparing 
#  preclinical AD and healthy cohorts. 

makeDietBoxplot <- function(data, subset, subset.category, wrap.var, y.val, nrows){
  
  if (subset == TRUE){
    data2 <- subset(data, Nutrient.Category == subset.category)
    print(dim(data2))
  } else {
    data2 <- data
  }
  
  diet.boxplot <- ggplot(data2, 
                         aes_string(x='amyloid.positive.AF', y=y.val, 
                                    color='amyloid.positive.AF'))+
    geom_hline(yintercept=100, color='darkgray')+
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(width=0.2, pch=19)+
    scale_color_brewer(palette='Dark2')+
    facet_wrap(as.formula(paste('~', wrap.var)), nrow=nrows, scales = 'free_y')+
    theme_classic()+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(size=10),
          axis.text.x = element_blank(),
          strip.text.x = element_text(size=10),
          panel.border = element_rect(fill=NA, colour="white", size=1))+
    ylab(subset.category)+
    scale_y_continuous(limits=c(0, NA))+ #, trans='log2', labels=comma)+
    scale_x_discrete(labels=c('healthy', 'preclinical'))
}



# Basic Components
p_per_dv_Basic <- makeDietBoxplot(diet.intake, subset=TRUE, 
                                  subset.category = 'Basic Components',
                                  wrap.var = 'Nutrients',
                                  y.val = 'PercentDV', nrows = 3)

# Vitamins
p_per_dv_Vitamins <- makeDietBoxplot(diet.intake, subset=TRUE, 
                                     subset.category = 'Vitamins',
                                     wrap.var = 'Nutrients',
                                     y.val = 'PercentDV', nrows = 4)

# Minerals
p_per_dv_Minerals <- makeDietBoxplot(diet.intake, subset=TRUE, 
                                     subset.category = 'Minerals',
                                     wrap.var = 'Nutrients',
                                     y.val = 'PercentDV', nrows = 4)

# PolyFats
p_per_dv_PolyFats <- makeDietBoxplot(diet.intake, subset=TRUE, 
                                     subset.category = 'PolyFats',
                                     wrap.var = 'Nutrients',
                                     y.val = 'PercentDV', nrows = 1)

# Other Nutrients
p_per_dv_Other <- makeDietBoxplot(diet.intake, subset=TRUE, 
                                  subset.category = 'Other Nutrients',
                                  wrap.var = 'Nutrients',
                                  y.val = 'PercentDV', nrows = 1)


# Total Calories by Source
p_total_cals <- makeDietBoxplot(total.cal, subset=FALSE, 
                                subset.category = NULL,
                                wrap.var = 'Source.of.Total.Calories',
                                y.val = 'Percent.Total.Calories', nrows = 1)


# Calories from Fats breakdown
p_fat_cals <- makeDietBoxplot(cal.fromfat, subset=FALSE, 
                              subset.category = NULL,
                              wrap.var = 'Source.of.Calories.from.Fat',
                              y.val = 'Percent.Calories.From.Fat', nrows = 1)






# RANDOM FOREST CLASSIFIER ANALYSES ############################################
# For those only interested in this analysis: 
#save(list=c('ps2.filt'), file='RData/ACS_AD/230109_STMRevision/230111_ps2filt.RData')

## Load packages and data-------------------------------------------------------
# Packages specific to these analyses: 
library(dplyr)     #v1.0.8    
library(tidyr)     #v1.2.0
library(ggplot2)   #v3.3.5 
library(phyloseq)  #v1.38.0
library(Boruta)    #v7.0.0
library(caret)     #v6.0.86
library(VIM)       #v6.1.0
library(stringr)   #v1.4.0
library(egg)       #v0.4.5
library(scales)    #v1.1.1

# The phyloseq object loaded in 230111_ps2filt.RData below is equivalent to 
# ps2.filt in the previous sections. Here it is provided separately for those 
# who are only interested in this section and would rather not recreate it. 
# Object ps2.filt contains clinical and demographic metadata for 164 healthy and 
# preclinical AD participants, as well as taxonomic relative abundance data that 
# has been filtered to exclude species with mean relative abundance < 0.1% 
# across all samples. 

load('path/to/230111_ps2filt.RData')


## Prepare data for training----------------------------------------------------

# Transform to counts, preserving precision to 0.00001%
ps2.filt.int <- transform_sample_counts(ps2.filt, function(x) trunc(x*100000))

# Access clinical metadata
meta.caret <- data.frame(sample_data(ps2.filt.int), check.names=FALSE)

# Access taxonomic abundance data (note transform of matrix)
taxa.caret <- t(data.frame(otu_table(ps2.filt.int), check.names=FALSE))


# Metadata variables to include in training
RF.vars <- c('amyloid.positive.AF',
             'combined.Centiloid', 'LUMIPULSE_CSF.ab42.ab40.ratio',
             'Tauopathy', 'LUMIPULSE_CSF_pTau', 
             'LUMIPULSE_CSF_tTau', 'WMH_volume', 'CortSig_Thickness', 'MR_LV_RV_HIPPOCAMPUS',
             'Polygenic.Risk.Score', 'APOE4.status', 
             'Age', 'Sex', 'Race', 'Years.Education', 'BMI', 'Hypertension', 'Diabetes', 
             'INT_PET', 'INT_CSF', 'INT_MRI') 

meta.caret2 <- meta.caret[, RF.vars]



# Update factor levels for ease of interpretation 
meta.caret2 <- meta.caret2 %>% 
                  mutate(APOE4.status = ifelse(APOE4.status == "e4+", "e4.pos", "e4.neg"),
                         Sex = ifelse(Sex == "1", "Male", "Female"),
                         Race = ifelse(Race == "4", "Black", 
                                       ifelse(Race == "5", "White", "Other")),
                         Hypertension = ifelse(Hypertension == "1", 'ht.yes', 'ht.no'),
                         Diabetes = ifelse(Diabetes == "1", 'dia.yes', 'dia.no'),
                         amyloid.positive.AF = ifelse(amyloid.positive.AF == "0", "healthy", "preclinical"))



## Impute missing Data----------------------------------------------------------

# Visualize missing data
par(mar = c(6,4.1,4.1,6))
meta.vimplot <- aggr(meta.caret2, numbers = TRUE, sortVars = TRUE,
                     col = c("#347B98", "#66B032", "#B2D732"),
                     labels = names(meta.caret2), cex.axis= 0.3, combined = FALSE)

# Impute using kNN algorithm. Preserve rownames (Participant ID), because this
# step will remove those. 
meta.caret2$Participant <- rownames(meta.caret2)
meta.caret.im <- kNN(meta.caret2, imp_var = FALSE)

# Replace rownames.
rownames(meta.caret.im) <- meta.caret.im$Participant
meta.caret.im$Participant <- NULL

# Update feature names for ease of interpretation.
meta.caret.im <- meta.caret.im %>% rename(PET_amyloid = combined.Centiloid, 
                                          CSF.ratio.ab42.ab40 = LUMIPULSE_CSF.ab42.ab40.ratio, 
                                          PET.Tau = Tauopathy,
                                          CSF.ptau = LUMIPULSE_CSF_pTau, 
                                          CSF.ttau = LUMIPULSE_CSF_tTau, 
                                          Cort.Signature = CortSig_Thickness,
                                          IntervalDays.PET = INT_PET, 
                                          IntervalDays.CSF = INT_CSF,
                                          IntervalDays.MRI = INT_MRI,
                                          Hippocampus.Volume = MR_LV_RV_HIPPOCAMPUS)

# VIM also removes variable class information (factors, etc.). Fix here by
# converting categorical variables to factors
factor.vars <- c('amyloid.positive.AF', 'APOE4.status', 'Sex', 'Race',
                 'Hypertension', 'Diabetes')

meta.caret.im[, factor.vars] <- lapply(meta.caret.im[, factor.vars], factor)


## Partition out training cohort for feature selection--------------------------
set.seed(42)
train_idx <- createDataPartition(meta.caret.im$amyloid.positive.AF, p = 0.6, list=FALSE)

# Access metadata (biomarkers) and taxanomic abundances for training cohort.
meta.train <- meta.caret.im[train_idx, ]
taxa.train <- taxa.caret[train_idx, ] 




## Iterative taxonomic feature selection on training cohort---------------------

# Merge sample class identity with taxonomic abundance data 
class <- subset(meta.train, select=c('amyloid.positive.AF'))
 
taxa.train.wclass <- merge(class, taxa.train, by='row.names')

rownames(taxa.train.wclass) <- taxa.train.wclass$Row.names
taxa.train.wclass$Row.names <- NULL



# FUNCTION: callBoruta helper function (called by runFeatureSelection() below)
# Required packages: Boruta, stringr
# Arguments:
#  taxa.data (dataframe)  = taxa.train.wclass (tax abundances for training cohort)
#  seed.Boruta (int)      = will be passed iteratively in defined range
#
# Return:
#  list of names of feature-selected taxa for current iteration

callBoruta <- function(taxa.data, seed.Boruta){
  set.seed(seed.Boruta)
  
  taxa.boruta <- Boruta(amyloid.positive.AF~., data = taxa.data, maxRuns=500, doTrace=0)
  taxa.boruta.fix <- TentativeRoughFix(taxa.boruta)
  
  taxa.boruta.df <- data.frame('boruta' = taxa.boruta.fix$finalDecision)
  
  taxanames <- str_replace_all(rownames(subset(taxa.boruta.df, 
                                               boruta=='Confirmed')), "`", "")
  
}


# FUNCTION: runFeatureSelection - iterative feature selection function
# Required packages: Boruta, stringr
# Arguments:
#  taxa.train.data (dataframe)   = taxa.train.wclass (passed to callBoruta())
#  seed.range (num list)         = range for iteration, e.g. 1:100.
#
# Return:
#  dataframe summarizing frequency at which unique taxa were feature-selected
#  across all iterations (random seeds).

runFeatureSelection <- function(taxa.train.data, seed.range){

  fs.taxa.train <- vector('list', length(seed.range))
  for (i in seed.range){
    fs.taxa.train[[i]] <- callBoruta(taxa.train.data, i)
  }
  
  fs.taxa.train.summ <- data.frame(table(unlist(fs.taxa.train)))
} 



# Carry out iterative feature selection (may take a little while..~20 min)
fs.taxa <- runFeatureSelection(taxa.train.wclass, 1:100)

# Recommended to save workspace at this point.

# Filter for taxa selected in > 25% of iterations.
fs.taxa.top <- subset(fs.taxa, Freq > 25, select='Var1')


# Subset taxa abundance data to these selected taxa (for entire cohort including
# train and validation sets -- will be re-partitioned later using same index). 
taxa.caret.boruta <- taxa.caret[ , colnames(taxa.caret) %in% fs.taxa.top$Var1]



## Plot relative abundances of feature-selected taxa----------------------------

# Subset melted phyloseq dataframe to the feature-selected taxa.
ps2.filt.int.boruta <- subset(ps2.filt.int.df, OTU %in% fs.taxa.top$Var1)


# Make facet labels vector 
# List of unique OTUs, split for just species name, but name the species vector 
# with full OTU name.
taxafeature.labels.df <- fs.taxa.top %>% separate(Var1, c(NA, 'species'), '\\|s__')
taxafeature.labels <- taxafeature.labels.df$species
names(taxafeature.labels) <- fs.taxa.top$Var1



# Plot relative abundances of feature-selected taxa.
# Note transformation back to percent scale, with addition of small pseudo-count
# to enable plotting on a log scale. 
p_taxafeature <- ggplot(ps2.filt.int.boruta, 
                        aes(x=amyloid.positive.AF, y=Abundance/100000+0.0001))+
  geom_violin(aes(color=amyloid.positive.AF), draw_quantiles = c(0.25,0.75))+
  geom_jitter(aes(color=amyloid.positive.AF),width=0.2, size=1, alpha=0.8, pch=21)+
  scale_color_brewer(palette ='Dark2')+
  theme_classic()+
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = "none",
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=1))+
  facet_wrap(~OTU, labeller = labeller(OTU = taxafeature.labels), nrow=1)+
  scale_y_continuous(trans="log10", labels=comma) +
  ylab('Abundance (%)')




## Create data subsets (i.e. omit biomarker categories)-------------------------

# FUNCTION: createSubsets - supply clinical metadata features to omit, then 
# merge remaining clinical metadata features with previously selected taxonomic 
# features.
# Required packages: none
# Arguments:
#  RFmetadata (dataframe)    = meta.caret.im (imputed metadata)
#  nullvars (chr list)       = list of vars to exclude from the model
#  Taxa (num matrix)         = taxa.caret.boruta (abundances of selected taxa)
#
# Return:
#  Named list: Base = values for selected metadata features, WithTaxa = values
#  for selected metadata features as well as feature-selected taxa. This will
#  enable testing of improvements in model performances with addition of
#  taxonomic feature data.

createSubsets <- function(RFmetadata, nullvars, Taxa){
  if (!all(nullvars %in% colnames(RFmetadata))){
    stop("At least one variable name not in dataframe")
  }
  
  #Non-taxonomic features only  
  metadata <- RFmetadata
  metadata[ , nullvars] <- list(NULL)
  
  #With taxonomic features
  metadata.wtax <- merge(metadata, Taxa, all = TRUE, by='row.names')
  rownames(metadata.wtax) <- metadata.wtax$Row.names
  metadata.wtax$Row.names <- NULL
  
  #RETURN
  all.out <- list("Base"= metadata, 
                  "WithTaxa"= metadata.wtax)
}



# Define metadata subsets (biomarker exclusion lists):

# Main text models: 
#  No exclusions
alpha.nullvars <- list()

#  Exclude amyloid only
beta.nullvars <- c('PET_amyloid', 'CSF.ratio.ab42.ab40')

#  Exclude all except genetics (amyloid, tau, and neurodegeneration)
omicron.nullvars <- c('PET_amyloid', 'CSF.ratio.ab42.ab40',
                      'PET.Tau', 'CSF.ptau',
                      'CSF.ttau', 'WMH_volume', 'Cort.Signature', 'Hippocampus.Volume')

#  Exclude all biomarkers
pi.nullvars <- c('PET_amyloid', 'CSF.ratio.ab42.ab40',
                 'PET.Tau', 'CSF.ptau',
                 'CSF.ttau', 'WMH_volume', 'Cort.Signature', 'Hippocampus.Volume',
                 'Polygenic.Risk.Score', 'APOE4.status')


# Supplementary models:
#  Exclude tau in addition to amyloid
zeta.nullvars <- c('PET_amyloid', 'CSF.ratio.ab42.ab40', 
                   'PET.Tau', 'CSF.ptau')

#  Exclude neurodegeneration in addition to amyloid
eta.nullvars <- c('PET_amyloid', 'CSF.ratio.ab42.ab40', 
                  'CSF.ttau', 'WMH_volume', 'Cort.Signature', 'Hippocampus.Volume')

#  Exclude genetics in addition addition to amyloid
theta.nullvars <- c('PET_amyloid', 'CSF.ratio.ab42.ab40', 
                    'Polygenic.Risk.Score', 'APOE4.status')

#  Exclude neurodegeneration and genetics in addition to amyloid
nu.nullvars <- c('PET_amyloid', 'CSF.ratio.ab42.ab40',
                 'CSF.ttau', 'WMH_volume', 'Cort.Signature', 'Hippocampus.Volume',
                 'Polygenic.Risk.Score', 'APOE4.status')

#  Exclude tau and genetics in addition to amyloid
xi.nullvars <- c('PET_amyloid', 'CSF.ratio.ab42.ab40',
                 'PET.Tau', 'CSF.ptau',
                 'Polygenic.Risk.Score', 'APOE4.status')



# Create Data Subsets
#  Main Text models:
alpha.data <- createSubsets(meta.caret.im, alpha.nullvars, Taxa=taxa.caret.boruta)
beta.data <- createSubsets(meta.caret.im, beta.nullvars, Taxa=taxa.caret.boruta)
omicron.data <- createSubsets(meta.caret.im, omicron.nullvars, Taxa=taxa.caret.boruta)
pi.data <- createSubsets(meta.caret.im, pi.nullvars, Taxa=taxa.caret.boruta)

#  Supplementary models:
theta.data <- createSubsets(meta.caret.im, theta.nullvars, Taxa=taxa.caret.boruta)
eta.data <- createSubsets(meta.caret.im, eta.nullvars, Taxa=taxa.caret.boruta)
zeta.data <- createSubsets(meta.caret.im, zeta.nullvars, Taxa=taxa.caret.boruta)
xi.data <- createSubsets(meta.caret.im, xi.nullvars, Taxa=taxa.caret.boruta)
nu.data <- createSubsets(meta.caret.im, nu.nullvars, Taxa=taxa.caret.boruta)


## Train Random Forest classifiers----------------------------------------------

# Define the test harness. Within the training cohort, will train using 10-fold
# cross-validation. 
cv10 <- trainControl(method='cv', number=10, classProbs=T, savePredictions =T)


# Define categorical variables that should not be normalized / scaled.
varsnot2norm <- c('Sex',
                  'Race',
                  'amyloid.positive.AF',  #class label
                  'APOE4.status',
                  'Diabetes',
                  'Hypertension')


# FUNCTION: train_rf_models trains a classifier on the provided training data  
# and test harness. It does this iteratively (100x) on random 80:20 partitions of
# the training cohort, testing on the entire validation cohort at each iteration. 
# Predictive results are collated.
# Required packages: caret
# Arguments:
#  data (dataframe)         = data subset generated by createSubsets()
#  control.harness          = cv10. Control harness generated  by caret::trainControl()
#  data.name.string (str)   = Identifier for the data subset provided
#  varsNot2Norm (chr list)  = List of categorical variables that shouldn't be
#                             normalized / scaled.
#  shuffle.class (logical)  = TRUE if class labels should be shuffled during
#                             model training to generate null performance 
#                             parameter distributions.
#
# Return: 
#  list(prediction results [on training cohort], prediction results [on 
#            validation cohort], variable importances)

train_rf_models <- function(data, control.harness, data.name.string, 
                            varsNot2Norm, shuffle.class) {
  #Cross Validation (within training cohort) 
  out <- list()
  varimportance <- list() 
  
  #Validation Set 
  out.val <- list()
  
  #Separate into train/test and validation subsets, using same index as before.
  set.seed(42)
  data1_idx <- createDataPartition(data$amyloid.positive.AF, p = 0.6, list=FALSE)
  data1 <- data[data1_idx, ]
  data.val <- data[-data1_idx, ]
  
  for (i in 1:100) {
    #Create random partition of training cohort (80:20). 
    set.seed(i)
    
    train_idx <- createDataPartition(data1$amyloid.positive.AF, p = 0.8, list=FALSE)
    
    data.train <- data1[train_idx, ]
    data.test <- data1[-train_idx, ]
    
    # Optionally shuffle class labels in training data (data.train).
    if (shuffle.class == TRUE) {
      data.train$amyloid.positive.AF <- sample(data.train$amyloid.positive.AF)
    }

    #Pre-process data (center and scale)
    preprocessrule <- preProcess(data.train[, !(colnames(data.train) %in% varsNot2Norm)], 
                                 method = c('center', 'scale'))
    data.train.p <- predict(preprocessrule, data.train)
    data.test.p <- predict(preprocessrule, data.test)
    
    data.val.p <- predict(preprocessrule, data.val)
    
    
    #train model
    set.seed(42)
    fit.rf <- train(amyloid.positive.AF~., data = data.train.p, method = 'rf',
                    metric = 'Accuracy', trControl = control.harness)
    
    #Make predictions for this iteration's test set and store performance measures.
    predictions.rf <- predict(fit.rf, data.test.p)
    cM <- confusionMatrix(predictions.rf, data.test.p$amyloid.positive.AF)
    pred.results <- c('healthy-healthy'=cM$table[1,1], 
                      'healthy-preclinical'=cM$table[2,1], 
                      'preclinical-healthy'=cM$table[1,2], 
                      'preclinical-preclinical'= cM$table[2,2], 
                      cM$overall, 
                      cM$byClass,
                      'Data' = data.name.string, 
                      'Seed' = i)
    
    out[[i]] <- pred.results
    
    #Make predictions for retained VALIDATION set and store performance measures
    predictions.val.rf <- predict(fit.rf, data.val.p)
    cM.val <- confusionMatrix(predictions.val.rf, data.val.p$amyloid.positive.AF)
    pred.results.val <- c('healthy-healthy'=cM.val$table[1,1], 
                          'healthy-preclinical'=cM.val$table[2,1], 
                          'preclinical-healthy'=cM.val$table[1,2], 
                          'preclinical-preclinical'= cM.val$table[2,2], 
                          cM.val$overall, 
                          cM.val$byClass,
                          'Data' = data.name.string, 
                          'Seed' = i)
    
    out.val[[i]] <- pred.results.val
    
    
    #Find important vars
    fit.rf.importance <- varImp(fit.rf, scale=FALSE)
    fit.rf.importance.df <- fit.rf.importance$importance
    var.importance <- fit.rf.importance.df$Overall
    names(var.importance) <- row.names(fit.rf.importance.df)
    
    varimportance[[i]] <- var.importance
  }
  
  #Return 
  out.df <- data.frame(do.call('rbind', out))
  out.val.df <- data.frame(do.call('rbind', out.val))
  varimportance.df <- data.frame(do.call('rbind', varimportance))
  
  allout <- list('Pred.Results.CV'=out.df, 
                 'Pred.Results.Val'=out.val.df, 
                 'Var.Importance'=varimportance.df)
}



# Train the models (without shuffling class labels)

# Alpha (all biomarkers)
rf.alpha <- train_rf_models(alpha.data$Base, cv10, 'Alpha', 
                            varsnot2norm, FALSE)
rf.alpha_T <- train_rf_models(alpha.data$WithTaxa, cv10, 'Alpha plus Taxa', 
                              varsnot2norm, FALSE)


# Beta (Excluding amyloid)
rf.beta <- train_rf_models(beta.data$Base, cv10, 'Beta', 
                           varsnot2norm, FALSE)
rf.beta_T <- train_rf_models(beta.data$WithTaxa, cv10, 'Beta plus Taxa', 
                             varsnot2norm, FALSE)


# Omicron (Excluding all biomarkers except genetics) 
rf.omicron <- train_rf_models(omicron.data$Base, cv10, 'Omicron', 
                              varsnot2norm, FALSE)
rf.omicron_T <- train_rf_models(omicron.data$WithTaxa, cv10, 'Omicron plus Taxa', 
                                varsnot2norm, FALSE)


# Pi (Excluding all biomarkers) 
rf.pi <- train_rf_models(pi.data$Base, cv10, 'Pi', 
                         varsnot2norm, FALSE)
rf.pi_T <- train_rf_models(pi.data$WithTaxa, cv10, 'Pi plus Taxa', 
                           varsnot2norm, FALSE)


# Zeta (No amyloid or tau)
rf.zeta <- train_rf_models(zeta.data$Base, cv10, 'Zeta', 
                           varsnot2norm, FALSE)
rf.zeta_T <- train_rf_models(zeta.data$WithTaxa, cv10, 'Zeta plus Taxa', 
                             varsnot2norm, FALSE)


# Eta (No amyloid or neurodegeneration)
rf.eta <- train_rf_models(eta.data$Base, cv10, 'Eta', 
                          varsnot2norm, FALSE)
rf.eta_T <- train_rf_models(eta.data$WithTaxa, cv10, 'Eta plus Taxa',
                            varsnot2norm, FALSE)

# Theta (No amyloid or genetics) 
rf.theta <- train_rf_models(theta.data$Base, cv10, 'Theta', 
                            varsnot2norm, FALSE)
rf.theta_T <- train_rf_models(theta.data$WithTaxa, cv10, 'Theta plus Taxa', 
                              varsnot2norm, FALSE)


#Nu (only including tau) 
rf.nu <- train_rf_models(nu.data$Base, cv10, 'Nu', 
                         varsnot2norm, FALSE)
rf.nu_T <- train_rf_models(nu.data$WithTaxa, cv10, 'Nu plus Taxa', 
                           varsnot2norm, FALSE)

#Xi (only including neurodegeneration) 
rf.xi <- train_rf_models(xi.data$Base, cv10, 'Xi', 
                         varsnot2norm, FALSE)
rf.xi_T <- train_rf_models(xi.data$WithTaxa, cv10, 'Xi plus Taxa', 
                           varsnot2norm, FALSE)



# Collate results: 
all.pred.CV <- rbind(rf.alpha$Pred.Results.CV, rf.alpha_T$Pred.Results.CV,
                     rf.beta$Pred.Results.CV, rf.beta_T$Pred.Results.CV,
                     rf.omicron$Pred.Results.CV, rf.omicron_T$Pred.Results.CV,
                     rf.pi$Pred.Results.CV, rf.pi_T$Pred.Results.CV,
                     rf.zeta$Pred.Results.CV, rf.zeta_T$Pred.Results.CV,
                     rf.eta$Pred.Results.CV, rf.eta_T$Pred.Results.CV,
                     rf.theta$Pred.Results.CV, rf.theta_T$Pred.Results.CV,
                     rf.nu$Pred.Results.CV, rf.nu_T$Pred.Results.CV,
                     rf.xi$Pred.Results.CV, rf.xi_T$Pred.Results.CV)



all.pred.Val <- rbind(rf.alpha$Pred.Results.Val, rf.alpha_T$Pred.Results.Val,
                      rf.beta$Pred.Results.Val, rf.beta_T$Pred.Results.Val,
                      rf.omicron$Pred.Results.Val, rf.omicron_T$Pred.Results.Val,
                      rf.pi$Pred.Results.Val, rf.pi_T$Pred.Results.Val,
                      rf.zeta$Pred.Results.Val, rf.zeta_T$Pred.Results.Val,
                      rf.eta$Pred.Results.Val, rf.eta_T$Pred.Results.Val,
                      rf.theta$Pred.Results.Val, rf.theta_T$Pred.Results.Val,
                      rf.nu$Pred.Results.Val, rf.nu_T$Pred.Results.Val,
                      rf.xi$Pred.Results.Val, rf.xi_T$Pred.Results.Val)




all.VarImp.list <- list(rf.alpha$Var.Importance, rf.alpha_T$Var.Importance,
                        rf.beta$Var.Importance, rf.beta_T$Var.Importance,
                        rf.omicron$Var.Importance, rf.omicron_T$Var.Importance,
                        rf.pi$Var.Importance, rf.pi_T$Var.Importance,
                        rf.zeta$Var.Importance, rf.zeta_T$Var.Importance,
                        rf.eta$Var.Importance, rf.eta_T$Var.Importance,
                        rf.theta$Var.Importance, rf.theta_T$Var.Importance,
                        rf.nu$Var.Importance, rf.nu_T$Var.Importance,
                        rf.xi$Var.Importance, rf.xi_T$Var.Importance)

all.VarImp <- all.VarImp.list %>% Reduce(function(d1, d2) full_join(d1, d2), .)



## Shuffle class labels and retrain for null distributions----------------------

# We'd like to compare importance of individual features to the performance of 
# the predictive models against their importance when class labels have been
# shuffled during model training (to generate random null distributions for 
# the variable importances). 

# Alpha (all biomarkers)
rfN.alpha <- train_rf_models(alpha.data$Base, cv10, 'Alpha', 
                            varsnot2norm, TRUE)
rfN.alpha_T <- train_rf_models(alpha.data$WithTaxa, cv10, 'Alpha plus Taxa', 
                              varsnot2norm, TRUE)


# Beta (Excluding amyloid)
rfN.beta <- train_rf_models(beta.data$Base, cv10, 'Beta', 
                           varsnot2norm, TRUE)
rfN.beta_T <- train_rf_models(beta.data$WithTaxa, cv10, 'Beta plus Taxa', 
                             varsnot2norm, TRUE)


# Omicron (Excluding all biomarkers except genetics) 
rfN.omicron <- train_rf_models(omicron.data$Base, cv10, 'Omicron', 
                              varsnot2norm, TRUE)
rfN.omicron_T <- train_rf_models(omicron.data$WithTaxa, cv10, 'Omicron plus Taxa', 
                                varsnot2norm, TRUE)


# Pi (Excluding all biomarkers) 
rfN.pi <- train_rf_models(pi.data$Base, cv10, 'Pi', 
                         varsnot2norm, TRUE)
rfN.pi_T <- train_rf_models(pi.data$WithTaxa, cv10, 'Pi plus Taxa', 
                           varsnot2norm, TRUE)


# Zeta (No amyloid or tau)
rfN.zeta <- train_rf_models(zeta.data$Base, cv10, 'Zeta', 
                           varsnot2norm, TRUE)
rfN.zeta_T <- train_rf_models(zeta.data$WithTaxa, cv10, 'Zeta plus Taxa', 
                             varsnot2norm, TRUE)


# Eta (No amyloid or neurodegeneration)
rfN.eta <- train_rf_models(eta.data$Base, cv10, 'Eta', 
                          varsnot2norm, TRUE)
rfN.eta_T <- train_rf_models(eta.data$WithTaxa, cv10, 'Eta plus Taxa',
                            varsnot2norm, TRUE)

# Theta (No amyloid or genetics) 
rfN.theta <- train_rf_models(theta.data$Base, cv10, 'Theta', 
                            varsnot2norm, TRUE)
rfN.theta_T <- train_rf_models(theta.data$WithTaxa, cv10, 'Theta plus Taxa', 
                              varsnot2norm, TRUE)


#Nu (only including tau) 
rfN.nu <- train_rf_models(nu.data$Base, cv10, 'Nu', 
                         varsnot2norm, TRUE)
rfN.nu_T <- train_rf_models(nu.data$WithTaxa, cv10, 'Nu plus Taxa', 
                           varsnot2norm, TRUE)

#Xi (only including neurodegeneration) 
rfN.xi <- train_rf_models(xi.data$Base, cv10, 'Xi', 
                         varsnot2norm, TRUE)
rfN.xi_T <- train_rf_models(xi.data$WithTaxa, cv10, 'Xi plus Taxa', 
                           varsnot2norm, TRUE)



# Collate results: 
NULL.all.pred.CV <- rbind(rfN.alpha$Pred.Results.CV, rfN.alpha_T$Pred.Results.CV,
                     rfN.beta$Pred.Results.CV, rfN.beta_T$Pred.Results.CV,
                     rfN.omicron$Pred.Results.CV, rfN.omicron_T$Pred.Results.CV,
                     rfN.pi$Pred.Results.CV, rfN.pi_T$Pred.Results.CV,
                     rfN.zeta$Pred.Results.CV, rfN.zeta_T$Pred.Results.CV,
                     rfN.eta$Pred.Results.CV, rfN.eta_T$Pred.Results.CV,
                     rfN.theta$Pred.Results.CV, rfN.theta_T$Pred.Results.CV,
                     rfN.nu$Pred.Results.CV, rfN.nu_T$Pred.Results.CV,
                     rfN.xi$Pred.Results.CV, rfN.xi_T$Pred.Results.CV)



NULL.all.pred.Val <- rbind(rfN.alpha$Pred.Results.Val, rfN.alpha_T$Pred.Results.Val,
                      rfN.beta$Pred.Results.Val, rfN.beta_T$Pred.Results.Val,
                      rfN.omicron$Pred.Results.Val, rfN.omicron_T$Pred.Results.Val,
                      rfN.pi$Pred.Results.Val, rfN.pi_T$Pred.Results.Val,
                      rfN.zeta$Pred.Results.Val, rfN.zeta_T$Pred.Results.Val,
                      rfN.eta$Pred.Results.Val, rfN.eta_T$Pred.Results.Val,
                      rfN.theta$Pred.Results.Val, rfN.theta_T$Pred.Results.Val,
                      rfN.nu$Pred.Results.Val, rfN.nu_T$Pred.Results.Val,
                      rfN.xi$Pred.Results.Val, rfN.xi_T$Pred.Results.Val)




NULL.all.VarImp.list <- list(rfN.alpha$Var.Importance, rfN.alpha_T$Var.Importance,
                        rfN.beta$Var.Importance, rfN.beta_T$Var.Importance,
                        rfN.omicron$Var.Importance, rfN.omicron_T$Var.Importance,
                        rfN.pi$Var.Importance, rfN.pi_T$Var.Importance,
                        rfN.zeta$Var.Importance, rfN.zeta_T$Var.Importance,
                        rfN.eta$Var.Importance, rfN.eta_T$Var.Importance,
                        rfN.theta$Var.Importance, rfN.theta_T$Var.Importance,
                        rfN.nu$Var.Importance, rfN.nu_T$Var.Importance,
                        rfN.xi$Var.Importance, rfN.xi_T$Var.Importance)

NULL.all.VarImp <- NULL.all.VarImp.list %>% Reduce(function(d1, d2) full_join(d1, d2), .)




## Plot predictive performance metrics------------------------------------------ 

# We'd like to plot the distributions of accuracy, sensitivity, and specificity
# of predictions on the validation cohort for each of the trained models. Also
# plot null distributions of these parameters from the class label shuffling. 

# For ease of faceting in plots, add columns DataSet indicating the base model
# (e.g. Alpha, Beta ..), and Taxa indicating if feature-selected taxa were 
# included as features in the model, to prediction results (all.pred.Val). 

all.pred.Val2 <- all.pred.Val %>% 
  separate(Data, into=c('DataSet', 'Taxa'), sep = ' plus ', 
           remove = FALSE, fill = 'right') %>% 
  mutate(Taxa = ifelse(is.na(Taxa), 'No Microbiome Data', 
                      'Including Selected Taxa'))


# FUNCTION: performancePlots will plot performance metrics for supplied models, 
# comparing models trained on data including or omitting taxonomic features.
# Required packages: ggplot2, tidyr, dplyr, forcats
# Arguments:
#  all.pred.df2 (dataframe)    = Combined output of train_rf_models() 
#                                (all.pred.CV2 or all.pred.Val2) after 
#                                execution on multiple data sets. 
#                                Should include columns 'DataSet' and 'Taxa'.  
#  mod.order (chr list)        = Order in which to plot models.
#  include.mod (chr list)      = Which models to plot. 
#  tax.order (chr list)        = Order in which to plot +/- taxonomic feature
#                                results.
#  pred.vars (chr list)        = Performance metrics to plot.
#  colors.fill (chr list)      = Color palette for box plot fill. Length should
#                                equal length(tax.order).
#  colors.col (chr list)       = Color palette for boxplot outline. Length 
#                                should equal length(tax.order).
# 
# Return: 
#  Plot          = boxplots of predictive performance metrics.
#  all.pred.long = long format data frame of performance metrics, REQUIRED FOR
#                  processAOVS() BELOW. 


performancePlots <- function(all.pred.df2, mod.order, include.mod, tax.order, 
                             pred.vars, colors.fill, colors.col){
  
  # Re-order DataSet (model) and Taxa features given provided arguments
  all.pred <- all.pred.df2 %>% mutate(DataSet = fct_relevel(DataSet, mod.order),
                                      Taxa = fct_relevel(Taxa, tax.order))
  
  # Gather performance metrics of interest and fix variable types.
  all.pred.long <- gather_(all.pred, 'Performance_Measure', 
                           'Performance_Measure_Value', pred.vars)

  all.pred.long$Performance_Measure <- as.factor(all.pred.long$Performance_Measure)
  all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)
    
  # Subset to models we want to include
  all.pred.long <- subset(all.pred.long, DataSet %in% include.mod)

  
  # Plotting function call
  p_rf_pred <- ggplot(all.pred.long, aes(x=DataSet, y=Performance_Measure_Value, 
                                         fill=Taxa, color=Taxa))+
    geom_boxplot(position=position_dodge(width=0.5),alpha=0.5, outlier.size=0.5)+
    stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5))+
    facet_wrap(~Performance_Measure, nrow =1)+
    scale_fill_manual(values = colors.fill)+
    scale_color_manual(values = colors.col)+
    theme_classic()+
    theme(axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", size=0.5))+
    scale_x_discrete(limits = rev)+
    ylim(0,NA)+   
    coord_flip()
  
  # RETURN PLOT OBJECT AND ALL.PRED.LONG
  all.out <- list('Plot' = p_rf_pred, 'all.pred.long' = all.pred.long)
}


# Plot model performances for main text and supplement: 
taxincl.order <- c('No Microbiome Data', 'Including Selected Taxa')
pred_vars <- c('Accuracy', 'Sensitivity', 'Specificity')
fillcols <- c('#476b6b', '#8de08d') 
linecols <- c('#1f2e2e', '#004d00')
model.order <- c('Alpha', 'Beta', 'Theta', 'Eta', 'Zeta','Omicron','Xi','Nu','Pi')
main.models <- c('Alpha', 'Beta', 'Omicron', 'Pi')
suppl.models <- c('Theta', 'Eta', 'Zeta', 'Xi', 'Nu')


pred_Val.Main <- performancePlots(all.pred.Val2, model.order, main.models, 
                                  taxincl.order, pred_vars, fillcols, linecols)

pred_Val.Supp <- performancePlots(all.pred.Val2, model.order, suppl.models, 
                                  taxincl.order, pred_vars, fillcols, linecols)




## Compare predictive performance of models (ANOVA)-----------------------------

# FUNCTION: modelAOVs compares predictive performance metrics of RF models, in this
# case models including or omitting taxonomic features. Can also handle comparisons
# between more than 2 groups (one-way ANOVA). If multiple models (trained on 
# different base datasets, e.g. Alpha, Beta, etc.) are supplied, will collate 
# results. Tukey's HSD is applied to each one-way ANOVA.
# Required packages: none other than default stats package
# Arguments:
#  model.order (chr list)          = list of models to test
#  predvars (chr list)             = list of performance metrics to compare
#  RFpredresults (dataframe)       = Combined output of train_rf_models() 
#                                    (all.pred.CV2 or all.pred.Val2) after 
#                                    execution on multiple data sets.  
# Return:
#  AOVresults (dataframe)     = Results of one-way ANOVA on performance metrics 
#                               collated across the models provided.
#  Tukeyresults (dataframe)   = TukeyHSD applied to each ANOVA.                         

modelAOVs <- function(model.order, predvars, RFpredresults){
  aov.out <- list()
  tk.out <- list()
  i = 0
  k = 0
  for (model in model.order) {
    model.data <- subset(RFpredresults, DataSet==model)
    for (var in predvars) {
      i = i + 1
      model.var.aov <- aov(model.data[[as.name(var)]] ~ model.data$Taxa)
      
      model.summ <- unlist(summary(model.var.aov))
      model.results <- c("Model" = model, 
                         "PerformanceParameter" = var, 
                         "AOV.pval" = model.summ['Pr(>F)1'])
      aov.out[[i]] <- model.results
      
      model.var.tk <- TukeyHSD(model.var.aov)
      print(model.var.tk)
      model.var.tk.df <- data.frame(model.var.tk$'model.data$Taxa')
      for (comparison in 1:dim(model.var.tk.df)[1]){
        k = k + 1
        tk.results <- c("Model" = model, "PerformanceParameter" = var,
                        "Comparison" = rownames(model.var.tk.df)[comparison],
                        "Diff" = model.var.tk.df$diff[comparison],
                        "CI95lower" = model.var.tk.df$lwr[comparison],
                        "CI95upper" = model.var.tk.df$upr[comparison],
                        "Tk.pval" = model.var.tk.df$p.adj[comparison])
        tk.out[[k]] <- tk.results
      }
      
    }
  }
  aov.out.df <- data.frame(do.call('rbind', aov.out))
  tk.out.df <- data.frame(do.call('rbind', tk.out))
  
  allout <- list('AOVResults'=aov.out.df, 'TukeyResults'=tk.out.df)
}


# Relevel 'Taxa' to set reference for contrasts. 
all.pred.Val2$Taxa <- relevel(as.factor(all.pred.Val2$Taxa), ref='No Microbiome Data')

# Execute ANOVA of predictive performance metrics with Tukey's post hoc test. 
ModelAOVs.Val <- modelAOVs(model.order = model.order, predvars=pred_vars, 
                                  RFpredresults=all.pred.Val2)


# Downstream processing of ANOVA results:

# FUNCTION: processAOVs takes the results of modelAOVs and applies correction
# for multiple hypothesis testing, subsetting to significant results only, and
# provides summary output tables for reporting.
# Required packages: dplyr, tidyr, forcats
# Arguments: 
#  modelAOV.object (list)       = output of modelAOVs()
#  padjust (str)                = p.adjust.method to apply (e.g.'bonferroni')
#  all.pred.df2 (dataframe)     = Output of train_models_rf(), collated after
#                                 executing on more than one dataset. 
#  mod.order (chr list)         = Models to include
#  tax.order (chr list)         = Order in which to include models trained 
#                                 with or without taxonomic freatures.
#  pred.vars (chr list)         = Performance parameters to compare.
# Return:
#  SigAOV (dataframe)      = summary of cases in which inclusion of taxonomic
#                            features significantly changed predictive performance
#                            of a model, adjusted for multiple hypothesis testing.
#  SigTK (dataframe)       = TukeyHSD post hoc comparisons for the tests summarized
#                            in SigAOV, with further p-value adjustment, and 
#                            parameter means and sd for the +/- taxonomic feature
#                            groups.

processAOVs <- function(modelAOV.object, padjust, all.pred.df2, 
                        mod.order, tax.order, pred.vars){
  # ACCESS ANOVA AND TUKEYHSD P-VALUES AND CONVERT TO NUMERIC.
  aov.results <- modelAOV.object$AOVResults
  tk.results <- modelAOV.object$TukeyResults
  
  colnames(aov.results)[3] <- 'AOV.pval'
  
  aov.results$AOV.pval <- as.numeric(aov.results$AOV.pval)
  tk.results$Tk.pval <- as.numeric(tk.results$Tk.pval)
  
  # CREATE UNIQUE MODEL/PERFORMANCE PARAMETER (MP) VARIABLES
  aov.results$MP <- paste0(aov.results$Model, '_', aov.results$PerformanceParameter)
  tk.results$MP <- paste0(tk.results$Model, '_', tk.results$PerformanceParameter)
  
  # ANOVA PVAL ADJUST
  aov.results$AOV.pval.adj <- p.adjust(aov.results$AOV.pval, method = padjust)
  aov.results.sig <- subset(aov.results, AOV.pval.adj <= 0.05)
  
  # SELECT (SIGNIFICANT) TUKEY RESULTS CORRESPONDING TO SIGNIFICANT ANOVA RESULTS
  tk.results.sig <- subset(tk.results, tk.results$MP %in% aov.results.sig$MP)
  tk.results.sig2 <- subset(tk.results.sig, Tk.pval <= 0.05)
  
  # SELECT SIGNIFICANT COMPARISONS OF INTEREST, WHERE ADDITION OF TAXONOMIC
  # TAXONOMIC FEATURES IMPROVED PERFORMANCE MEASURES. 
  tk.results.sig2 <- tk.results.sig2 %>% separate(Comparison, into=c('A','B'), 
                                                  sep='-', remove=FALSE )
  tk.results.sig2 <- subset(tk.results.sig2, !(B=='No Microbiome Data' & Diff < 0))
  
  # TUKEY PVAL ADJUST (AND ADD STARS)
  # Stars for unadjusted Tukey P-values
  tk.results.sig2 <- tk.results.sig2 %>% 
    mutate(Stars = if_else(Tk.pval <= 0.05 & Tk.pval > 0.01, '*',
                           if_else(Tk.pval <= 0.01 & Tk.pval > 0.001, '**', 
                                   if_else(Tk.pval <= 0.001,'***', 'NS'))))
  # Stars for adjusted Tukey P-values
  tk.results.sig2$Tk.pval.adj <- p.adjust(tk.results.sig2$Tk.pval, method = padjust)
  tk.results.sig2 <- tk.results.sig2 %>% 
    mutate(Adj_Stars = if_else(Tk.pval.adj <= 0.05 & Tk.pval.adj > 0.01, '*',
                               if_else(Tk.pval.adj <= 0.01 & Tk.pval.adj > 0.001, '**', 
                                       if_else(Tk.pval.adj <= 0.001, '***', 'NS'))))
  
  
  
  # FOR REPORTING - ADD PARAMETER GROUP MEANS FOR SIGNIFICANT AOV COMPARISONS
  # Re-order DataSet (model) and Taxa features given provided arguments
  all.pred <- all.pred.df2 %>% mutate(DataSet = fct_relevel(DataSet, mod.order),
                                      Taxa = fct_relevel(Taxa, tax.order))
  
  # Gather performance metrics of interest and fix variable types.
  all.pred.long <- gather_(all.pred, 'Performance_Measure', 
                           'Performance_Measure_Value', pred.vars)
  
  all.pred.long$Performance_Measure <- as.factor(all.pred.long$Performance_Measure)
  all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)
  
  
  
  # SUMMARISE PERFORMANCE METRICS (MEAN, SD), BY MODEL AND WHETHER TAXA WERE
  # INCLUDED AS FEATURES.
  all.pred.summ <- all.pred.long %>% group_by(DataSet, Taxa, Performance_Measure) %>% 
    summarise(mean = mean(Performance_Measure_Value), sd = sd(Performance_Measure_Value))
  
  
  # SUBSET TO AND MERGE WITH SIGNIFICANT TUKEY RESULTS.
  tk.results.sig3 <- merge(tk.results.sig2, all.pred.summ, 
                           by.x=c('Model', 'PerformanceParameter', 'A'),
                           by.y=c('DataSet', 'Performance_Measure', 'Taxa'))
  tk.results.sig3 <- tk.results.sig3 %>% rename(mean_A = mean, sd_A = sd)
  
  tk.results.sig3 <- merge(tk.results.sig3, all.pred.summ, 
                           by.x=c('Model', 'PerformanceParameter', 'B'),
                           by.y=c('DataSet', 'Performance_Measure', 'Taxa'))
  tk.results.sig3 <- tk.results.sig3 %>% rename(mean_B = mean, sd_B = sd)
  
  
  tk.results.sig3.wide <- spread(tk.results.sig3, PerformanceParameter, Diff)
  
  
  # RETURN
  all.out <- list('SigAOV' = aov.results.sig, 'SigTK' = tk.results.sig3.wide)
}


model.order <- c('Alpha', 'Beta', 'Theta', 'Eta', 'Zeta','Omicron','Xi','Nu','Pi')
taxa.order <- c('No Microbiome Data', 'Including Selected Taxa')
pred_vars <- c('Accuracy', 'Sensitivity', 'Specificity')

processedAOVS.Val <- processAOVs(ModelAOVs.Val, 'bonferroni', all.pred.Val2, 
                               model.order, taxa.order, pred_vars)



## Plot empirical and null variable importances---------------------------------

# Add DataSet and Taxa variables to the VARIMP dataframes. These are inherited
# from all.pred.Val2 since the data were generated simultaneously in 
# train_rf_models(). 
all.VarImp2 <- all.VarImp
all.VarImp2$DataSet <- all.pred.Val2$DataSet
all.VarImp2$Taxa <- all.pred.Val2$Taxa
all.VarImp2$X <- NULL

all.VarImp_NULL2 <- NULL.all.VarImp
all.VarImp_NULL2$DataSet <- all.pred.Val2$DataSet
all.VarImp_NULL2$Taxa <- all.pred.Val2$Taxa
all.VarImp_NULL2$X<- NULL



# Summarise Variable Importances by Model (DataSet), +/- GM features (Taxa)
all.VarImp.mean <- all.VarImp2 %>% group_by(DataSet, Taxa) %>%
  summarise_all(mean) 

all.VarImp.sd <- all.VarImp2 %>% group_by(DataSet, Taxa) %>%
  summarise_all(sd) 

all.VarImpNULL.mean <- all.VarImp_NULL2 %>% group_by(DataSet, Taxa) %>%
  summarise_all(mean)

all.VarImpNULL.sd <- all.VarImp_NULL2 %>% group_by(DataSet, Taxa) %>%
  summarise_all(sd)



# Convert to long and merge into a single data frame
VarImp.mean.long <- gather(all.VarImp.mean, Feature, Feature.MeanImp, 
                                  PET_amyloid:X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruminococcus.s__Ruminococcus_lactaris.,
                                  factor_key=TRUE)

VarImp.sd.long <- gather(all.VarImp.sd, Feature, Feature.SDImp, 
                                PET_amyloid:X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruminococcus.s__Ruminococcus_lactaris.,
                                factor_key=TRUE)

VarImp.NULL.mean.long <- gather(all.VarImpNULL.mean, Feature, Feature.NULLMeanImp, 
                                PET_amyloid:X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruminococcus.s__Ruminococcus_lactaris.,
                                factor_key=TRUE)

VarImp.NULL.sd.long <- gather(all.VarImpNULL.sd, Feature, Feature.NULLSDImp, 
                         PET_amyloid:X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruminococcus.s__Ruminococcus_lactaris.,
                         factor_key=TRUE)


VarImp.summ <- merge(VarImp.mean.long, VarImp.sd.long)
VarImp.summ <- merge(VarImp.summ, VarImp.NULL.mean.long)
VarImp.summ <- merge(VarImp.summ, VarImp.NULL.sd.long)



# We are only interested in plotting variable importance when taxonomic 
# features were included in the models. 
# Subset to omit the 'No Microbiome Data' models.
VarImp.summ2 <- subset(VarImp.summ, Taxa != 'No Microbiome Data')

# Feature Race, level 'Other' has near 0 importance in all models (only one
# individual had Race - 'Other' in the whole cohort).
VarImp.summ2 <- subset(VarImp.summ2, Feature != 'RaceOther')


# Split out to main text models and supplementary models
VarImp.summ.main <- subset(VarImp.summ2, DataSet %in% main.models)
VarImp.summ.suppl <- subset(VarImp.summ2, DataSet %in% suppl.models)


# Fix feature names for plotting
varimp.labels <- c("Age" = "Age",                                                                                                                                          
                   "APOE4.statuse4.pos" = "APOE ??4",                                                                                                                                  
                   "BMI" = "BMI",                                                                                                                                          
                   "Cort.Signature" = "Cortical Signature",                                                                                                                               
                   "CSF.ptau" = "CSF p-tau-181", 
                   "CSF.ratio.ab42.ab40" = "CSF A??42/A??40",                                                                                                                          
                   "CSF.ttau" = "CSF t-tau",                                                                                                                                     
                   "Diabetesdia.yes" = "Diabetes",                                                                                                                              
                   "Years.Education" = "Education (Years)",                                                                                                                                         
                   "Hippocampus.Volume" = "Hippocampus Volume",                                                                                                                           
                   "Hypertensionht.yes" = "Hypertension",                                                                                                                           
                   "IntervalDays.CSF" = "CSF - Stool Interval (Days)",                                                                                                                             
                   "IntervalDays.MRI" = "MRI - Stool Interval (Days)",                                                                                                                             
                   "IntervalDays.PET" = "PET - Stool Interval (Days)",                                                                                                                             
                   "PET.Tau" = "PET tau",                                                                                                                                      
                   "PET_amyloid" = "PET amyloid",                                                                                                                                  
                   "Polygenic.Risk.Score" = "Polygenic Risk Score",                                                                                                                             
                   #   "race2Other" = "Race - Other",                                                                                                                                 
                   "RaceWhite" = "Race",                                                                                                                                   
                   "SexMale" = "Sex",                                                                                                                                      
                   "WMH_volume" = "WMH volume",                                                                                                                                   
                   "X.k__Archaea.p__Euryarchaeota.c__Methanobacteria.o__Methanobacteriales.f__Methanobacteriaceae.g__Methanosphaera.s__Methanosphaera_stadtmanae." = "Methanosphaera_stadtmanae",
                   "X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Lachnospiraceae.g__Anaerostipes.s__Anaerostipes_hadrus." = "Anaerostipes_hadrus",                       
                   "X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Lachnospiraceae.g__Coprococcus.s__Coprococcus_catus." = "Coprococcus_catus",                         
                   "X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Lachnospiraceae.g__Dorea.s__Dorea_formicigenerans." = "Dorea_formicigenerans",                            
                   "X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Oscillospiraceae.g__Oscillibacter.s__Oscillibacter_sp_57_20." = "Oscillibacter_sp_57_20",                 
                   "X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii." = "Faecalibacterium_prausnitzii",           
                   "X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruminococcus.s__Ruminococcus_lactaris." = "Ruminococcus_lactaris")




# FUNCTION plotVarImp takes long-format summary of variable importances (mean
# and SD) for a given model, both the empiricial variable importances and
# those from training after class label shuffling (to find null distributions of
# variable importances). It plots variable importances (bar plots with SD error 
# bars) ranked from most to least important, and overlays null mean and SD variable
# importances from the permutational analysis (class label shuffling).
# Required packages: ggplot2
# Arguments:
#  VarImp.df (dataframe)       = Mean and SD of variable importances for a model,
#                                and optionally multiple models. Should include
#                                column names: 'Feature.MeanImp', 'Feature.SDImp',
#                                'Feature.NULLMeanImp' and 'Feature.NULLSDImp',  
#                                as well as 'DataSet' (model names).
#  Data.Set (str)              = Model (e.g. 'Alpha') to which to subset, if 
#                                data collated for more than one model are provided.
#  var.labels (named chr list) = Feature labels for plotting.
# Return:
#  Plot of ranked variable importances for a given model, both empirical and 
#  null from class label shuffling. 

plotVarImp <- function(VarImp.df, Data.Set, var.labels){
  
  p_varimp <- ggplot(subset(VarImp.df, DataSet== Data.Set & !is.na(Feature.MeanImp)), 
                                aes(x=reorder(Feature, Feature.MeanImp), y=Feature.MeanImp))+
    geom_col(fill=fillcols[2])+
    geom_point(color='black')+
    geom_errorbar(aes(ymin = Feature.MeanImp - Feature.SDImp, 
                      ymax = Feature.MeanImp + Feature.SDImp), width=0.7, size=0.8)+
    geom_point(aes(x=reorder(Feature, Feature.MeanImp), y=Feature.NULLMeanImp), 
               color='cyan', alpha=0.80)+
    geom_errorbar(aes(ymin = Feature.NULLMeanImp - Feature.NULLSDImp, 
                      ymax = Feature.NULLMeanImp + Feature.NULLSDImp), width=0.5, size=0.6, 
                  color='cyan', alpha=0.80)+
    theme_classic()+
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y = element_blank(),
          plot.title = element_blank(),
          panel.border = element_rect(fill=NA, colour="black", size=0.5))+
    facet_wrap(~DataSet, nrow=2)+
    scale_x_discrete(labels = var.labels)+
    coord_flip()
}


# Main text models
pvar_alpha <- plotVarImp(VarImp.summ.main, 'Alpha', varimp.labels)
pvar_beta <- plotVarImp(VarImp.summ.main, 'Beta', varimp.labels)
pvar_omicron  <- plotVarImp(VarImp.summ.main, 'Omicron', varimp.labels)
pvar_pi <- plotVarImp(VarImp.summ.main, 'Pi', varimp.labels)

pvarimp_main <- ggarrange(pvar_alpha, pvar_beta,
                           pvar_omicron, pvar_pi, nrow=2)



# Supplementary models
pvar_theta <- plotVarImp(VarImp.summ.suppl, 'Theta', varimp.labels)
pvar_eta <- plotVarImp(VarImp.summ.suppl, 'Eta', varimp.labels)
pvar_zeta  <- plotVarImp(VarImp.summ.suppl, 'Zeta', varimp.labels)
pvar_xi <- plotVarImp(VarImp.summ.suppl, 'Xi', varimp.labels)
pvar_nu <- plotVarImp(VarImp.summ.suppl, 'Nu', varimp.labels)

pvarimp_suppl <- ggarrange(pvar_theta, pvar_eta, pvar_zeta,
                          pvar_xi, pvar_nu, nrow=3)


## Compare empirical and null variable importances------------------------------

# Combine empirical and null variable importance data frames, subsetting to 
# models that include taxonomic features
all.VarImp3 <- subset(all.VarImp2, Taxa == 'Including Selected Taxa')
all.VarImp3$Group <- 'Empirical'

all.VarImp_NULL3 <- subset(all.VarImp_NULL2, Taxa == 'Including Selected Taxa')
all.VarImp_NULL3$Group <- 'Null'

all.VarImp.combined <- rbind(all.VarImp3, all.VarImp_NULL3)
all.VarImp.combined$Group <- as.factor(all.VarImp.combined$Group)

all.VarImp.combined$RaceOther <- NULL


# For each model (e.g. Alpha, Beta..), for each feature, perform t-test for
# differences in importance by whether it was null (class label-shuffled) or 
# empirical version of the model. 

# Set up empty data frame to collate t-test results.
varimp.tests <- data.frame(matrix(ncol=5, nrow=195))
colnames(varimp.tests) <- c('DataSet', 'Feature', 'P.val', 'Mean.Emp', 'Mean.Null')

# Iterate through each model, and features included therein 
models <- unique(all.VarImp.combined$DataSet)
count = 0
for (i in 1:length(models)){
  varimp.data <- subset(all.VarImp.combined, DataSet == models[i])
  
  varimp.data.long <- gather(varimp.data, Feature, Feature.Importance, 
                                PET_amyloid:X.k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruminococcus.s__Ruminococcus_lactaris.,
                                factor_key=TRUE)
  varimp.data.long <- na.omit(varimp.data.long)
  
  features <- as.character(unique(varimp.data.long$Feature))
  
  for (j in 1:length(features)) {
    count = count + 1
    print(count)
    students.t <- t.test(Feature.Importance ~ Group, 
                         data = subset(varimp.data.long, Feature == features[j]))
    
    newrow <- c(DataSet = models[i], Feature = features[j], 
                P.val = students.t$p.value, Mean.Emp = students.t$estimate[[1]],
                Mean.Null = students.t$estimate[[2]])
    
    varimp.tests[count, ] <- newrow
  }
}


# Subset to comparisons in which mean Empirical is greater than mean Null
varimp.tests2 <- subset(varimp.tests, Mean.Emp > Mean.Null)

# Apply BH (FDR) correction for multiple hypothesis testing and subset to 
# significant comparisons. Add star nomenclature.
varimp.tests2$P.val.adj <- p.adjust(varimp.tests2$P.val, method='BH')
varimp.tests2 <- subset(varimp.tests2, P.val.adj < 0.05)

varimp.tests2 <- varimp.tests2 %>% 
  mutate(Adj_Stars = if_else(P.val.adj <= 0.05 & P.val.adj > 0.01, '*',
                             if_else(P.val.adj <= 0.01 & P.val.adj > 0.001, '**', 
                                     if_else(P.val.adj <= 0.001, '***', 'NS'))))



## You have reached the end of this script !! ##  
