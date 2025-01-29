# Trait Analysis for Moisture Effects: Length

# Clear R workspace
rm(list=ls())

# Loading R Packages
pacman::p_load(metaAidR, metafor, cowplot, hrbrthemes, corrplot, ggpubr, ggplot2, MCMCglmm, ape, phytools, stats, rotl, readr, MuMIn, clubSandwich, dplyr, viridis, patchwork, orchaRd)

# Load Length Dataset
L_dataRR <- read.csv("3_trait data/length_data.csv")

# Check that continuous variables are not categorical
range(L_dataRR$waterpot_diff)
hist(L_dataRR$waterpot_diff)

range(L_dataRR$waterpot_mid)
hist(L_dataRR$waterpot_mid)

range(L_dataRR$T)
hist(L_dataRR$T)

range(L_dataRR$egg_mass_g)
hist(L_dataRR$egg_mass_g)

# Make sure that appropriate variables are coded as factors
L_dataRR$order <- as.factor(L_dataRR$order)

# Check data types in dataset
str(L_dataRR)


#### 1. Calculate effect size ####
L_dataRR<-escalc(measure = "ROM", n1i = N.2, n2i =N, m1i =mean.2, m2i =mean, sd1i =sd.2, sd2i =sd, data = L_dataRR)

# Checks variable types (integer, character, numeric)
str(L_dataRR)

# Sample sizes for table 1
length(L_dataRR$yi)
length(unique(L_dataRR$paper_no))
length(unique(L_dataRR$species))


#### 2. Create covariance matrix with RVE ####
Vmat = impute_covariance_matrix(L_dataRR$vi, cluster = L_dataRR$paper_no, r= 0.6)


#### 3. Heterogeneity: calculating and checking the I2 ####
# Fit the meta-analysis that accounts for sampling covariance
model_test <- rma.mv(yi ~ 1, V=Vmat , random= ~ 1 | row_count, data=L_dataRR)
summary(model_test)
sigma2 <- 0.0029 

# Calculate I2
i2 <- orchaRd::i2_ml(model_test, boot = 1000) 
i2 
# I2_Total; I2 = 92.144, 95% CI = 89.875, 93.666


#### 4. Generating the phylogenetic tree and matrix ####
# Load the species list data 
speciesList <- read_csv("4_species data/length_species.csv")

## Phylogeny
## Access taxon relationships from Open Tree of Life; Needed for phylogenetic control
L_dataRR$animal = L_dataRR$genus_species
L_dataRR$animal<-as.character(L_dataRR$animal)

# Match species in dataset
tree2 = tnrs_match_names(as.character(unique(speciesList$Species)), context_name = "Animals")

# Create a tree based on itt id's found on the open tree of life
tl2 <- tol_induced_subtree(ott_ids=na.omit(tree2$ott_id)) 

# Remove ott labels on end to make sure to matches species in dataset
tl2$tip.label = strip_ott_ids(tl2$tip.label, remove_underscores=FALSE)

# Getting a phylogenetic correlaiton matrix from the tree
tl2_brlen <- compute.brlen(tl2, method = "Grafen", power = 0.5) 

# Check trees are ultrametric
is.ultrametric(tl2_brlen)

# Plot the tree
plot(tl2, cex = 0.8, label.offset = .1, no.margin = TRUE)

# Generate the phylogenetic matrix
tl2_brlen$node.label <- NULL
R_phylo <- vcv(tl2_brlen, corr = TRUE)
print(R_phylo)

# Checking the spelling etc matches due to error generating the phylogenetic matrix
unique(L_dataRR$genus_species)
rownames(R_phylo) %in% unique(L_dataRR$genus_species)
rownames(R_phylo)[!rownames(R_phylo) %in% unique(L_dataRR$genus_species)]


#### 6. Meta-regressions to determine effect of moderators ####
# Get MuMIn to work with metafor
eval(metafor:::.MuMIn)

# z-transforming (mean centering) the continuous variables so they are on the same scale and able to be compared
L_dataRR$egg_z <-scale(L_dataRR$egg_mass_g) #z-scaled
L_dataRR$T_scaled <- (L_dataRR$T - 25) # what this will do is make the value of 0 mean that the intercept is 25 C
L_dataRR$waterpotdiff_scaled <- (L_dataRR$waterpot_diff - 320) # what this will do is make the moisture difference of 0 mean a value of 320 KPA 

# Full model including all random effects and moderators
mod <- rma.mv(yi = yi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = ~ order + egg_z + T_scaled + waterpotdiff_scaled:T_scaled + waterpotdiff_scaled, data = L_dataRR, method = 'ML')
summary(mod)
robustmod <- robust(mod, cluster = L_dataRR$paper_no)
summary(robustmod)

# AICc model selection
AICcAlt <- dredge(robustmod, trace=2)
subset(AICcAlt, delta <= 2, recalc.weights=FALSE)
AICcAlt
-417.1-(-416.6)

## Final model based on AICc selection top model (null model has at least temperature and water potential diff included)
final.model <- rma.mv(yi = yi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = ~ T_scaled + waterpotdiff_scaled, data = L_dataRR, method = 'REML')
final.modelb <- robust(final.model, cluster = L_dataRR$paper_no)
summary(final.modelb)

## Calculating the meta-analytic mean (overall effect size) 
meta_all = rma.mv(yi = yi, V = Vmat, random = list(~1|animal, ~1|paper_no, ~1|row_count), R = list(animal=R_phylo), data=L_dataRR, method = 'REML')
meta_allb <- robust(meta_all, cluster = L_dataRR$paper_no)
summary(meta_allb)

## MEAN MLMA 
M <- coef(meta_all)
M

## Calculate the prediction intervals for the MLMA
predictionsMLMR <- predict(meta_all)
print(predictionsMLMR)


#### 7. Publication bias ####

# Calculating weighted average (or precision)
L_dataRR$wi<-  1/(L_dataRR$vi) 

# z-transform pub_year
L_dataRR$pub_year <-scale(L_dataRR$pub_year) 

# Model for publication bias includes all moderators and random effects, as well as year to account for time lag and precision (wi)
pub.model_z <-rma.mv(yi ~ waterpotdiff_scaled + T_scaled + waterpotdiff_scaled:T_scaled + pub_year + vi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = , data = L_dataRR, method = 'ML')
summary(pub.model_z)
