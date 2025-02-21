# Trait Analysis for Moisture Effects: Incubation Duration

# Clear R workspace
rm(list=ls())

# Loading R Packages
pacman::p_load(metaAidR, metafor, cowplot, hrbrthemes, corrplot, ggpubr, ggplot2, MCMCglmm, ape, phytools, stats, rotl, readr, MuMIn, clubSandwich, dplyr, viridis, patchwork, orchaRd)

# Load incubation duration Dataset
ID_dataRR <- read.csv("3_trait data/incubation_data.csv")

# Check that continuous variables are not categorical
range(ID_dataRR$waterpot_diff)
hist(ID_dataRR$waterpot_diff)

range(ID_dataRR$waterpot_mid)
hist(ID_dataRR$waterpot_mid)

range(ID_dataRR$T)
hist(ID_dataRR$T)

range(ID_dataRR$egg_mass_g)
hist(ID_dataRR$egg_mass_g)

# Make sure that appropriate variables are coded as factors
ID_dataRR$order <- as.factor(ID_dataRR$order)

# Check data types in dataset
str(ID_dataRR)


#### 1. Calculate effect size ####
ID_dataRR <- escalc(measure = "ROM", n1i = N.2, n2i =N, m1i =mean.2, m2i =mean, sd1i =sd.2, sd2i =sd, data = ID_dataRR)

# Checks variable types (integer, character, numeric)
str(ID_dataRR)

# Sample sizes for table 1
length(ID_dataRR$yi)
length(unique(ID_dataRR$paper_no))
length(unique(ID_dataRR$species))


#### 2. Create covariance matrix with RVE ####
Vmat = impute_covariance_matrix(ID_dataRR$vi, cluster = ID_dataRR$paper_no, r= 0.6)


#### 3. Heterogeneity: calculating and checking the I2 ####
# Fit the meta-analysis that accounts for sampling covariance
model_test <- rma.mv(yi ~ 1, V=Vmat , random= ~ 1 | row_count, data=ID_dataRR)
summary(model_test)
sigma2 <- 0.0010 

# Calculate I2
i2 <- orchaRd::i2_ml(model_test, boot = 1000) 
i2
# I2_Total; I2 = 92.52, 95% CI = 89.827, 94.512


#### 4. Generating the phylogenetic tree and matrix ####
# Load the species list data 
speciesList <- read_csv("4_species data/incubation_species.csv")

## Phylogeny
## Access taxon relationships from Open Tree of Life; Needed for phylogenetic control
ID_dataRR$animal = ID_dataRR$genus_species
ID_dataRR$animal <- as.character(ID_dataRR$animal)

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
unique(ID_dataRR$genus_species)
rownames(R_phylo) %in% unique(ID_dataRR$genus_species)
rownames(R_phylo)[!rownames(R_phylo) %in% unique(ID_dataRR$genus_species)]


#### 6. Meta-regressions to determine effect of moderators ####
# Get MuMIn to work with metafor
eval(metafor:::.MuMIn)

# z-transforming egg mass and mean centering temperature and water potential difference
ID_dataRR$egg_z <- scale(ID_dataRR$egg_mass_g) #z-scaled
ID_dataRR$T_scaled <- (ID_dataRR$T - 25) # what this will do is make the value of 0 mean that the intercept is 25 C
ID_dataRR$waterpotdiff_scaled <- (ID_dataRR$waterpot_diff - 320) # what this will do is make the moisture difference of 0 mean a value of 320 KPA

# Full model including all random effects and moderators 
mod <- rma.mv(yi = yi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = ~ order + egg_z + waterpotdiff_scaled + T_scaled + waterpotdiff_scaled*T_scaled, data = ID_dataRR, method = 'ML')
robustmod <- robust(mod, cluster = ID_dataRR$paper_no)
summary(robustmod)

## AICc model selection
AICcAlt <- dredge(robustmod, trace=2)
subset(AICcAlt, delta <= 2, recalc.weights=FALSE)
466.4-465.3
AICcAlt

## Final model based on AICc selection top model (null model has at least temperature and water potential diff included)
final.model <- rma.mv(yi = yi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = ~ waterpotdiff_scaled + T_scaled + waterpotdiff_scaled:T_scaled, data = ID_dataRR, method = 'REML')
final.modelb <- robust(final.model, cluster = ID_dataRR$paper_no)
summary(final.modelb)

## Calculating the meta-analytic mean (overall effect size) 
meta_all = rma.mv(yi = yi, V = Vmat, random = list(~1|animal, ~1|paper_no, ~1|row_count), R = list(animal=R_phylo), data=ID_dataRR, method = 'REML')
meta_allb <- robust(meta_all, cluster = ID_dataRR$paper_no)
summary(meta_allb)

## MEAN MLMA 
M <- coef(meta_allb)
M

## Calculate the prediction intervals for the MLMA
predictionsMLMR <- predict(meta_allb)
print(predictionsMLMR)


#### 7. Publication bias ####

# Calculating weighted average (or precision)
ID_dataRR$wi <-  1/(ID_dataRR$vi) 

# z-transform publication year
ID_dataRR$pub_year <- scale(ID_dataRR$pub_year) 

# Model for publication bias includes all moderators and random effects, as well as year to account for time lag and precision (wi)
pub.model_z <-rma.mv(yi ~ waterpotdiff_scaled + T_scaled + waterpotdiff_scaled:T_scaled + pub_year + vi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = , data = ID_dataRR, method = 'ML')
summary(pub.model_z)
