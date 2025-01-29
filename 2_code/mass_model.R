# Trait Analysis for Moisture Effects: Mass

# Clear R workspace
rm(list=ls())

# Loading R Packages
pacman::p_load(metaAidR, metafor, cowplot, hrbrthemes, corrplot, ggpubr, ggplot2, MCMCglmm, ape, phytools, stats, rotl, readr, MuMIn, clubSandwich, dplyr, viridis, patchwork, orchaRd)

# Load Length Dataset
M_dataRR <- read.csv("3_trait data/mass_data.csv")

# Check that continuous variables are not categorical
range(M_dataRR$waterpot_diff)
hist(M_dataRR$waterpot_diff)
mean(M_dataRR$waterpot_diff)

range(M_dataRR$waterpot_mid)
hist(M_dataRR$waterpot_mid)

range(M_dataRR$T)
hist(M_dataRR$T)
mean(M_dataRR$T)

range(M_dataRR$egg_mass_g)
hist(M_dataRR$egg_mass_g)

# Make sure that appropriate variables are coded as factors
M_dataRR$order <- as.factor(M_dataRR$order)

# Check data types in dataset
str(M_dataRR)


#### 1. Calculate effect size ####
M_dataRR <- escalc(measure = "ROM", n1i = N.2, n2i =N, m1i =mean.2, m2i =mean, sd1i =sd.2, sd2i =sd, data = M_dataRR)

# Checks variable types (integer, character, numeric)
str(M_dataRR)

# Sample sizes for table 1
length(M_dataRR$yi)
length(unique(M_dataRR$paper_no))
length(unique(M_dataRR$species))


#### 2. Create covariance matrix with RVE ####
Vmat = impute_covariance_matrix(M_dataRR$vi, cluster = M_dataRR$paper_no, r= 0.6)

#### 3. Heterogeneity: calculating and checking the I2 ####
# Fit the meta-analysis that accounts for sampling covariance
model_test <- rma.mv(yi ~ 1, V=Vmat , random= ~ 1 | row_count, data=M_dataRR)
summary(model_test)
sigma2 <- 0.0128 

# Calculate I2. You want I2 total. You can also bootstrap.
i2 <- orchaRd::i2_ml(model_test, boot = 1000) 
i2
# I2_Total; I2 = 93.888, 95% CI = 92.264, 95.331


#### 4. Generating the phylogenetic tree and matrix ####
# Load the species list data 
speciesList <- read_csv("4_species data/mass_species.csv")

## Phylogeny
## Access taxon relationships from Open Tree of Life; Needed for phylogenetic control
M_dataRR$animal = M_dataRR$genus_species
M_dataRR$animal<-as.character(M_dataRR$animal)

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
unique(M_dataRR$genus_species)
rownames(R_phylo) %in% unique(M_dataRR$genus_species)
rownames(R_phylo)[!rownames(R_phylo) %in% unique(M_dataRR$genus_species)]


#### 6. Meta-regressions to determine effect of moderators ####
# Get MuMIn to work with metafor
eval(metafor:::.MuMIn)

# z-transforming (mean centering) the continuous variables so they are on the same scale and able to be compared
M_dataRR$egg_z <- scale(M_dataRR$egg_mass_g) #z-scaled
M_dataRR$T_scaled <- (M_dataRR$T - 25) # what this will do is make the value of 0 mean that the intercept is 25 C
M_dataRR$waterpotdiff_scaled <- (M_dataRR$waterpot_diff - 320) # what this will do is make the moisture difference of 0 mean a value of 320 KPA

# Full model including all random effects and moderators 
mod <- rma.mv(yi = yi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = ~ order + egg_z + waterpotdiff_scaled + T_scaled + waterpotdiff_scaled:T_scaled, data = M_dataRR, method = 'ML')
robustmod <- robust(mod, cluster = M_dataRR$paper_no)
summary(robustmod)

## AICc model selection
AICcAlt <- dredge(robustmod, trace=2)
subset(AICcAlt, delta <= 2, recalc.weights=FALSE)
AICcAlt
(-283.0)-(-281.8)

## Final model based on AICc selection top model (null = models with temp and water diff in them)
final.model <- rma.mv(yi = yi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = ~ waterpotdiff_scaled + T_scaled, data = M_dataRR, method = 'REML')
final.modelb <- robust(final.model, cluster = M_dataRR$paper_no)
summary(final.modelb)

## Calculating the meta-analytic mean (overall effect size) 
meta_all = rma.mv(yi = yi, V = Vmat, random = list(~1|animal, ~1|paper_no, ~1|row_count), R = list(animal=R_phylo), data=M_dataRR, method = 'REML')
meta_allb <- robust(meta_all, cluster = M_dataRR$paper_no)
summary(meta_allb)

## MEAN MLMA 
M <- coef(meta_all)
M

## Calculate the prediction intervals for the MLMA
predictionsMLMR <- predict(meta_all)
print(predictionsMLMR)


#### 7. Publication bias ####

# Calculating weighted average (or precision)
M_dataRR$wi <-  1/(M_dataRR$vi) 

# z-transform publication year
M_dataRR$pub_year_z <- scale(M_dataRR$pub_year) 

# Model for publication bias includes all moderators and random effects, as well as year to account for time lag and precision (wi)
pub.model_z <-rma.mv(yi ~ waterpotdiff_scaled + T_scaled + waterpotdiff_scaled:T_scaled + pub_year + vi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = , data = M_dataRR, method = 'ML')
summary(pub.model_z)
