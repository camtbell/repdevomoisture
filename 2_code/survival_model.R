# Trait Analysis for Moisture Effects: Hatching Success

# Clear R workspace
rm(list=ls())

# Loading R Packages
pacman::p_load(metaAidR, metafor, cowplot, hrbrthemes, corrplot, ggpubr, ggplot2, MCMCglmm, ape, phytools, stats, rotl, readr, MuMIn, clubSandwich, dplyr, viridis, patchwork, orchaRd)

# Load Length Dataset
Sur_dataOR <- read_csv("3_trait data/survival_data.csv")

# Check that continuous variables are not categorical
range(Sur_dataOR$waterpot_diff)
hist(Sur_dataOR$waterpot_diff)

range(Sur_dataOR$waterpot_mid)
hist(Sur_dataOR$waterpot_mid)

range(Sur_dataOR$T)
hist(Sur_dataOR$T)

range(Sur_dataOR$egg_mass_g)
hist(Sur_dataOR$egg_mass_g)

# Make sure that appropriate variables are coded as factors
Sur_dataOR$order <- as.factor(Sur_dataOR$order)

# Check data types in dataset
str(Sur_dataOR)


#### 1. Calculate effect size ####
Sur_dataOR <- escalc(measure="OR", ai=alive.2, bi=dead.2, ci =alive.1, di =dead.1 , data=Sur_dataOR)

# Sample sizes for table 1
length(Sur_dataOR$yi)
length(unique(Sur_dataOR$paper_no))
length(unique(Sur_dataOR$species))


#### 2. Create covariance matrix with RVE ####
Vmat = impute_covariance_matrix(Sur_dataOR$vi, cluster = Sur_dataOR$paper_no, r= 0.6)


#### 3. Heterogeneity: calculating and checking the I2 ####
# Fit the meta-analysis that accounts for sampling covariance
model_test <- rma.mv(yi ~ 1, V=Vmat , random= ~ 1 | row_count, data=Sur_dataOR)
summary(model_test)
sigma2 <- 0.7572

# Calculate I2
i2 <- orchaRd::i2_ml(model_test, boot = 1000, data = Sur_dataOR) 
i2
# I2_Total; I2 = 50.463, 95% CI = 37.079, 64.572


#### 4. Generating the phylogenetic tree and matrix ####
# Load the species list data 
speciesList <- read_csv("4_species data/survival_species.csv")

## Phylogeny
## Access taxon relationships from Open Tree of Life; Needed for phylogenetic control
Sur_dataOR$animal = Sur_dataOR$Genus_species
Sur_dataOR$animal<-as.factor(Sur_dataOR$animal)

# Match species in dataset
tree2 = tnrs_match_names(as.character(unique(speciesList$Species)), context = "Animals")

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
unique(Sur_dataOR$Genus_species)
rownames(R_phylo) %in% unique(Sur_dataOR$Genus_species)
rownames(R_phylo)[!rownames(R_phylo) %in% unique(Sur_dataOR$Genus_species)]


#### 6. Meta-regressions to determine effect of moderators ####
# Get MuMIn to work with metafor
eval(metafor:::.MuMIn)

# z-transforming egg mass and mean centering temperature and water potential difference
Sur_dataOR$egg_z <-scale(Sur_dataOR$egg_mass_g) #z-scaled
Sur_dataOR$T_scaled <- (Sur_dataOR$T - 25) # what this will do is make the value of 0 mean that the intercept is 25 C
Sur_dataOR$waterpotdiff_scaled <- (Sur_dataOR$waterpot_diff - 320) # what this will do is make the moisture difference of 0 mean a value of 320 KPA

# Full model including all random effects and moderators 
full.model_z <- rma.mv(yi = yi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = ~ order + T_scaled + waterpotdiff_scaled:T_scaled + waterpotdiff_scaled + egg_z, data =Sur_dataOR, method = 'ML')
summary(full.model_z)
robustmod <- robust(full.model_z, cluster = Sur_dataOR$paper_no)
summary(robustmod)

## AICc model selection
AICcAlt <- dredge(robustmod, trace=2)
subset(AICcAlt, delta <= 2, recalc.weights=FALSE)
AICcAlt
643.2-642.4

## Final model based on AICc selection top model (null = models with temp and water diff in them)
final.model <-rma.mv(yi = yi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = ~ waterpotdiff_scaled + T_scaled, data =Sur_dataOR, method = 'REML')
final.modelb <- robust(final.model, cluster = Sur_dataOR$paper_no)
summary(final.modelb)

## Calculating the meta-analytic mean (overall effect size)
metamean <- rma.mv(yi ~ 1,V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), data = Sur_dataOR, method = 'REML')
metameanb <- robust(metamean, cluster = Sur_dataOR$paper_no)
summary(metameanb)

#Conversion Trial Block
predict(metamean, transf = exp, digits = 2)

## MEAN MLMA 
M <- coef(metamean)
M

## Calculate the prediction intervals for the MLMA
predictionsMLMR <- predict(metamean)
print(predictionsMLMR)


#### 7. Publication bias ####

# Calculating weighted average (or precision)
Sur_dataOR$wi<-  1/(Sur_dataOR$vi) 

#z-transforming pub_year
Sur_dataOR$pub_year <-scale(Sur_dataOR$pub_year) 

# Model for publication bias includes all moderators and random effects, as well as year to account for time lag and precision (wi)
pub.model_z <-rma.mv(yi ~ waterpotdiff_scaled + T_scaled + waterpotdiff_scaled:T_scaled + pub_year + vi + wi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = , data = Sur_dataOR, method = 'ML')
summary(pub.model_z)
