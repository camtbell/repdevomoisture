# 1. Clear R workspace and load data
	rm(list=ls())
	pacman::p_load(interp, metaAidR, metafor, cowplot, hrbrthemes, corrplot, ggpubr, ggplot2, MCMCglmm, ape, phytools, stats, rotl, readr, MuMIn, clubSandwich, dplyr, viridis, patchwork, orchaRd, mgcv, interp, akima, fields)

	Sex_dataOR <- read_csv("3_trait data/sex_data.csv")

	Sex_dataOR <- escalc(measure="OR", ai=male.2, bi=female.2, ci =male, di=female , data=Sex_dataOR)

	Sex_dataOR$T_scaled <- (Sex_dataOR$T - 25) 
	Sex_dataOR$waterpotdiff_scaled <- (Sex_dataOR$waterpot_diff - 320)

#### 2. Create covariance matrix with RVE ####
	Vmat = impute_covariance_matrix(Sex_dataOR$vi, cluster = Sex_dataOR$paper_no, r=0.6)

#### 3. Generating the phylogenetic tree and matrix ####
	# Load the species list data 
	speciesList <- read_csv("4_species data/sex_species.csv")

	## Phylogeny
	## Access taxon relationships from Open Tree of Life; Needed for phylogenetic control
	Sex_dataOR$animal = Sex_dataOR$Genus_species
	Sex_dataOR$animal <- as.factor(Sex_dataOR$animal)

	# Match species in dataset
	tree2 = tnrs_match_names(as.character(unique(speciesList$Species)), context = "Animals")

	# Create a tree based on itt id's found on the open tree of life
	tl2 <- tol_induced_subtree(ott_ids=na.omit(tree2$ott_id)) 

	# Remove ott labels on end to make sure to matches species in dataset
	tl2$tip.label = strip_ott_ids(tl2$tip.label, remove_underscores=FALSE)

	# Getting a phylogenetic correlation matrix from the tree
	tl2_brlen <- compute.brlen(tl2, method = "Grafen", power = 0.5) 

	# Generate the phylogenetic matrix
	tl2_brlen$node.label <- NULL
	R_phylo <- vcv(tl2_brlen, corr = TRUE)

## 4. Final model based on AICc selection top model - see .R file for more details. We'll drop order through so predcitions can be marginalised over order
final.model <-rma.mv(yi = yi, V = Vmat, random = list(~1|paper_no, ~1|animal, ~1|row_count), R = list(animal=R_phylo), mod = ~ T_scaled + waterpotdiff_scaled:T_scaled + waterpotdiff_scaled, data =Sex_dataOR, method = 'REML')
final.modelb <- robust(final.model, cluster = Sex_dataOR$paper_no)
summary(final.modelb)

# Create a dataframe used for predcitions. Note that these MUST be names exactly the same as the data and be in the sam units as the data used to fit the model. So if you centred then they must be centered also. 
newdata <- as.matrix(data.frame(T_scaled = seq(min(Sex_dataOR$T_scaled), max(Sex_dataOR$T_scaled), length.out = 100), 
					  waterpotdiff_scaled = rep(seq(min(Sex_dataOR$waterpotdiff_scaled), max(Sex_dataOR$waterpotdiff_scaled), length = 100), each = 100))  %>% mutate(`T_scaled:waterpotdiff_scaled` = T_scaled*waterpotdiff_scaled))

# Then make predcitions
preds <- data.frame(predict(final.model, newmods = newdata, transf = exp, digits = 2, addx=TRUE))


## First method
fld <- with(preds, interp(x = X.T_scaled, y = X.waterpotdiff_scaled, z = preds$pred))

filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                   colorRampPalette(c("white", "blue")),
               xlab = "Temperature (centered)",
               ylab = "Moisture Difference (Centered)",
               main = "Sex Ratio",
               key.title = title(main = "Odds Ratio", cex.main = 0.9), cex.axis = 2)

# Second method
s <- interp(x = preds$X.T_scaled, y = preds$X.waterpotdiff_scaled, z = preds$pred)
image.plot(s, xlab = "Temperature (centered)", ylab = "Moisture Difference (Centered)", las = 1, col = viridis(option = "magma", 50), main = "Sex Ratio", cex.main = 2, cex.axis = 1, axis.args = list(cex.axis = 1.5), cex.lab = 1.8)
contour(s, add = TRUE, col = "white")
points(y = Sex_dataOR$T_scaled, x = Sex_dataOR$waterpotdiff_scaled, pch = 16, col = "white")


###############################dataframe code ref###########################################################
#Will be based off sex ratio; moderators being Temperature vs. kPa sig interaction.
newdata <- data.frame(temp_diff = rep(seq(15, 36, by = 1), each = 12), moist_diff = seq(100, 320, by = 20))
#Mean Centring, taken from survival_model.R
############################################################################################################

# Figure 6 - Model predictions - traits - raw
#------------------------------------------------------------------
png(height = 8.511013, width = 15.674009, res = 600, units = "in", file = "./figure/fig.6_R2.png")
	pred.dat7.6<- read.csv("./output/modelResults/pred.dat7.6")

	labcex = 0.8
	par(mfrow = c(2,3), cex.lab = 2, mar = c(6, 6, 5, 6))

	# No fluctuation in Temperature, behaviour
	pred.datBehav <- subset(pred.dat7.6, Trait_Cat == "Behaviour")
	s <- interp(x = pred.datBehav$midTemp, y = pred.datBehav$diffTemp, z = pred.datBehav$pred)
	image.plot(s, xlab = "", ylab = expression("Temperature difference"~(Delta~degree~C)), las = 1, col = viridis(option = "magma", 50), main = "Behavioural traits", cex.main = 2, cex.axis = 1.5, axis.args = list(cex.axis = 1.5))
	contour(s, add = TRUE, labcex = labcex)
	points(y = tincNoInc$T_diff[tincNoInc$Trait_Cat == "Behaviour"], x = tincNoInc$T_mid[tincNoInc$Trait_Cat == "Behaviour"], pch = 16)

	# No fluctuation in temperature, development
	pred.datDev <- subset(pred.dat7.6, Trait_Cat == "Development")
	s <- interp(x = pred.datDev$midTemp, y = pred.datDev$diffTemp, z = pred.datDev$pred)
	image.plot(s, xlab = "", ylab = "", las = 1, col = viridis(option = "magma", 50), main = "Development traits", cex.main = 2, cex.axis = 1.5, axis.args = list(cex.axis = 1.5))
	contour(s, add = TRUE, labcex = labcex)
	points(y = tincNoInc$T_diff[tincNoInc$Trait_Cat == "Development"], x = tincNoInc$T_mid[tincNoInc$Trait_Cat == "Development"], pch = 16)

	# No fluctuation in temperature, development
	pred.datMor <- subset(pred.dat7.6, Trait_Cat == "Morphology")
	s <- interp(x = pred.datMor$midTemp, y = pred.datMor$diffTemp, z = pred.datMor$pred)
	image.plot(s, xlab = "", ylab = "", las = 1, col = viridis(option = "magma", 50), main = "Morphological traits", cex.main = 2, cex.axis = 1.5, axis.args = list(cex.axis = 1.5))
	contour(s, add = TRUE, labcex = labcex)
	points(y = tincNoInc$T_diff[tincNoInc$Trait_Cat == "Morphology"], x = tincNoInc$T_mid[tincNoInc$Trait_Cat == "Morphology"], pch = 16)

	# No fluctuation in temperature, performance
	pred.datPer <- subset(pred.dat7.6, Trait_Cat == "Performance")
	s <- interp(x = pred.datPer$midTemp, y = pred.datPer$diffTemp, z = pred.datPer$pred)
	image.plot(s, xlab = expression("Mid-temperature"~(degree~C)), ylab = expression("Temperature difference"~(Delta~degree~C)), las = 1, col = viridis(option = "magma", 50), main = "Performance traits", cex.main = 2, cex.axis = 1.5, axis.args = list(cex.axis = 1.5))
	contour(s, add = TRUE, labcex = labcex)
	points(y = tincNoInc$T_diff[tincNoInc$Trait_Cat == "Performance"], x = tincNoInc$T_mid[tincNoInc$Trait_Cat == "Performance"], pch = 16)

	# No fluctuation in temperature, performance
	pred.datPhys <- subset(pred.dat7.6, Trait_Cat == "Physiology")
	s <- interp(x = pred.datPhys$midTemp, y = pred.datPhys$diffTemp, z = pred.datPhys$pred)
	image.plot(s, xlab = expression("Mid-temperature"~(degree~C)), ylab = "", las = 1, col = viridis(option = "magma", 50), main = "Physiological traits", cex.main = 2, cex.axis = 1.5, axis.args = list(cex.axis = 1.5))
	contour(s, add = TRUE, labcex = labcex)
	points(y = tincNoInc$T_diff[tincNoInc$Trait_Cat == "Physiology"], x = tincNoInc$T_mid[tincNoInc$Trait_Cat == "Physiology"], pch = 16)

	# No fluctuation in temperature, performance
	pred.datSurv <- subset(pred.dat7.6, Trait_Cat == "Survival")
	s <- interp(x = pred.datSurv$midTemp, y = pred.datSurv$diffTemp, z = pred.datSurv$pred)
	image.plot(s, xlab = expression("Mid-temperature"~(degree~C)), ylab = "", las = 1, col = viridis(option = "magma", 50), main = "Survival", cex.main = 2, cex.axis = 1.5, axis.args = list(cex.axis = 1.5))
	contour(s, add = TRUE, labcex = labcex)
	points(y = tincNoInc$T_diff[tincNoInc$Trait_Cat == "Survival"], x = tincNoInc$T_mid[tincNoInc$Trait_Cat == "Survival"], pch = 16)
dev.off()
