library(ggplot)
library(gmcontour)
library(mcgv)

# Create a dataframe used for predcitions. Note that these MUST be names exactly the same as the data and be in the sam units as the data used to fit the model. So if you centred then they must be centered also. 
newdata <- data.frame(temp_diff = rep(seq(15, 36, by = 1), each = 12), moist_diff = seq(100, 320, by = 20))

# Then make predcitions
preds <- predict(metamean, newdata = newdata, transf = exp, digits = 2)

# Then you use this griant dataframe with x (temperature), y (kPa diff) and z (maybe effect size?) values to make a contour plot

###############################Cam's Attempt.###############################################################
#Will be based off sex ratio; moderators being Temperature vs. kPa sig interaction.
#Therefore:
Sex_dataOR <- read_csv("3_trait data/sex_data.csv")

newdata <- data.frame(temp_diff = rep(seq(15, 36, by = 1), each = 12), moist_diff = seq(100, 320, by = 20))
#Mean Centring, taken from survival_model.R
Sur_dataOR$T_scaled <- (Sur_dataOR$T - 25) 
Sur_dataOR$waterpotdiff_scaled <- (Sur_dataOR$waterpot_diff - 320)
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
