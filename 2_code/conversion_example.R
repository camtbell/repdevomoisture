#(For hatching success and sex ratio (sig effect from moist diff and temp vars.))
## Script to show the conversions. Let say we also put high moistute on top and low moisture on bottom
	x1 <- log(10/5)  # logRR of 0.69
	x2 <- log(5/10)  # logRR of -0.69, obviously exact opposite of x1
	x3 <- log(10/10) # logRR of 0

## Now, from the above, we know what the % difference between the numerator and denominator are, so we can calculate this
	x1_r <- (exp(x1)-1)*100 # Numerator is 100% larger than denominator. Or, high moisture increases trait by 100%
	x2_r <- (exp(x2)-1)*100 # Numerator is 50% of the denominator. Or, high moisture decreases trait by 50%, 

	
### THIS IS JLR------------------------------------------------------------------------------------------------------------------------------------
# MLMA for both non-sig and sig effects
#MLMA FOR lnOR use link below, the metafor link!
		
#MLMA FOR LENGTH (lnRR)
	mlma_length = 0.015
	(exp(mlma_length)-1)*100
	#Numerator is 1.5% larger than the denominator. Or, the high moisture treatment increase reptile length by 1.5%.
	
#MLMA FOR MASS (lnRR)
	mlma_mass = 0.031
	(exp(mlma_mass)-1)*100
	#Numerator is 3.1% larger than the denominator. Or, the high moisture treatment increases reptile mass by 3.1%
	
#MLMA FOR INCUBATION DURATION (lnRR)
	mlma_incub = 0
	(exp(mlma_incub)-1)*100
	#Numerator is 0% larger. There is no difference between the high and low moisture treatment!
	
#MLMA FOR HATCHING SUCCESS (lnOR)
	#From "survival_model.R" (predict(metamean, transf = exp, digits = 2))
	#Numerator is 0.97% larger than the denominator. High moisture increases survival by 0.97%
	
#MLMA FOR SEX RATIO (lnOR)
	#From similar formula in "sex_model.R"
	#Numerator is 0.59% larger than the denominator. High moisture increases proportion of male hatchlings by 0.59%
###----------------------------------------------------------------------------------------------------------------------------------------------------------------
		
	
### INTERPRETING META_REGRESSION COEFFICIENTS---------------------------------------------------------------------------------------------------------------
# only sig coefficients
## Slope from the model for survival (is this okay if lnOR?)
  b_temp = 0.038
  (exp(b_temp)-1)*100 
  3.8731 # When we have a moisture difference of 150kPA (isnt it 320kPa?), a 1C increase in temperature increases the trait in high moisture treatment by 3.8% relative to the low moisture treatment.

#Slope from the model for survival (moisture difference)
  b_moist = 0.001
  (exp(b_moist)-1)*100 
  0.1001 #A 1kPa increase in moisture difference between treatments increases the trait in high moisture treatment by 0.1% relative to low moisture

 #Slope from the model for sex ratio (moisture difference*temperature) 
  b_moisttemp = 0.001
  (exp(b_moisttemp)-1)*100
  0.1001 #A 1 kPa increase in moisture difference, and 1C increase in temp DECREASES the trait in high moisture by 0.1% relative to low moisture
###----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
    
 data %>% group_by(Order) %>% summarise(mean = mean(T_scaled), sd = sd(T_scaled), n = n())


# Log odds is fairly straight forward too. Check out here: https://www.metafor-project.org/doku.php/tips:assembling_data_or

