#(For hatching success and sex ratio (sig effect from moist diff and temp vars.))
## Script to show the conversions. Let say we also put high moistute on top and low moisture on bottom
	x1 <- log(10/5)  # logRR of 0.69
	x2 <- log(5/10)  # logRR of -0.69, obviously exact opposite of x1
	x3 <- log(10/10) # logRR of 0

## Now, from the above, we know what the % difference between the numerator and denominator are, so we can calculate this
	x1_r <- (exp(x1)-1)*100 # Numerator is 100% larger than denominator. Or, high moisture increases trait by 100%
	x2_r <- (exp(x2)-1)*100 # Numerator is 50% of the denominator. Or, high moisture decreases trait by 50%, 

	
### THIS IS JLR
# MLMA for both non-sig and sig effects
#MLMA FOR lnOR use link below
		
#MLMA FOR LENGTH (lnRR)
	mlma_length = 0.015
	(exp(mlma_length)-1)*100
	#Numerator is 1.5% larger than the denominator. Or, the high moisture treatment increase reptile length by 1.5%.
	
### INTERPRETING META_REGRESSION COEFFICIENTS
# only sig coefficients
## Slope from the model for survival (is this okay if lnOR?)
  b_temp = 0.038
 (exp(b_temp)-1)*100 # When we have a moisture difference of 150kPA, a 1C increase in temperature increases the trait in high moisture treatment by 3.8% relative to the low moisture treatment.

  
  
  
 # A 1C increase in temperature increases the trait in the 1500kpA treatment by 3% relative to the 1(500-150)kPA moisture treatment.


 # 
 data %>% group_by(Order) %>% summarise(mean = mean(T_scaled), sd = sd(T_scaled), n = n())


# Log odds is fairly straight forward too. Check out here: https://www.metafor-project.org/doku.php/tips:assembling_data_or

