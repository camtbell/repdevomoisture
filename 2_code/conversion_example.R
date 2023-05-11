
## Script to show the conversions. Let say we also put high moistute on top and low moisture on bottom
	x1 <- log(10/5)  # logRR of 0.69
	x2 <- log(5/10)  # logRR of -0.69, obviously exact opposite of x1
	x3 <- log(10/10) # logRR of 0

## Now, from the above, we know what the % difference between the numerator and denominator are, so we can calculate this
	x1_r <- (exp(x1)-1)*100 # Numerator is 100% larger than demoninator. Or, high moisture increases trait by 100%
	x2_r <- (exp(x2)-1)*100 # Numerator is 50% of the demoninator. Or, high moisture decreases trait by 50%, 


# Log odds is fairly straight forward too. Check out here: https://www.metafor-project.org/doku.php/tips:assembling_data_or

