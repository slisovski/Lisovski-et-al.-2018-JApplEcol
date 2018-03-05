sirs <- function(p, mod, aggr = TRUE) {
	
	## Time ---
	x <- 1:(365*10)
  	years <- ceiling(max(x)/365)
  	yday <- rep(1:365, years)
  	t <- seq(1, length(yday))
  	
    ## Fixed parameters
	
		## Epidemiology
		h = 5e-3
			
		## Demographie
		start_pop = 600
		departure_day = 60
		mort = 0.315/365

	beta = as.numeric(p[1])
	beta_r = as.numeric(p[2])
	beta_m = as.numeric(p[3])
	gamma = as.numeric(p[4])
	sigma = as.numeric(p[5])
	
	hatchlings = ifelse(unlist(strsplit(mod, split = ""))[2]=="1", 4, 0.63)
	breed_mean = as.numeric(p[6])
	breed_sd = as.numeric(p[7])
	mort_j <- ifelse(unlist(strsplit(mod, split = ""))[2]=="1", mu_j(breed_mean = as.numeric(breed_mean), breed_sd = as.numeric(breed_sd)), NA)
	
	Pr_migrants = as.numeric(p[8])
	arrival_mean = as.numeric(p[9])
	arrival_sd = as.numeric(p[10])
	turn_mean = as.numeric(p[11])
	turn_amp = as.numeric(p[12])
	turn_slope = as.numeric(p[13])
	turn_kurt = as.numeric(p[14])


	# Forcing functions ---
	f_Birth <- approxfun(
                    x = t,
                    y = rep(dnorm(1:365, breed_mean, breed_sd), max(ceiling(t/365)))[t],
                    rule = 2
                    )
	f_Migr  <-  approxfun(
                    x = t,
                    y = rep(dnorm(1:365, arrival_mean, arrival_sd), max(ceiling(t/365)))[t],
                    rule = 2
                    )
    
    if(unlist(strsplit(mod, split = ""))[6]=="1") {

    	 	turn_curve0 <- c(turn_amp*(exp(-((turn_mean - (1:turn_mean))/turn_slope)^turn_kurt)))
  			turn_curve1 <- c(turn_curve0, rev(turn_curve0)[2:length(turn_curve0)])[1:365]

  	f_Turn  <- approxfun(
  					 x = t,
                     y = rep(turn_curve1, max(ceiling(t/365)))[t], rule = 2
                     )
    }
    
    
    event_times <- t[which(yday==departure_day)]
    
    # Model ini ----                
	  if(unlist(strsplit(mod, split = ""))[2]=="1") {
	  xstart <- c(Sa = start_pop, Ia = 0, Ra = 0, Sj = 0, Ij = 0, Rj = 0, Sm = 0, Im = 0, Rm = 0)	
	  } else {
	  xstart <- c(Sa = start_pop, Ia = 0, Ra = 0, Sm = 0, Im = 0, Rm = 0)		
	  }
	
	  parms0 <- c(start_pop = start_pop, hatchlings = hatchlings, mort = mort, mort_j = mort_j,
		  		Pr_migrants = Pr_migrants, beta = beta, beta_r = beta_r, beta_m = beta_m, 
			  	gamma = gamma, sigma = sigma, h = h)
	  parms  <- parms0[!is.na(parms0)]
	
	
	if(unlist(strsplit(mod, split = ""))[6]=="0") {
		out <- ode(y = xstart, t, func = "modSIRS", parms = parms, 
			   dllname = mod, initforc = "initforc",
			   forcings = list(cbind(t, f_Birth(t)), cbind(t, f_Migr(t))),
			   initfunc = "initmod", nout = 0,
			   events = list(func = "event", times = event_times),
			   method = rkMethod("rk34f"))
		} else {
		out <- ode(y = xstart, t, func = "modSIRS", parms = parms, 
			   dllname = mod, initforc = "initforc",
			   forcings = list(cbind(t, f_Birth(t)), cbind(t, f_Migr(t)), cbind(t, f_Turn(t))),
			   initfunc = "initmod", nout = 0,
			   events = list(func = "event", times = event_times),
			   method = rkMethod("rk34f"))
		}
	
	out <- tail(out, 365)
	out[,1] <- 1:365
	
out
}