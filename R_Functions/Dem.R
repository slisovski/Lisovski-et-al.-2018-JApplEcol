library(doMC)
	registerDoMC(4)

mu_j <- function(my_a = 0.315/360, start_pop = 600, hatchlings = 4, breed_mean = 150, breed_sd = 10, trans_day = 93) {
	
	mod <- function(mu, b) {
	
		t <- 1:(365*3)
		xstart <- c(A = start_pop, J = 0)
		
		
		f_B  <- approxfun(
	                     x = t,
	                     y = rep(dnorm(1:365, b, breed_sd), max(ceiling(t/365)))[t],
	                     rule = 2
						 )
		
		
		event_times <- t[which(rep(1:365, max(ceiling(t/365)))[1:length(t)] == trans_day)]
		eventfun <- function(t, y, parms) {
	  			with(as.list(y), { 
	   					A <- A+J
	    				J <- 0
	    				c(A, J)
	  					})
			    	}
	
		unlist(foreach(i = mu) %dopar% {
			parms <-  c(mu_a = my_a, mu_j = i, pop = start_pop, hatch = hatchlings) 
			out1 <-  ode(y = xstart, t, func = "popM", parms = parms, 
				   		 dllname = "dem", initforc = "initforc",
				   		 forcings = list(cbind(t, f_B(t))),
				   		 initfunc = "initmod", nout = 0,
				   		 events = list(func = "event", times = event_times),
				   		 method = rkMethod("rk34f"))	
			diff(apply(out1[,2:3], 1, sum)[event_times+3])[2]
		})
	}	
	
uniroot.all(f = mod, interval = c(0.315/360, 1/20), b = breed_mean)
}
