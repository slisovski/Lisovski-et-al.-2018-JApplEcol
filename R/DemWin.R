########################
### Parallel ###########
########################


parfun <- function(x, xstart,
                      t,
                      mu_a, 
                      start_pop, 
                      hatchlings, 
                      forc, 
                      event_times) {
  
  library(deSolve)
  dyn.load("CCode/dem.dll")
  
  
  parms <-  c(mu_a = mu_a, mu_j = x, pop = start_pop, hatch = hatchlings) 
  out1 <-  ode(y = xstart, t, func = "popM", parms = parms,
                 dllname = "dem", initforc = "initforc",
                 forcings = forc,
                 initfunc = "initmod", nout = 0,
                 events = list(func = "event", times = event_times),
                 method = rkMethod("rk34f"))
  diff(apply(out1[,2:3], 1, sum)[event_times+3])[2]
  
}
  


mod <- function(x, 
                mu_a,
                start_pop,
                hatchlings,
                breed_mean,
                breed_sd,
                trans_day) {
  

  t <- 1:(365*3)
  xstart <- c(A = start_pop, J = 0)
  
  
  f_B  <- approxfun(
    x = t,
    y = rep(dnorm(1:365, breed_mean, breed_sd), max(ceiling(t/365)))[t],
    rule = 2
  )
  
  event_times <- t[which(rep(1:365, max(ceiling(t/365)))[1:length(t)] == trans_day)]
  
  parSapply(cl, X = x, FUN = parfun, xstart =  xstart, t = t, mu_a = mu_a,  start_pop = start_pop, hatchlings = hatchlings,  
                                     forc = list(cbind(t, f_B(t))), 
                                     event_times = event_times)
  
}



mu_j <- function(my_a = 0.315/360, start_pop = 600, hatchlings = 4, breed_mean = 150, breed_sd = 10, trans_day = 93) {

uniroot.all(f = mod, interval = c(0.315/360, 1/20), mu_a = my_a,
                                                    start_pop = start_pop,
                                                    hatchlings = hatchlings,
                                                    breed_mean = breed_mean,
                                                    breed_sd = breed_sd,
                                                    trans_day = trans_day)
}