################################################################
#  Project ID:  Mallard Simulation Study                	   #
#  Date: 08/04/2015                                            #
#  Author: Simeon Lisovski                                     #
################################################################

library(deSolve)
library(FME)
library(msir)

source("R_Functions/Dem.R")
source("SIRS.R")

dir.create("/Data/SimResults")


  binom.low <- function(x) ifelse(!is.na(x[2]), binom.test(x[1], x[2])$conf.int[1], NA)
  binom.upp <- function(x) ifelse(!is.na(x[2]), binom.test(x[1], x[2])$conf.int[2], NA)

Data0  <- read.csv("Data/Mallard_NL.csv")
Data   <- data.frame(Data0[,1], Data0[c(3:12, 1:2), 2:ncol(Data0)])
  n    <- apply(Data[,6:9], 1, sum, na.rm = T)
  pos  <- apply(Data[,2:5], 1, sum, na.rm = T)
  prev <- pos/n

matplot(Data[,1], cbind(prev, apply(cbind(pos,n), 1, binom.low), apply(cbind(pos,n), 1, binom.upp)), 
			type="l", lwd = c(2,1,1), lty = c(1,2,2), col = c("black", "grey40", "grey40"),
			xlab = "", ylab = "prevalence", xaxt = "n")
axis(1, at = 1:12, labels = format(seq(as.POSIXct("2012-03-01"), as.POSIXct("2013-02-01"), by = "month"), "%b"))
			
	
  n_r    <- apply(Data[,c("N_juv", "N_ad")], 1, sum, na.rm = T)
  pos_r  <- apply(Data[,c("Pos_juv", "Pos_ad")], 1, sum, na.rm = T)

  n_m    <- Data[,"N_mig"]
  pos_m  <- Data[,"Pos_mig"]
  
  n_u    <- Data[,"N_u"]
  pos_u  <- Data[,"Pos_u"]	

Mod <- read.csv("Data/Parameters.csv")
out <- data.frame(Mod = rep(NA, nrow(Mod)), Nr.success = rep(NA, nrow(Mod)), LogLik = rep(NA, nrow(Mod)))


#### Load .c functions
dyn.load("C_Functions/dem.so")
#######################


#_________________________________________________________________________
#------------------------------------------------------------------------#
#--      Simulations                                           		   --# 
#------------------------------------------------------------------------#
for(i in 1:nrow(Mod)) {
	system(paste("R CMD SHLIB C_Functions/", paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = ""), ".c", sep = ""))
	}

for(i in 1:nrow(Mod)) {
	out[i,1] <- paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = "")
}


for(i in 1:nrow(Mod)) {
	
	if(is.na(out[i,3])) {
		
	dyn.load(paste("C_Functions/", paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = ""), ".so", sep = ""))
	
	parId <- as.numeric(unlist(strsplit(as.character(Mod$parId[i]), ".", fixed = T)))
	
			 		 
	mod_cost_nll <- function(lpars) {
		pars <- rep(NA, 14)
		pars[parId]  <- lpars
		out <- sirs(p = pars, mod = paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = ""))

			mday  <- which(as.numeric(format(as.POSIXlt(as.POSIXct("2012-01-01 00:00:00")+((out[,1]-1)*24*60*60)), "%d"))==15)
		
		if(ncol(out)==10) {
			Prev_a <- apply(out[mday,c(3,6,9)], 1, sum)/apply(out[mday, 2:10],1,sum)
			Prev_r <- apply(out[mday,c(3,6)], 1, sum)/apply(out[mday, 2:7],1,sum)
			Prev_m <- ifelse(apply(out[mday,8:10],1,sum)<1, 5e-3, out[mday,9])/ifelse(apply(out[mday,8:10],1,sum)<1, 1, apply(out[mday,8:10],1,sum))
		} else {
			Prev_a <- apply(out[mday,c(3,6)], 1, sum)/apply(out[mday, 2:7],1,sum)
			Prev_r <- ifelse(apply(out[mday,2:4],1,sum)<1, 5e-3, out[mday,3])/ifelse(apply(out[mday,2:4],1,sum)<1, 1, apply(out[mday,2:4],1,sum))
			Prev_m <- ifelse(apply(out[mday,5:7],1,sum)<1, 5e-3, out[mday,6])/ifelse(apply(out[mday,5:7],1,sum)<1, 1, apply(out[mday,5:7],1,sum))
		}
		
		-2*sum(dbinom(pos_r, n_r, Prev_r[c(3:12,1:2)], log = T),
		   	   dbinom(pos_m, n_m, Prev_m[c(3:12,1:2)], log = T),
		   	   dbinom(pos_u, n_u, Prev_a[c(3:12,1:2)], log = T), na.rm = T)     
	}
	
	
	start <- Mod[i, which(unlist(lapply(strsplit(names(Mod), "_"), function(x) any(x=="s"))))[parId]]
			 names(start) <- apply(cbind(names(start)), 1, function(x) substr(x, 1, nchar(x)-2))
	upper <- Mod[i, which(unlist(lapply(strsplit(names(Mod), "_"), function(x) any(x=="u"))))[parId]]
			 names(upper) <- names(start)
	lower <- Mod[i, which(unlist(lapply(strsplit(names(Mod), "_"), function(x) any(x=="l"))))[parId]]
			 names(lower) <- names(start)


	mcmc <- modMCMC(f = mod_cost_nll, p = as.numeric(start), upper = as.numeric(upper), lower = as.numeric(lower),
					niter = 15000, updatecov = 100, ntrydr = 5, wvar0 = NULL, outputlength = 5000)

	save(mcmc, file = paste("Data/SimResults/mcmc_", paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = "") ,".RData", sep = ""))
	
	out[i,2] <- mcmc$naccepted
	out[i,3] <- mod_cost_nll(mcmc$bestpar)
	
	## Sensitivity analaysis
	parRange <- cbind(min = as.numeric(lower), 
    	              max = as.numeric(upper))
	rownames(parRange) <- names(start)
	
	sens <- modCRL(fun = mod_cost_nll, parRange = parRange, num = 5000)
	save(sens, file = paste("Data/SimResults/sens_", paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = ""),".RData", sep = ""))

	### Confidence of fit
	best_sir <- function(lpars) {
		pars <- rep(NA, 15)
		pars[parId]  <- lpars
		out <- sirs(p = pars, mod = paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = ""))
	
		mday  <- which(as.numeric(format(as.POSIXlt(as.POSIXct("2012-01-01 00:00:00")+((out[,1]-1)*24*60*60)), "%d"))==15)

	
		if(ncol(out)==10) {
			cbind(1:12, apply(out[mday,c(3,6,9)], 1, sum)/apply(out[mday, 2:10],1,sum))
			} else {
			cbind(1:12, apply(out[mday,c(3,6)], 1, sum)/apply(out[mday, 2:7],1,sum))
			}	
	}		

	sR <- sensRange(func = best_sir, parms = NULL, parInput = mcmc$par)
	save(sR, file = paste("Data/SimResults/sR_", paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = ""),".RData", sep = ""))
	
	dyn.unload(paste("C_Functions/", paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = ""), ".so", sep = ""))
	write.csv(out, "Data/SimResults/Results2.csv", row.names = F)
	}
}


for(i in 1:nrow(Mod)) {
	
	dyn.load(paste("C_Functions/", paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = ""), ".so", sep = ""))
	
	parId <- as.numeric(unlist(strsplit(as.character(Mod$parId[i]), ".", fixed = T)))
	
	mod <- paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = "")
	
load(paste("Data/SimResults/mcmc_", mod ,".RData", sep = ""))
load(paste("Data/SimResults/sens_", mod ,".RData", sep = ""))
load(paste("Data/SimResults/sR_", mod ,".RData", sep = ""))


    pars <- rep(NA, 14)
  	pars[parId]  <- summary(mcmc)["mean",]
		out <- sirs(p = pars, mod = mod)

save(out, file = paste("Data/SimResults/out_", mod,".RData", sep = ""))
dyn.unload(paste("C_Functions/", paste(Mod[i,which(names(Mod)%in%letters[1:6])], collapse = ""), ".so", sep = ""))
}



#### unload .c functions
dyn.unload("C_Functions/dem.so")
########################