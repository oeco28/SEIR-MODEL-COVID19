##############################
# These are the simulations for COVID-19
#############################

install.packages("deSolve")
library(deSolve)


#############################################
## Simple model 
############################################


N=1000000

parameters <- c(b = 0.186,   # or 146 for the 17.16 case. # 0.186 for the 12.39
	s = 0.192,      #1/5.2
	g = 0.0807,     #average 1/12.39=0.0807  or ceiling 1/17.16=0.0583
	m = 0.0004)

state <- c(X = 1000000,
	Y = 1,
	Z = 0,
	Z0 = 0,
	W = 0)

seir_covid <-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		dX <- (-b*X*Z)/(X+Y+Z+W)
		dY <- (b*X*Z)/(X+Y+Z+W) - s*Y
		dZ <- s*Y - (g + m)*Z
		dZ0 <- m*Z
		dW <- g*Z
		# return the rate of change
		list(c(dX, dY, dZ, dZ0, dW))
	}) # end with(as.list ...
}



times <- seq(0, 120, by = 1)


out_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
head(out_1)

par(oma = c(0, 0, 4, 0))
plot(out_1, xlab = "time (days)", ylab = "Number of individuals")
plot(out_1[, "Susceptible"], out_1[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "SEIR 120 days post introduction", cex = 1.5)

out_1[time=121,]/N*100                 # Percentages of each class at a year after introduction
(100 - (out_1[time=121,]/N*100)[2])    # The Percent of cumulative number of infected individuals at the end of the epidemic
(1 - (out_1[time=121,]/N)[2])*N        # cumulative number of infected individuals at the end of the epidemic



times <- seq(0, 365, by = 1)


out_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
head(out_1)

par(oma = c(0, 0, 4, 0))
plot(out_1, xlab = "time (days)", ylab = "Number of individuals")
plot(out_1[, "Susceptible"], out_1[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "SEIR 365 days post introduction", cex = 1.5)

out_1[time=366,]/N*100                 # Percentages of each class at a year after introduction
(100 - (out_1[time=366,]/N*100)[2])    # The Percent of cumulative number of infected individuals at the end of the epidemic
(1 - (out_1[time=366,]/N)[2])*N        # cumulative number of infected individuals at the end of the epidemic


out_1.1 <- as.data.frame(out_1)
out_1.1[which(out_1.1[,"Infected"] == max(out_1.1[,"Infected"])),]   # Number of infected individuals at the peak of the epidemic! This is important because it will tell us something about the burden on the health care system




#####################################
# transmission reduction efforts 2/3*b after 30 days of introduction for case E=1 We assume that the number of individuals in each class is that observed from the previous dynamic with no reduction in transmission
####################################



parameters <- c(b = (2/3)*0.186,   # or 146 for the 17.16 case. # 0.186 for the 12.39
	s = 0.192,      #1/5.2
	g = 0.0807,     #average 1/12.39=0.0807  or ceiling 1/17.16=0.0583
	m = 0.0004)


state <- c(X = 9.999947e+05,
	Y = 1.546402e+00,
	Z = 2.614416e+0,
	Z0 = 1.448262e-02,
	W = 2.110842e+00)

seir_covid <-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		dX <- (-b*X*Z)/(X+Y+Z+W)
		dY <- (b*X*Z)/(X+Y+Z+W) - s*Y
		dZ <- s*Y - (g + m)*Z
		dZ0 <- m*Z
		dW <- g*Z
		# return the rate of change
		list(c(dX, dY, dZ, dZ0, dW))
	}) # end with(as.list ...
}



times <- seq(0, 90, by = 1)

out_m_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_m_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
head(out_m_1)

par(oma = c(0, 0, 4, 0))
plot(out_m_1, xlab = "time (days after mitigation)", ylab = "Number of individuals")
plot(out_m_1[, "Susceptible"], out_m_1[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "mitigation early after outbreak", cex = 1.5)


out_m_1[time=91,]/N*100
(100 - (out_m_1[time=91,]/N*100)[2]) 
(1 - (out_m_1[time=91,]/N)[2])*N 
 


times <- seq(0, 630, by = 1)

out_m_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_m_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
head(out_m_1)

par(oma = c(0, 0, 4, 0))
plot(out_m_1, xlab = "time (days after mitigation)", ylab = "Number of individuals")
plot(out_m_1[, "Susceptible"], out_m_1[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "mitigation full epidemic", cex = 1.5)


out_m_1[time=631,]/N*100
(100 - (out_m_1[time=631,]/N*100)[2]) 
(1 - (out_m_1[time=631,]/N)[2])*N 
 

out_m_1.1 <- as.data.frame(out_m_1)
out_m_1.1[which(out_m_1.1[,"Infected"] == max(out_m_1.1[,"Infected"])),]




######################
# comparison with and without mitigation
#####################


par(mar = c(4, 4, .5, .5),mfrow=c(1,1))
plot(out_1.1$time, out_1.1$Infected, type="l", lwd=2, xlab = "time (days)", ylab = "Infected individuals", 
    col=rgb(0,0,0),ylim=c(0,(out_1.1[which(out_1.1[,"Infected"] == max(out_1.1[,"Infected"])),4] + 100)), xlim=c(0,650))
polygon(c(out_1.1$time,rev(out_1.1$time)),c(out_1.1$Infected,rep(0,length(out_1.1$time))),
    col=rgb(0,0,0,0.5),border=NA)
lines(out_m_1.1$time, out_m_1.1$Infected, lwd=2, col=rgb(0,0,1))
polygon(c(out_m_1.1$time,rev(out_m_1.1$time)),c(out_m_1.1$Infected,rep(0,length(out_m_1.1$time))),
   col=rgb(0,0,1,0.5),border=NA)
    

text(250,100000,"no protective measures",cex=0.7,pos=4)
text(430,45000,"mitigation to 2/3 of transmission rate",cex=0.7,col=rgb(0,0,1),pos=4)


