##############################
# These are the simultaions for COVID-19
##############################

library(deSolve)



mitigation <- data.frame(matrix(NA,500,7))
for(i in 1:500){
mitigation[i,1] <- 1 - i*0.001
}

colnames(mitigation) <- c("contact","Infected_all","Deceased_all","Infected_lowOnly","Deceased_lowOnly","Infected_highOnly","Deceased_highOnly")



####################
#### Mitigation impact if we all contribute!
####################

N <- 1000000

parameters <- function(case) {
    for(i in 1:500){
    if(case==i) {
        t_I_H <- 1 - i*0.001
        t_I_L <- 1 - i*0.001
        t_S_H <- 1 - i*0.001
        t_S_L <- 1 - i*0.001
}
}
pars <- c(beta = 0.146,
            theta_I_H = t_I_H,
            theta_I_L = t_I_L,
            theta_S_H = t_S_H,
            theta_S_L = t_S_L,
            sigma_H = 0.192,      #1/5.2
            sigma_L = 0.192,
            gamma_H = 0.0583,     #1/17.16 ceiling or an average of 1/12.39=0.0807
            gamma_L = 0.0583,
            mu_H = 0.002,         #0.0004/.2
            mu_L = 0.0004)
  return(pars)
}

state <- c(S_H = 0.2*N,
           S_L = 0.8*N,
           E_H = 0,
           E_L = 100,
           I_H = 0,
           I_L = 0,
           R_H = 0,
           R_L = 0,
           D_H = 0,
           D_L = 0)


seir_2_risk_covid <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dS_H <- (-beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_H * S_H))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L)
    dS_L <- (-beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_L * S_L))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L)
    dE_H <- (beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_H * S_H))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L) - (sigma_H * E_H)
    dE_L <- (beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_L * S_L))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L) - (sigma_L * E_L)
    dI_H <- (sigma_H * E_H) - (gamma_H + mu_H)*I_H
    dI_L <- (sigma_L * E_L) - (gamma_L + mu_L)*I_L
    dR_H <- (gamma_H * I_H)
    dR_L <- (gamma_L * I_L)
    dD_H <- (mu_H * I_H)
    dD_L <- (mu_L * I_L)
    # return the rate of change
    list(c(dS_H, dS_L, dE_H, dE_L, dI_H, dI_L, dR_H, dR_L, dD_H, dD_L))
  }) # end with(as.list ...
}

simulate <- function(case,total.time) {
  times <- seq(0, total.time, by = 1)
  out <- ode(y = state, times = times, func = seir_2_risk_covid, parms = parameters(case))
  colnames(out) <- c("time","Susceptible_H","Susceptible_L",
                       "Exposed_H","Exposed_L","Infected_H","Infected_L",
                       "Recovered_H","Recovered_L","Deceased_H","Deceased_L")
  epi.df<-as.data.frame(out)
  return(epi.df)
}


for(i in 1:500) {
    epi<-simulate(i,3000)
    mitigation[i,2] <- (100 - ((epi[901,2]+epi[901,3])/N*100))
    mitigation[i,3] <- (epi[901,10]+epi[901,11])
}



####################
#### Mitigation impact if we only low risk contribute !
####################

N <- 1000000

parameters <- function(case) {
    for(i in 1:500){
    if(case==i) {
        t_I_H <- 1
        t_I_L <- 1 - i*0.001
        t_S_H <- 1
        t_S_L <- 1 - i*0.001
}
}
pars <- c(beta = 0.146,
            theta_I_H = t_I_H,
            theta_I_L = t_I_L,
            theta_S_H = t_S_H,
            theta_S_L = t_S_L,
            sigma_H = 0.192,      #1/5.2
            sigma_L = 0.192,
            gamma_H = 0.0583,     #1/17.16 ceiling or an average of 1/12.39=0.0807
            gamma_L = 0.0583,
            mu_H = 0.002,         #0.0004/.2
            mu_L = 0.0004)
  return(pars)
}

state <- c(S_H = 0.2*N,
           S_L = 0.8*N,
           E_H = 0,
           E_L = 100,
           I_H = 0,
           I_L = 0,
           R_H = 0,
           R_L = 0,
           D_H = 0,
           D_L = 0)


seir_2_risk_covid <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dS_H <- (-beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_H * S_H))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L)
    dS_L <- (-beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_L * S_L))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L)
    dE_H <- (beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_H * S_H))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L) - (sigma_H * E_H)
    dE_L <- (beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_L * S_L))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L) - (sigma_L * E_L)
    dI_H <- (sigma_H * E_H) - (gamma_H + mu_H)*I_H
    dI_L <- (sigma_L * E_L) - (gamma_L + mu_L)*I_L
    dR_H <- (gamma_H * I_H)
    dR_L <- (gamma_L * I_L)
    dD_H <- (mu_H * I_H)
    dD_L <- (mu_L * I_L)
    # return the rate of change
    list(c(dS_H, dS_L, dE_H, dE_L, dI_H, dI_L, dR_H, dR_L, dD_H, dD_L))
  }) # end with(as.list ...
}

simulate <- function(case,total.time) {
  times <- seq(0, total.time, by = 1)
  out <- ode(y = state, times = times, func = seir_2_risk_covid, parms = parameters(case))
  colnames(out) <- c("time","Susceptible_H","Susceptible_L",
                       "Exposed_H","Exposed_L","Infected_H","Infected_L",
                       "Recovered_H","Recovered_L","Deceased_H","Deceased_L")
  epi.df<-as.data.frame(out)
  return(epi.df)
}


for(i in 1:500) {
    epi<-simulate(i,3000)
    mitigation[i,4] <- (100 - ((epi[901,2]+epi[901,3])/N*100))
    mitigation[i,5] <- (epi[901,10]+epi[901,11])
}




####################
#### Mitigation impact if we only high risk contribute !
####################

N <- 1000000

parameters <- function(case) {
    for(i in 1:500){
    if(case==i) {
        t_I_H <- 1 - i*0.001
        t_I_L <- 1
        t_S_H <- 1 - i*0.001
        t_S_L <- 1
}
}
pars <- c(beta = 0.146,
            theta_I_H = t_I_H,
            theta_I_L = t_I_L,
            theta_S_H = t_S_H,
            theta_S_L = t_S_L,
            sigma_H = 0.192,      #1/5.2
            sigma_L = 0.192,
            gamma_H = 0.0583,     #1/17.16 ceiling or an average of 1/12.39=0.0807
            gamma_L = 0.0583,
            mu_H = 0.002,         #0.0004/.2
            mu_L = 0.0004)
  return(pars)
}

state <- c(S_H = 0.2*N,
           S_L = 0.8*N,
           E_H = 0,
           E_L = 100,
           I_H = 0,
           I_L = 0,
           R_H = 0,
           R_L = 0,
           D_H = 0,
           D_L = 0)


seir_2_risk_covid <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dS_H <- (-beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_H * S_H))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L)
    dS_L <- (-beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_L * S_L))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L)
    dE_H <- (beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_H * S_H))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L) - (sigma_H * E_H)
    dE_L <- (beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_L * S_L))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L) - (sigma_L * E_L)
    dI_H <- (sigma_H * E_H) - (gamma_H + mu_H)*I_H
    dI_L <- (sigma_L * E_L) - (gamma_L + mu_L)*I_L
    dR_H <- (gamma_H * I_H)
    dR_L <- (gamma_L * I_L)
    dD_H <- (mu_H * I_H)
    dD_L <- (mu_L * I_L)
    # return the rate of change
    list(c(dS_H, dS_L, dE_H, dE_L, dI_H, dI_L, dR_H, dR_L, dD_H, dD_L))
  }) # end with(as.list ...
}

simulate <- function(case,total.time) {
  times <- seq(0, total.time, by = 1)
  out <- ode(y = state, times = times, func = seir_2_risk_covid, parms = parameters(case))
  colnames(out) <- c("time","Susceptible_H","Susceptible_L",
                       "Exposed_H","Exposed_L","Infected_H","Infected_L",
                       "Recovered_H","Recovered_L","Deceased_H","Deceased_L")
  epi.df<-as.data.frame(out)
  return(epi.df)
}


for(i in 1:500) {
    epi<-simulate(i,3000)
    mitigation[i,6] <- (100 - ((epi[901,2]+epi[901,3])/N*100))
    mitigation[i,7] <- (epi[901,10]+epi[901,11])
}




#########################
## Plotting
#########################

plot(-mitigation$contact, mitigation$Infected_all, type="l", lwd=2, xlab = "degree of proactive measures", ylab = "Strengh of epidemic", col=rgb(1,0,1),ylim=c(0,100), xaxt="none", yaxt="none")
lines(-mitigation$contact, mitigation$Infected_lowOnly, lwd=2, col="blue")
lines(-mitigation$contact, mitigation$Infected_highOnly, lwd=2, col="red")
#axis(side=1, at=c(-0.5,-0.6,-0.7,-0.8,-0.9,-1), labels=c("0.5*beta","0.6*beta","0.7*beta","0.8*beta","0.9*beta","beta"))
axis(side=1, at=c(-0.5,-1), labels=c("half transmission","full transmission"))
axis(side=2, at=c(20,85), labels=c("weak","strong"), las=2)
mtext(outer = TRUE, side = 3, "Impact of mitigation epidemic", cex = 1.5)
text(-0.9,65,"We all contribute",col=rgb(1,0,1),cex=0.7,pos=4)
text(-0.8,73,"only low risk groups contribute",cex=0.7,col="blue",pos=4)
text(-0.85,89,"only high risk groups contribute",cex=0.7,col="red",pos=4)




plot(-mitigation$contact, mitigation$Deceased_all/10, type="l", lwd=2, xlab = "degree of proactive measures", ylab = "Expected number of deaths per 100,000 individuals", col=rgb(1,0,1), xaxt="none")
lines(-mitigation$contact, mitigation$Deceased_lowOnly/10, lwd=2, col="blue")
lines(-mitigation$contact, mitigation$Deceased_highOnly/10, lwd=2, col="red")
axis(side=1, at=c(-0.5,-1), labels=c("half transmission","full transmission"))
 mtext(outer = TRUE, side = 3, "Impact of mitigation epidemic", cex = 1.5)
text(-0.96,900,"We all contribute",col=rgb(1,0,1),cex=0.7,pos=4)
text(-0.75,830,"only low risk groups contribute",cex=0.7,col="blue",pos=4)
text(-0.85,1045,"only high risk groups contribute",cex=0.7,col="red",pos=4)



#############################################
########################
### This is a comparison of a reduction of to 0.9, 0.8 and 0.5 of the original value of beta
########################
############################################

N <- 1000000

parameters <- function(case) {
  if(case==1) {
    t_I_H <- 1.0
    t_I_L <- 1.0
    t_S_H <- 1.0
    t_S_L <- 1.0
  }
  if(case==2) {
    t_I_H <- .9
    t_I_L <- .9
    t_S_H <- .9
    t_S_L <- .9
  }
  if(case==3) {
    t_I_H <- .8
    t_I_L <- .8
    t_S_H <- .8
    t_S_L <- .8
  }
  if(case==4) {
    t_I_H <- .65
    t_I_L <- .65
    t_S_H <- .65
    t_S_L <- .65
  }
  if(case==5) {
    t_I_H <- .55
    t_I_L <- .55
    t_S_H <- .55
    t_S_L <- .55
  }
  pars <- c(beta = 0.146,
            theta_I_H = t_I_H,
            theta_I_L = t_I_L,
            theta_S_H = t_S_H,
            theta_S_L = t_S_L,
            sigma_H = 0.192,      #1/5.2
            sigma_L = 0.192,
            gamma_H = 0.0583,     #1/17.16
            gamma_L = 0.0583,
            mu_H = 0.002,         #0.0004/.2
            mu_L = 0.0004)
  return(pars)
}

state <- c(S_H = 0.2*N,
           S_L = 0.8*N,
           E_H = 0,
           E_L = 0,
           I_H = 100,
           I_L = 0,
           R_H = 0,
           R_L = 0,
           D_H = 0,
           D_L = 0)

seir_2_risk_covid <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dS_H <- (-beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_H * S_H))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L)
    dS_L <- (-beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_L * S_L))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L)
    dE_H <- (beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_H * S_H))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L) - (sigma_H * E_H)
    dE_L <- (beta*(theta_I_H * I_H + theta_I_L * I_L)*(theta_S_L * S_L))/(S_H + S_L + E_H + E_L + I_H + I_L + R_H + R_L) - (sigma_L * E_L)
    dI_H <- (sigma_H * E_H) - (gamma_H + mu_H)*I_H
    dI_L <- (sigma_L * E_L) - (gamma_L + mu_L)*I_L
    dR_H <- (gamma_H * I_H)
    dR_L <- (gamma_L * I_L)
    dD_H <- (mu_H * I_H)
    dD_L <- (mu_L * I_L)
    # return the rate of change
    list(c(dS_H, dS_L, dE_H, dE_L, dI_H, dI_L, dR_H, dR_L, dD_H, dD_L))
  }) # end with(as.list ...
}


simulate <- function(case,total.time) {
  times <- seq(0, total.time, by = 1)
  out <- ode(y = state, times = times, func = seir_2_risk_covid, parms = parameters(case))
  colnames(out) <- c("time","Susceptible_H","Susceptible_L",
                       "Exposed_H","Exposed_L","Infected_H","Infected_L",
                       "Recovered_H","Recovered_L","Deceased_H","Deceased_L")
  epi.df<-as.data.frame(out)
  return(epi.df)
}


add.curve <- function(x, y, first, rgbv, ylabel, ylimit) {
  if(first) {
    par(mar = c(4, 4, .5, .5),mfrow=c(1,1))
    plot(x, y, type="l", lwd=2, xlab = "time (days)", ylab = ylabel,
         col=rgb(rgbv[1],rgbv[2],rgbv[3]),ylim=c(ylimit[1],ylimit[2]))
         #polygon(c(x,rev(x)),c(y,rep(0,length(x))),
         #col=rgb(rgbv[1],rgbv[2],rgbv[3],0.5),border=NA)
    lines(x, y, lwd=2, col=rgb(rgbv[1],rgbv[2],rgbv[3]))
  } else {
      #polygon(c(x,rev(x)),c(y,rep(0,length(x))),
      #col=rgb(rgbv[1],rgbv[2],rgbv[3],0.5),border=NA)
    lines(x, y, lwd=2, col=rgb(rgbv[1],rgbv[2],rgbv[3]))
  }
}

color.vec <- function(case) {
  if(case==1) {
    return(c(0,0,0))
  }
  if(case==2) {
    return(c(1,0,0))
  }
  if(case==3) {
    return(c(0,0,1))
  }
  if(case==4) {
    return(c(1,0,1))
  }
  if(case==5) {
      return(c(1,.5,1))
   }
}

#Draw total infected cases plot
for(i in 1:5) {
  epi<-simulate(i,4500)
  add.curve(epi$time, log(epi$Infected_H + epi$Infected_L), i==1,
            color.vec(i), "infected cases", c(0,log(175000)))
}
text(150,12,"no mitigation",cex=0.7,pos=4)
text(200,11.5,"0.9*beta",cex=0.7,col="red",pos=4)
text(350,10.8,"0.8*beta",cex=0.7,col="blue",pos=4)
text(2500,6.2,"0.65*beta",cex=0.7,col=rgb(1,0,1),pos=4)
text(30,4.2,"0.55*beta",cex=0.7,col=rgb(1,0.5,1),pos=4)


#Draw high-risk infected cases plot
for(i in 1:5) {
  epi<-simulate(i,1500)
  add.curve(epi$time, log(epi$Infected_H), i==1,
            color.vec(i), "high-risk infected cases", c(0,log(40000)))
}
text(180,33000,"no protective measures",cex=0.7,pos=4)
text(200,28000,"protective measures only for high-risk individuals",cex=0.7,col="red",pos=4)
text(230,19600,"protective measures only for low-risk individuals",cex=0.7,col="blue",pos=4)
text(360,11750,"full protective measures",cex=0.7,col=rgb(1,0,1),pos=4)


#Draw deceased plot
for(i in 1:4) {
  epi<-simulate(i,500)
  add.curve(epi$time, epi$Deceased_H, i==1,
            color.vec(i), "cumulative deaths", c(0,6500))
}
