# Susceptible, Exposed, Infected and Removed (SEIR) model for Coronavirus spread. The Pullman introduction case

#### by: Omar E. Cornejo, PhD (School of Biological Sciences). March 11, 2020

Most models I have seen to this date focusing on SARS-CoV-2 are of the SEIR type. Something to be recognized is that there is a lot of variability in the estimates of $R_0$ ranging from 2 to 3.5 Most of the estimates seem to converge towards $R_0=2.3$.  So I am going to assume that this is a reasonable estimate.  A modified SEIR model from influenza has been used to estimate the number of cases progression in Wuhan$^1$:
$$
\begin{align}
\dot{S}=-\beta SI/N \\
\dot{E}=\beta SI/N-\sigma E \\
\dot{I}=\sigma E- \gamma I \\
\dot{R}=\gamma I \\
N=S+E+I+R \\
\end{align}
$$
where $\beta$ is the transmission rate, $\sigma$ is the infection rate (development of infection) and it is usually estimated as the inverse of the incubation period, and $\gamma$ is the recovery rate estimated as the inverse of the infectious period.

For this we need to make reasonable assumptions about incubation period and the infectious period. It has been assumed that $\sigma=1/5.2$ under the estimation that the incubation period is 5.2 days$^2$. In Wuhan, it was estimated that hospitalization period of recovered patients was 12.39$\pm$ 4.77 days$^{3,4}$. For simplicity we assume the average of the infectious period to be then 12.39 and $\gamma=1/12.39=0.0807$.  I had originally considered the ceiling of the recovery time, but given the general health status of our population, it seems reasonable to make the estimations with the average.

$R_0$ is usually estimated as the leading eigenvalue from the $FV^{-1}$ matrix that describes the terms for the acquisition of the infections and the transition between states. In the simplest SEIR model that does not consider the age structure we can then simply describe it as:

$R_0=\frac{\beta\sigma}{\sigma\gamma}$ 

The transmission rate then can be estimated from the relationship $\beta=R_0\gamma$ from this model. Given that the mean $R_0$ has been estimated to be 2.3 across multiple studies$^{1,2,5,6}$, we can make a conservative estimate of $\beta$ assuming that recovery takes the ceiling of the distribution (17.16 days), that gives us $\beta=2.3/17.16=0.147$.

All scenarios were examined using simulations in R 3.6.4 using the *deSolve* package to explore numerically the dynamics. 

To set the starting point for the simulations we assume that the Susceptible population of Pullman follows the census data available at [https://datausa.io/profile/geo/pullman-wa](https://datausa.io/profile/geo/pullman-wa) (Figure 1). We will assume that the total number of individuals is 32,382. We perform simulations under two different introduction scenarios, all assuming that exposed individuals travel from the west side to Pullman. The two scenarios assume: i) that 1 exposed individual ($E_{b}=1$) introduces the disease to town;and  ii) that 10 exposed individuals ($E_{b}=10$) introduce the disease to town. (I considered the third scenario that 100 exposed individuals, $E_{b}=100$, introduce the disease to town, but it seem to bad).

![census data](https://tva1.sinaimg.cn/large/00831rSTgy1gcpwcmny5cj30s00b9t97.jpg)

Figure 1. Census data from the town of Pullman

## Wang model (assuming no mortality and no age population structure)

### Scenario 1: 1 exposed (E) individual arrives to Pullman

Under the assumption that 1 exposed individual arrives to Pullman, we examine the expected dynamics after 120 days and 420 days after introduction. Figure 2 shows that after 120 days after the introduction of a single infected individual we would expect that approximately 118 individuals would have been exposed (0.36% of the Pullman population). Out of these,  we would expect to see 31 infected individuals (0.1% of the Pullman population) and 74 recovered individuals (0.22% of the Pullman population) at day 120 post introduction.

![initial_spread](https://tva1.sinaimg.cn/large/00831rSTgy1gcsy5kzjlfj30xp0u0n3d.jpg)

Figure 2. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals during the intial phase of the outbreak caused by 1 exposed individual.

If we allow the simulations to continue for 420 days (14 months) we predict that at the peak of the epidemic we will have 1,068 infected individuals (at 236 days after the introduction of the exposed individual), giving you an idea of the potential burden of the number of disease individuals to our hospital resources in town. Now, this could be an unreasonable situation because even in Wuhan, it has been observed that as the epidemic progresses, $R_0$ has been estimated to decay. The multi-phase analysis of the epidemics has suggested that in these places $R_0=1.3$ after the initial introduction$^1$. It is unclear to me at the moment if this decay in $R_0$ is due to mitigation effects of public health measurements or changes in the intrinsic dynamic of the epidemic. In this scenario, after 14 months of transmission 21,522 individuals ($66.5\%$) of the population of Pullman could have been infected.

![prolongued_nonmitigated_epidemic](https://tva1.sinaimg.cn/large/00831rSTgy1gcsy5urm7mj30xp0u0451.jpg)

Figure 3. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals for 10 months after introduction after 1 Exposed individual.

### Scenario 2: 10 exposed (E) individuals arrive to Pullman

Under the assumption that 10 exposed individuals arrive to Pullman after Spring Break, we examine the expected dynamics after 120 days and 300 days after introduction. Figure 4 shows that after 120 days after the introduction of ten  infected individuals we would expect that approximately 1,165 individuals would have been exposed and or infected (3.6% of the Pullman population). Out of these,  we would expect to see 291 infected individuals (0.9 % of the Pullman population) and 715 recovered individuals (2.21% of the Pullman population) at day 120 post introduction.

![10_initial_spread](https://tva1.sinaimg.cn/large/00831rSTgy1gcsym7rpnoj30xp0u0age.jpg)Figure 4. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals during the intial phase of the outbreak caused by 10 exposed individuals.

If we allow the simulations to continue for 420 days (14 months) we predict that at the peak of the epidemic we will have 2,106 infected individuals (at 221 days after the introduction of the exposed individual, 15 days earlier than that observed for the case of a single introduction). In this scenario, after 14 months of transmission 21,762 individuals ($67.2\%$) of the population of Pullman could have been infected (Figure 5).

![10_prolongued](https://tva1.sinaimg.cn/large/00831rSTgy1gcsymi0zouj30xp0u0jxz.jpg)Figure 5. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals for 10 months after introduction after 10 Exposed individuals.

## Considerations about the mortality rate associated to SARS-CoV-2 infections

At this point the mortality rate by SARS-CoV-2 seems to be remarkably different between individuals of different ages$^7$.  In this study, following 72,314 cases in Wuhan 44,672 (62%) were classified as confirmed cases of COVID-19, the number of suspected cases was 16 186 (22%), the number of diagnosed cases 10,567 (15%) and the number of asymptomatic 889 (1%). I am going to assume the distribution of fatality rate from:c

  1. 2.3% (1023 of 44 672 confirmed cases)
  2. 14.8% in patients aged ≥80 years (208 of 1408)
  3. 8.0% in patients aged 70-79 years (312 of 3918)
  4. 49.0% in critical cases (1023 of 2087)

From these, I am going to assume the mortality rate ($\mu_{o}$) of individuals > 70 year old to be 8% or 11.4% in a 2 months period to model the spread of disease in an age structured population for individuals, and consider the mortality rate ($\mu_{b}$) for the rest of the population to be 2.3% in a 2 months period in which these observations were done (60 days). The mortality rate for >70 yo would then be $\mu_{o}=0.0013 $ *deaths/day* or $\mu_{o}=0.0019$ *deaths/day* and the mortality rate of rest of the population is $\mu_{b}=0.0004$.

### SEIR model including mortality rate as part of the dynamics

I have already established in the previous section what the mortality rate would be for individuals younger than 70 years old,  $\mu_{b}=0.0004$. This is assuming that 2.3% of the cases result in fatalities as has been estimated from the longest cohort to the moment $^7$. We are going to assume this as an overall mortality rate under the assumption that the population at risk is only a small fraction of the population. We will extend the model later to consider the population at risk separately

Our back of the envelop calculations make a very pessimistic prediction about the number of fatalities. I intend to determine in this part of the analysis what is the impact of including the mortality rate in the dynamics of the epidemic. Under this assumptions I extended the model from Wang et al.$^1$ to include mortality
$$
\begin{align}\dot{S}=-\beta SI/N \\\dot{E}=\beta SI/N-\sigma E \\\dot{I}=\sigma E-(\gamma+\mu)I \\\dot{R}=\gamma I \\\end{align}
$$
In this case we performed the analysis considering again two possible scenarios: i) the introduction of one infected individual into Pullman ($E=1$); and ii) the introduction of 10 infected individuals ($E=10$).

In this case:

$R_0=\frac{\beta\sigma}{\sigma(\gamma + \mu)}$ 

### Scenario 1 with mortality explicitely included in the dynamics: 1 exposed (E) individual arrives to Pullman

Under the assumption that 1 exposed individual arrives to Pullman, we examine the expected dynamics after 120 days and 420 days after introduction. Figure 6 shows that after 120 days after the introduction of a single infected individual we would expect that approximately 119 individuals would have been exposed/infected (0.37% of the Pullman population). Out of these,  we would expect to see 30 infected individuals (0.09% of the Pullman population) and 72 recovered individuals (0.22% of the Pullman population), and 0 deaths at day 120 post introduction.  

![1_introduction](https://tva1.sinaimg.cn/large/00831rSTgy1gcszcvmnqgj30xp0u0n47.jpg)Figure 6. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals during the intial phase of the outbreak caused by 1 exposed individual assuming an overall mortality rate of 2.3% over a two month period ($\mu=0.0004$ *per day*).

If we allow the simulations to continue for 420 days (ten months) we predict that at the peak of the epidemic we will have 2,062 infected individuals (at 292 days after the introduction of the exposed individual), giving you an idea of the potential burden of the number of disease individuals to our hospital resources in town. The mortality rate does not change this estimate significantly.  In this scenario, after 14 months of transmission 21,373 individuals ($66\%$) of the population of Pullman could have been infected. The number of deaths at the end of the epidemic is predicted to be 104 individuals (See Figure 7).

![1_prolongued](https://tva1.sinaimg.cn/large/00831rSTgy1gcszmo4aarj30xp0u0ahc.jpg)Figure 7. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals for 10 months after introduction after 1 Exposed individuals and a mortality rate $\mu=0.0004$. As you can see the number of deaths at the end of the epidemic could reach 100 individuals or more.

### Scenario 2 with mortality explicitely included in the dynamics: 10 exposed (E) individuals arrive to Pullman

Under the assumption that 10 exposed individual arrives to Pullman, we examine the expected dynamics after 120 days and 420 days after introduction. Figure 8 shows that after 120 days after the introduction of ten infected individuals we would expect that approximately 1,139 individuals would have been exposed (3.5% of the Pullman population). Out of these,  we would expect to see 283 infected individuals (0.87% of the Pullman population), 698 recovered individuals (2.15% of the Pullman population), and 3-4 deaths (0.01% of the Pullman population) at day 120 post introduction.  Yet, you can see the effec that this can have on the mortality at 120 days post introduction.

![10_introduction](https://tva1.sinaimg.cn/large/00831rSTgy1gcsztaloosj30xp0u045l.jpg)Figure 8. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals during the intial phase of the outbreak caused by 10 exposed individual assuming an overall mortality rate of 2.3% over a two month period ($\mu=0.0004$ *per day*).

If we allow the simulations to continue for 420 days (ten months) we predict that at the peak of the epidemic we will have 2,069 infected individuals (at 222 days after the introduction of the exposed individuals which is 70 days earlier than what estimated for $E=1$), giving you an idea of the potential burden of the number of disease individuals to our hospital resources in town (Figure 9). The mortality rate does not change this estimate by a large amount.  In this scenario, after 14 months of transmission 21,631 individuals ($66.8\%$) of the population of Pullman could have been infected. The number of deaths at the end of the epidemic is predicted to be 107 individuals.

![10_prolongued](https://tva1.sinaimg.cn/large/00831rSTgy1gcszxgtv0lj30xp0u0wlv.jpg)Figure 9. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals for 10 months after introduction after 10 Exposed individuals and a mortality rate $\mu=0.0004$. As you can see the number of deaths at the end of the epidemic could reach 100 individuals or more.

#### *My overall warning from these simulations is if we fail to prevent the introduction of multiple cases early in the epidemic, we will have a larger number of potential deaths early on. The overall number in the long run of the epidemic after mitigation might not be that different, something that is yet to be determined after models that incorporate the effect of mitigation are appropriately addressed.*

## Let us consider mitigation strategies. What about if we implement measures that reduce the rate of contact and thus transmission

In this case I am goint to examine the same two scenarios that consider the introduction of $E=1$ exposed individual or $E=10$ exposed individuals under the model in which mortality is explicitely taken into account in the dynamics. For this, I would like to look at an intervention strategies that would result in the reduction of the transmission rate to 2/3 of the original value ($\beta_{m}=2/3\beta$).  Because of time constraints I will leave other scenarios for another update of the document.

I am going to assume that intervention does not occur but after 30 days of the introduction of the exposed individual(s) and will examine the consequences of the mitigation strategy in the first 200 days after mitigation.  For this, I will use as starting conditions for the simulation the numbers of Susceptible, Infected, Exposed and Recovered that were obtained from the previous dynamics with $\beta=0.146$.  My intention in this case is to examine not just how mitigation will affect the dynamics of the disease but also the final numbers of infected, recovered and dead individuals in the time frames we have examined all previous scenarios.

The starting conditions, 30 days after introduction, under the scenario where the disease is introduced by a single individual are:  $S_{30}=32,345$, $E_{30}=0.91$, $I_{30}=1.53$, $R_{30}=2.2$, and $deads=1.09e-02$.

The starting conditions, 30 days after introduction, under the scenario where the disease is introduced by a single individual are:  $S_{30}=32,345$, $E_{30}=9.09$, $I_{30}=15.29.6$, $R_{30}=21.97$, and $deads=1.09e-01$. 

### Reduction of transmission rate to $2/3\beta$ 

In this case the underlying model is the same as the one described by equations 6-9 but we assume a reduced rate of spread.

Our simulations indicate that reducing the transmission rate to $2/3\beta$ would greatly reduce the expected mortality and number of infected individuals at 230 days post introduction of the exposed individual in the population. Figure 10 shows the dynamics after 200 days of mitigation (230 days post introduction of the disease in the community). According to this analysis, the number of people that would have been exposed/infected is 59 (0.18%), the number of people exposed would be 2.3 (0.007%), 5 infected individuals (0.016%) and 52 recovered individuals (0.16%). It is possible that no deaths would have occurred ($Dead_{120}=0.25$). This would represent a significant reduction of the number of cases to 17.6 % of what we predict without an intervention (21/119*100) at day 120 post introduction of the virus in the community.

![E1_2thirdsBeta_200days](https://tva1.sinaimg.cn/large/00831rSTgy1gct15v1pl4j30y10u0wlh.jpg)Figure 10. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals during the intial phase of the outbreak caused by 1 exposed individual assuming an overall mortality rate of 2.3% over a two month period ($\mu=0.0004$ *per day*). In this case we see the progression of the dynamics after mitigation, implemented 30 days after introduction of the disease, reduces transmission to 2/3 of its original value.

Now, this estimates change considerably if we now examine the same mitigation policy in a case in which 10 exposed individuals introduce the disease in the community. At 200 days post mitigation implementation (230 days post introduction of disease in the community) we would expect that 566 individuals would have been infected (1.75%, Figure 11). At this point in time we would have 21 people exposed (0.07%), 48 infected individuals (0.15%), 3 deads (0.008%) and 505 recovered individuals (1.56%). This means that having 10 early introductions make the mitigation taken 30 days post-introduction 10 times less effective than teh case with a single introduction (566 people infected vs 59). This is equivalent to say that if more cases are expected, then earlier action to reduce transmission is necessary for a better outcome.

![E10_2thirdsBeta_200days](https://tva1.sinaimg.cn/large/00831rSTgy1gct1dnpscqj30y10u0jyr.jpg)Figure 11. Dynamics of Susceptible, Exposed, Infected, and Recovered individuals during the intial phase of the outbreak caused by 10 exposed individuals assuming an overall mortality rate of 2.3% over a two month period ($\mu=0.0004$ *per day*). In this case we see the progression of the dynamics after mitigation, implemented 30 days after introduction of the disease, reduces transmission to 2/3 of its original value.

With this simple model, the consequences of reducing the transmission rate to this magnitude is that we expect the epidemic to extend much beyond the previous estimates. We expect that the number of infected individuals won't come down until many months after introduction of the disease in the community if it is not warm sensitive.

## Structured SEIR model for the spread of SARS-CoV-2

Up to this point I have not started examining simulations of more elaborate models, I am going to introduce the age structured model I have been started to develop. The motivation for this is my interest in assessing what is going to be the relative impact of the non-risk population to the outcome of the epidemic for those at risk. 

This is just a reminder of information already presented. At this point the mortality rate by SARS-CoV-2 seems to be remarkably different between individuals of different ages$^7$.  In this study, following 72,314 cases in Wuhan 44,672 (62%) were classified as confirmed cases of COVID-19, the number of suspected cases was 16 186 (22%), the number of diagnosed cases 10,567 (15%) and the number of asymptomatic 889 (1%). I am going to assume the distribution of fatality rate from$^7$:

  1. 2.3% (1023 of 44 672 confirmed cases)
  2. 14.8% in patients aged ≥80 years (208 of 1408)
  3. 8.0% in patients aged 70-79 years (312 of 3918)
  4. 49.0% in critical cases (1023 of 2087)

From these, We are going to separate the population in two strata: low risk (l) and high risk (h) and have a separate mortality rate ($\mu_{h}$) for individuals > 70 year old to be 8% or 11.4% in a 2 months period to model the spread of disease in an age structured population for individuals, and consider the mortality rate ($\mu_{l}$) for the rest of the population to be 2.3% in a 2 months period. The mortality rate for >70 yo would then be $\mu_{h}=0.0013$ or $\mu_{h}=0.019$ and the mortality rate of rest of the population is $\mu_{l}=0.0004$. To simplify matters, we will assume that for a period of three months (120 days) the mortality in the population younger than 70 can be neglected.

We then modified the previous model with equations 6-9 to include mortality and two different categories of individuals (higher risk with sub-index *h* and low risk with sub-index *l*):

$$
\begin{align}\dot{S_{h}}=-\beta S_{h}(I_{h}+I_{l})/N \\\dot{S_{l}}=-\beta S_{l}(I_{h}+I_{l})/N \\\dot{E_{h}}=\beta S_{h}(I_{h}+I_{h})/N-\sigma E_{h} \\\dot{E_{l}}=\beta S_{l}(I_{h}+I_{l})/N-\sigma E_{l} \\\dot{I_{h}}=\sigma E_{h}-(\gamma+\mu_{h})I_{h} \\\dot{I_{l}}=\sigma E_{l}-(\gamma + \mu_{l})I_{l}\\\dot{R_{h}}=\gamma I_{h} \\\dot{R_{l}}=\gamma I_{l}\\ \dot{deaths_{h}}=\mu_{h}I_{h} \\ \dot{deaths_{l}}=\mu_{l}I_{l}\\\end{align}
$$
Following the census information from Figure 1, we establish that for Pullman, the number of susceptible individuals older than 70 years is $S_{o}=1,653$ and the number of susceptible individuals younger than 70 is $S_{b}=30,729$.

## <font color = blue>I have yet to set the simulations for this model and establish different outcomes I want to assess</font>



You can see the code for the simulations below.

```{lang}Simple model 

##############################
# These are the simultaions for COVID-19
#############################

library(deSolve)


# some notes
# for looking at the numbers I am assuming that at anyt point in time t (let us say 300 days) I can invoke the values at the desired time, I can estimate the percent of population that are in each compartment by simply:


out_m_1[time=300,]/N*100

# and the percent of people that were ever infected are:

(100 - (out_m_1[time=300,]/N*100)[2])

# and the total number of people ever infected is

(1 - (out_m_1[time=300,]/N)[2])*N 


#############################################
## Simple model 
############################################


# case to consider sqrt(2.5*0.0583*(0.0583+0.0004)) = 0.0925. This is a note for Omar 
#  original estimation was b = 0.147

parameters <- c(b = 0.134,
	s = 0.192,      #1/5.2
	g = 0.0807,     #average 1/12.39=0.0807  or ceiling 1/17.16=0.0583
	m = 0.0004)

state <- c(X = 32382,
	Y = 1,
	Z = 0,
	W = 0)

N=32382

seir_covid <-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		dX <- (-b*X*Z)/N
		dY <- (b*X*Z)/N - s*Y
		dZ <- s*Y - g*Z
		dW <- g*Z
		# return the rate of change
		list(c(dX, dY, dZ, dW))
	}) # end with(as.list ...
}


times <- seq(0, 120, by = 1)


out_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_1) <- c("time","Susceptible","Exposed","Infected","Recovered")
head(out_1)

par(oma = c(0, 0, 4, 0))
plot(out_1, xlab = "time (days)", ylab = "Number of individuals")
#plot(out[, "Susceptible"], out[, "Infected"], pch = ".")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)


out_1[time=120,]/N*100
(100 - (out_1[time=120,]/N*100)[2])
(1 - (out_1[time=120,]/N)[2])*N


times <- seq(0, 420, by = 1)


out_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_1) <- c("time","Susceptible","Exposed","Infected","Recovered")
head(out_1)

par(oma = c(0, 0, 4, 0))
plot(out_1, xlab = "time (days)", ylab = "Number of individuals")
#plot(out[, "Susceptible"], out[, "Infected"], pch = ".")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)

out_1[time=421,]/N*100
(100 - (out_1[time=421,]/N*100)[2])
(1 - (out_1[time=421,]/N)[2])*N

out_1.1 <- as.data.frame(out_1)
out_1.1[which(out_2.1[,"Infected"] == max(out_2.1[,"Infected"])),]

#######################
##. 10 exposed
#######################

parameters <- c(b = 0.134,
	s = 0.192,      #1/5.2
	g = 0.0807,     #average 1/12.39=0.0807  or ceiling 1/17.16=0.0583
	m = 0.0004)

state <- c(X = 32382,
	Y = 10,
	Z = 0,
	W = 0)

N=32382

seir_covid <-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		dX <- (-b*X*Z)/N
		dY <- (b*X*Z)/N - s*Y
		dZ <- s*Y - g*Z
		dW <- g*Z
		# return the rate of change
		list(c(dX, dY, dZ, dW))
	}) # end with(as.list ...
}


times <- seq(0, 120, by = 1)


out_2 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_2) <- c("time","Susceptible","Exposed","Infected","Recovered")
head(out_2)

par(oma = c(0, 0, 4, 0))
plot(out_2, xlab = "time (days)", ylab = "Number of individuals")
#plot(out[, "Susceptible"], out[, "Infected"], pch = ".")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)

out_2[time=121,]/N*100
(100 - (out_2[time=121,]/N*100)[2])
(1 - (out_2[time=121,]/N)[2])*N


times <- seq(0, 420, by = 1)


out_2 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_2) <- c("time","Susceptible","Exposed","Infected","Recovered")
head(out_2)

par(oma = c(0, 0, 4, 0))
plot(out_2, xlab = "time (days)", ylab = "Number of individuals")
#plot(out[, "Susceptible"], out[, "Infected"], pch = ".")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)

out_2[time=421,]/N*100
(100 - (out_2[time=421,]/N*100)[2])
(1 - (out_2[time=421,]/N)[2])*N


out_2.1 <- as.data.frame(out_2)
out_2.1[which(out_2.1[,"Infected"] == max(out_2.1[,"Infected"])),]


##########################
# overall mortality included (not separated by age class yet)
#########################



######
# E = 1
#####

parameters <- c(b = 0.134,
	s = 0.192,      #1/5.2
	g = 0.0807,     #average 1/12.39=0.0807  or ceiling 1/17.16=0.0583
	m = 0.0004)


state <- c(X = 32382,
	Y = 1,
	Z = 0,
	Z0 = 0,
	W = 0)

N=32382

seir_covid <-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		dX <- (-b*X*Z)/N
		dY <- (b*X*Z)/N - s*Y
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
#plot(out[, "Susceptible"], out[, "Infected"], pch = ".", ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)

out_1[time=121,]/N*100
(100 - (out_1[time=121,]/N*100)[2])
(1 - (out_1[time=121,]/N)[2])*N


times <- seq(0, 420, by = 1)


out_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
head(out_1)

par(oma = c(0, 0, 4, 0))
plot(out_1, xlab = "time (days)", ylab = "Number of individuals")
#plot(out_1[, "Susceptible"], out_1[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)

out_1[time=421,]/N*100
(100 - (out_1[time=421,]/N*100)[2])
(1 - (out_1[time=421,]/N)[2])*N


out_1.1 <- as.data.frame(out_1)
out_1.1[which(out_1.1[,"Infected"] == max(out_1.1[,"Infected"])),]


######
# E = 10
#####

parameters <- c(b = 0.134,
	s = 0.192,      #1/5.2
	g = 0.0807,     #average 1/12.39=0.0807  or ceiling 1/17.16=0.0583
	m = 0.0004)


state <- c(X = 32382,
	Y = 10,
	Z = 0,
	Z0 = 0,
	W = 0)

N=32382

seir_covid <-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		dX <- (-b*X*Z)/N
		dY <- (b*X*Z)/N - s*Y
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
#plot(out[, "Susceptible"], out[, "Infected"], pch = ".", ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)


out_1[time=121,]/N*100
(100 - (out_1[time=121,]/N*100)[2])
(1 - (out_1[time=121,]/N)[2])*N


times <- seq(0, 420, by = 1)


out_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
head(out_1)

par(oma = c(0, 0, 4, 0))
plot(out_1, xlab = "time (days)", ylab = "Number of individuals")
#plot(out[, "Susceptible"], out[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)


out_1[time=421,]/N*100
(100 - (out_1[time=421,]/N*100)[2])
(1 - (out_1[time=421,]/N)[2])*N


out_1.1 <- as.data.frame(out_1)
out_1.1[which(out_1.1[,"Infected"] == max(out_1.1[,"Infected"])),]



#####################################
# transmission reduction efforts 2/3*b after 30 days of introduction for case E=1
####################################



parameters <- c(b = (2/3)*0.134,
	s = 0.192,      #1/5.2
	g = 0.0807,     #average 1/12.39=0.0807  or ceiling 1/17.16=0.0583
	m = 0.0004)


state <- c(X = 32378.35,
	Y = 0.9103652,
	Z = 1.5296350,
	Z0 = 1.089002e-02,
	W = 2.19706235)

N=32382


seir_covid <-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		dX <- (-b*X*Z)/N
		dY <- (b*X*Z)/N - s*Y
		dZ <- s*Y - (g + m)*Z
		dZ0 <- m*Z
		dW <- g*Z
		# return the rate of change
		list(c(dX, dY, dZ, dZ0, dW))
	}) # end with(as.list ...
}


times <- seq(0, 200, by = 1)


out_m_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_m_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
head(out_m_1)

par(oma = c(0, 0, 4, 0))
plot(out_m_1, xlab = "time (days after mitigation)", ylab = "Number of individuals")
#plot(out_m_1[, "Susceptible"], out_m_1[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)


out_m_1[time=201,]/N*100
(100 - (out_m_1[time=201,]/N*100)[2])
(1 - (out_m_1[time=201,]/N)[2])*N


#times <- seq(0, 1900, by = 1)


#out_m_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
#colnames(out_m_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
#head(out_m_1)

#par(oma = c(0, 0, 4, 0))
#plot(out_m_1, xlab = "time (days after mitigation)", ylab = "Number of individuals")
#plot(out_m_1[, "Susceptible"], out_m_1[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
#mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)


#out_m_1[time=1900,]/N*100
#(100 - (out_m_1[time=1900,]/N*100)[2])
#(1 - (out_m_1[time=1900,]/N)[2])*N


#out_1.1 <- as.data.frame(out_1)
#out_1.1[which(out_1.1[,"Infected"] == max(out_1.1[,"Infected"])),]




#####################################
# transmission reduction efforts 2/3*b after 30 days of introduction for case E=10
####################################



parameters <- c(b = (2/3)*0.134,
	s = 0.192,      #1/5.2
	g = 0.0807,     #average 1/12.39=0.0807  or ceiling 1/17.16=0.0583
	m = 0.0004)


state <- c(X = 32345.55,
	Y = 9.092178,
	Z = 15.285794,
	Z0 = 0.1088733830,
	W = 21.96520502)

N=32382

seir_covid <-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		dX <- (-b*X*Z)/N
		dY <- (b*X*Z)/N - s*Y
		dZ <- s*Y - (g + m)*Z
		dZ0 <- m*Z
		dW <- g*Z
		# return the rate of change
		list(c(dX, dY, dZ, dZ0, dW))
	}) # end with(as.list ...
}


times <- seq(0, 200, by = 1)


out_m_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
colnames(out_m_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
head(out_m_1)

par(oma = c(0, 0, 4, 0))
plot(out_m_1, xlab = "time (days after mitigation)", ylab = "Number of individuals")
#plot(out_m_1[, "Susceptible"], out_m_1[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)


out_m_1[time=201,]/N*100
(100 - (out_m_1[time=201,]/N*100)[2])
(1 - (out_m_1[time=201,]/N)[2])*N

#times <- seq(0, 540, by = 1)


#out_m_1 <- ode(y = state, times = times, func = seir_covid, parms = parameters)
#colnames(out_m_1) <- c("time","Susceptible","Exposed","Infected","Dead","Recovered")
#head(out_m_1)

#par(oma = c(0, 0, 4, 0))
#plot(out_m_1, xlab = "time (days after mitigation)", ylab = "Number of individuals")
#plot(out_m_1[, "Susceptible"], out_m_1[, "Infected"], pch = ".",ylab="Infected",xlab="Susceptibles")
#mtext(outer = TRUE, side = 3, "SEIR COVID-19", cex = 1.5)
`
```

### References

$^1$Wang H, Wang Z, Dong Y, et al. Phase-adjusted estimation of the number of Coronavirus Disease 2019 cases in Wuhan, China. _Cell Discov_. 2020;6:10. Published 2020 Feb 24. doi:10.1038/s41421-020-0148-0
$^2$Li, Q. et al. Early transmission dynamics in Wuhan, China, of novel coronavirus–infected pneumonia. _N. Engl. J. Med._  [https://doi.org/10.1056/NEJMoa2001316](https://doi.org/10.1056/NEJMoa2001316) (2020).
$^3$Chen, N. et al. Epidemiological and clinical characteristics of 99 cases of 2019 novel coronavirus pneumonia in Wuhan, China: a descriptive study. _Lancet_  [https://doi.org/10.1016/S0140-6736(20)30211-7](https://doi.org/10.1016/S0140-6736(20)30211-7).
$^4$Yang, Y. et al. Epidemiological and clinical features of the 2019 novel coronavirus outbreak in China. _medRxiv._  [https://doi.org/10.1101/2020.02.10.20021675](https://doi.org/10.1101/2020.02.10.20021675) (2020).
$^5$Read, J. M., Bridgen, J. R. E., Cummings, D. A. T., Ho, A. & Jewell, C. P. Novel coronavirus 2019-nCoV: early estimation of epidemiological parameters and epidemic predictions. _medrxiv_. [https://www.medrxiv.org/content/10.1101/2020.01.23.20018549v1.full.pdf](https://www.medrxiv.org/content/10.1101/2020.01.23.20018549v1.full.pdf) (2020).
$^6$Imai, N. et al. _Report 3: Transmissibility of 2019-nCoV._  [https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/news–wuhan-coronavirus/](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/news%E2%80%93wuhan-coronavirus/) (2020).
$^7$Wu Z, McGoogan JM. Characteristics of and Important Lessons From the Coronavirus Disease 2019 (COVID-19) Outbreak in China: Summary of a Report of 72 314 Cases From the Chinese Center for Disease Control and Prevention [published online ahead of print, 2020 Feb 24]. _JAMA_. 2020;10.1001/jama.2020.2648. doi:10.1001/jama.2020.2648