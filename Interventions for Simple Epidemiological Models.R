#Adding interventions into simple deterministic epidemiological models in R - Vaccination

#----Load Packages into R----
require(deSolve)
library("ggplot2")

#----Exercise 3: Adapt the R code to represent the model that you have described above. Let’s assume the scenario where you have the introduction of a new strain of virus into the human population. Assume that this pathogen has an R0 of 2.0. In reality, it can take quite a long time for new vaccine to be produced and made available for public distribution. Assume that it takes 25 days before vaccine is available. Assume you have two different scenarios. In one case, you could only vaccinate 10% of the population with a vaccine that is not very effective (VE = 10%) and in the second case, you could vaccinate 70% of the population with a vaccine that is 40% effective. Generate a figure to show the difference in the epidemic trajectory for each scenario and discuss your findings----

#----Define parameters, initial conditions, time, vaccinated proportion, and vaccine effectiveness----

#The following parameters are defined:
R0 <- 2.0 #The basic reproductive number of a new pathogen
N <- 1 #The population size
gamma <- 1/5 #The recovery rate in days
beta <- R0*(gamma) #The transmission rate (From S -> I)

#The proportion of susceptibles that are vaccinated is defined:
vax.init <- 0.0 #The initial proportion of susceptibles vaccinated is 0 because vaccine isn't introduced until day = 25.
vax.s1 <- 0.10 #In this scenario, a small proportion of the population is vaccinated.
vax.s2 <- 0.70 #In this scenario, a large proportion of the population is vaccinated. 

#The effectiveness of the vaccine is defined:
eff.init <- 0.0 #The initial effectiveness of the vaccine is 0 because it isn't introduced until day = 25
eff.s1 <- 0.10 #In this scenario, the vaccine is not very effective.
eff.s2 <- 0.40 #In this scenario, the vaccine is more effective than in scenario 1. 

#The initial conditions and time of the SIRV model is defined:
times <- seq(0,120,by=5) #Returns a sequence
xstart <- c(S=9999/10000,I=1/10000,R=0,V=0) #The initial conditions 
params.init <- c(beta=beta,gamma=gamma,vax=vax.init,eff=eff.init) #Parameters to be used for the initial conditions 
params.s1 <- c(beta=beta,gamma=gamma,vax=vax.s1, eff=eff.s1) #Parameters to be used for scenario one
params.s2 <- c(beta=beta,gamma=gamma,vax=vax.s2, eff=eff.s2) #Parameters to be used for scenario two. 
params <- c(beta=beta,gamma=gamma)

#----Set up the Closed SIRV Model----
#This is the closed SIRV model that will be used to simulate the model trajectory. It includes the addition of a V compartment, for the vaccinated individuals, and two new parameters: vaccination rate (vax) and vaccine effectiveness (eff). 
sirv.model <- function (t, x, params) {
  S <- x[1] #create local variable S
  I <- x[2] #create local variable I
  R <- x[3] #create local variable R
  V <- x[4] #create local variable V
  
  with( 
    as.list(params),
    { 
      dS <- -beta*S*I - (vax*eff*S)
      dI <- beta*S*I-gamma*I
      dR <- gamma*I
      dV <- vax*eff*S
      dx <- c(dS,dI,dR,dV)
      list(dx) 
    }
  )
}

#----Simulate Disease Trajectory: Initial Conditions 0-25 days No Vaccine----

out1 <- ode(xstart,seq(0,25,by=1),sirv.model,params.init,method='ode45',rtol=1e-7)
#Demonstrates the initial conditions before the vaccine has been implemented. This represents day 0 to 25 without vaccine. 

#----Simulate Disease Trajectory: Vaccine Intervention At Day = 25----

out2 <- ode(tail(out1,1)[2:5],seq(25,365,by=1),sirv.model,params.s1,method='ode45',rtol=1e-7)
out2 <-out2[-1,]
#The vaccine has been implemented at day 25 and the disease trajectory of scenario one is simulated.

out2.s2 <- ode(tail(out1,1)[2:5],seq(25,365,by=1),sirv.model,params.s2,method='ode45',rtol=1e-7)
out2.s2 <- out2.s2[-1,]
#The vaccine has been implemented at day 25 and the disease trajectory of scenario two is simulated.

scenario.1 <- as.data.frame(rbind(out1, out2))
scenario.2 <- as.data.frame(rbind(out1, out2.s2))
#Create a dataframe of the two scenario to be used later for plotting purposes. 

#----Plotting Disease Trajectory For Scenario 1 and 2----

plot(I~time,data=scenario.1,type='b') 
plot(I~time,data=scenario.2,type='b') 
#The plots show the effects these scenarios have on the Infected over time. 

sirv.1 <- ggplot() + 
  geom_line(size = 1.2, data=scenario.1, aes(x=time, y=S, colour = "S")) +
  geom_line(size = 1.2, data=scenario.1, aes(x=time, y=I, colour = "I")) +
  geom_line(size = 1.2, data=scenario.1, aes(x=time, y=R, colour = "R")) +
  geom_line(size = 1.2, data=scenario.1, aes(x=time, y=V, colour = "V")) +
  labs(y="Population") + labs(x="Time") + ggtitle("Scenario 1: SIRV Rates over 1 year (365 days)")
sirv.1 
#The disease trajectory for scenario one is plotted with SIRV curves. 

sirv.2 <- ggplot() + 
  geom_line(size = 1.2, data=scenario.2, aes(x=time, y=S, colour = "S")) +
  geom_line(size = 1.2, data=scenario.2, aes(x=time, y=I, colour = "I")) +
  geom_line(size = 1.2, data=scenario.2, aes(x=time, y=R, colour = "R")) +
  geom_line(size = 1.2, data=scenario.2, aes(x=time, y=V, colour = "V")) +
  labs(y="Population") + labs(x="Time") + ggtitle("Scenario 2: SIRV Rates over 1 year (365 days)")
sirv.2 
#The disease trajectory for scenario two is plotted with SIRV curves. 

#----Exercise 4. Does your answer to the question above change if the vaccine is available later (e.g. 50 days after the introduction of the pathogen instead of 25days)? If yes, how and why?----

out3 <- ode(xstart,seq(0,50,by=1),sirv.model,params.init,method='ode45',rtol=1e-7)
#Demonstrates the initial conditions before the vaccine has been implemented. This represents day 0 to 50 without vaccine. 

out4 <- ode(tail(out3,1)[2:5],seq(50,365,by=1),sirv.model,params.s1,method='ode45',rtol=1e-7)
out4 <-out4[-1,]
#The vaccine has been implemented at day 50 and the disease trajectory of scenario one is simulated.

out4.s2 <- ode(tail(out3,1)[2:5],seq(50,365,by=1),sirv.model,params.s2,method='ode45',rtol=1e-7)
out4.s2 <- out4.s2[-1,]
#The vaccine has been implemented at day 50 and the disease trajectory of scenario two is simulated.

scenario.1.ex4 <- as.data.frame(rbind(out3, out4))
scenario.2.ex4 <- as.data.frame(rbind(out3, out4.s2))
#Create a dataframe of the two scenario to be used later for plotting purposes. 

plot(I~time,data=scenario.1.ex4,type='b') 
plot(I~time,data=scenario.2.ex4,type='b') 
#The plots show the effects these scenarios have on the Infected over time. 

sirv.3 <- ggplot() + 
  geom_line(size = 1.2, data=scenario.1.ex4, aes(x=time, y=S, colour = "S")) +
  geom_line(size = 1.2, data=scenario.1.ex4, aes(x=time, y=I, colour = "I")) +
  geom_line(size = 1.2, data=scenario.1.ex4, aes(x=time, y=R, colour = "R")) +
  geom_line(size = 1.2, data=scenario.1.ex4, aes(x=time, y=V, colour = "V")) +
  labs(y="Population") + labs(x="Time") + ggtitle("Scenario 1: SIRV Vaccination at Day = 50")
sirv.3 
#The disease trajectory for scenario one, where vaccine is introduced at day = 50, is plotted with SIRV curves. 

sirv.4 <- ggplot() + 
  geom_line(size = 1.2, data=scenario.2.ex4, aes(x=time, y=S, colour = "S")) +
  geom_line(size = 1.2, data=scenario.2.ex4, aes(x=time, y=I, colour = "I")) +
  geom_line(size = 1.2, data=scenario.2.ex4, aes(x=time, y=R, colour = "R")) +
  geom_line(size = 1.2, data=scenario.2.ex4, aes(x=time, y=V, colour = "V")) +
  labs(y="Population") + labs(x="Time") +ggtitle("Scenario 2: SIRV Vaccination at Day = 50")
sirv.4
#The disease trajectory for scenario two, where vaccine is introduced at day = 50, is plotted with SIRV curves. 

#----Exercise 5:Does your answer change if the pathogen is much more transmissible (e.g. R0 = greater than 2.0)? If yes, how and why?----

R0_2 <- 3.0 #New basic reproductive number that is greater than 2.0
beta_2 <- R0_2*(gamma) #New beta value using the new R0

params.init.2 <- c(beta=beta_2,gamma=gamma,vax=vax.init,eff=eff.init)
params.s1.2 <- c(beta=beta_2,gamma=gamma,vax=vax.s1, eff=eff.s1)
params.s2.2 <- c(beta=beta_2,gamma=gamma,vax=vax.s2, eff=eff.s2)
params.2 <- c(beta=beta_2,gamma=gamma)
#Set up the new parameters with the new value of beta.

out5 <- ode(xstart,seq(0,25,by=1),sirv.model,params.init.2,method='ode45',rtol=1e-7)
#Demonstrates the initial conditions before the vaccine has been implemented. This represents day 0 to 25 without vaccine and an increased R0. 

out6 <- ode(tail(out5,1)[2:5],seq(25,365,by=1),sirv.model,params.s1.2,method='ode45',rtol=1e-7)
out6 <-out6[-1,]
#The vaccine has been implemented at day 25 and the disease trajectory of scenario one is simulated, where R0 = 3.00

out6.s2 <- ode(tail(out5,1)[2:5],seq(25,365,by=1),sirv.model,params.s2.2,method='ode45',rtol=1e-7)
out6.s2 <- out6.s2[-1,]
#The vaccine has been implemented at day 25 and the disease trajectory of scenario two is simulated, where R0 = 3.00

scenario.5 <- as.data.frame(rbind(out5, out6))
scenario.6 <- as.data.frame(rbind(out5, out6.s2))
#Save the scenarios into a dataframe to be used later for plotting purposes. 

plot(I~time,data=scenario.5,type='b') 
plot(I~time,data=scenario.6,type='b') 
#New plots to check I over Time if the R0 was greater than 2. 

sirv.5 <- ggplot() + 
  geom_line(size = 1.2, data=scenario.5, aes(x=time, y=S, colour = "S")) +
  geom_line(size = 1.2, data=scenario.5, aes(x=time, y=I, colour = "I")) +
  geom_line(size = 1.2, data=scenario.5, aes(x=time, y=R, colour = "R")) +
  geom_line(size = 1.2, data=scenario.5, aes(x=time, y=V, colour = "V")) +
  labs(y="Population") + labs(x="Time") + ggtitle("Scenario 1: SIRV rates over 1 year (365 days) where R0 = 3")
sirv.5 
#The disease trajectory for scenario one is plotted with SIRV curves, where R0 = 3 

sirv.6 <- ggplot() + 
  geom_line(size = 1.2, data=scenario.6, aes(x=time, y=S, colour = "S")) +
  geom_line(size = 1.2, data=scenario.6, aes(x=time, y=I, colour = "I")) +
  geom_line(size = 1.2, data=scenario.6, aes(x=time, y=R, colour = "R")) +
  geom_line(size = 1.2, data=scenario.6, aes(x=time, y=V, colour = "V")) +
  labs(y="Population") + labs(x="Time") + ggtitle("Scenario 2: SIRV rates over 1 year (365 days) where R0 = 3")
sirv.6
#The disease trajectory for scenario two is plotted with SIRV curves, where R0 = 3 
