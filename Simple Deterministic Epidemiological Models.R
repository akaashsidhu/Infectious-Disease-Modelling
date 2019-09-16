#Coding and solving simple deterministic epidemiological models in R

#-------Load package to do numerical integration----
require(deSolve)
require(ggplot2)

#-------Setup a closed SIR model---------------
sir.model.closed <- function (t, x, params) { #here we begin a function with three arguments
   S <- x[1] #create local variable S, the first element of x
   I <- x[2] #create local variable I
   R <- x[3] #create local variable R
   with( #we can simplify code using "with"
     as.list(params), #this argument to "with" lets us use the variable names
     { #the system of rate equations
       dS <- -beta*S*I
       dI <- beta*S*I-gamma*I
       dR <- gamma*I
       dx <- c(dS,dI,dR) #combine results into a single vector dx
       list(dx) #return result as a list
       }
     )
   }

#--------Assign times, parameters and initial conditions--------------
times <- seq(0,120,by=5) #function seq returns a sequence
params <- c(beta=0.3,gamma=1/7) #function c "c"ombines values into a vector
xstart <- c(S=9999/10000,I=1/10000,R=0) #initial conditions

#---------Simulate a model trajectory------------

out <- as.data.frame(lsoda(xstart,times,sir.model.closed,params)) #result stored in dataframe

#----------Plot the model results-------------

plot(I~time,data=out,type='b') 
#plot the I variable against time Notice: this is a simple plot not using ggplot

#Exercise 1. Explore the dynamics of the system for different values of the β and γ parameters by simulating and plotting trajectories as time ----

times <- seq(0,120,by=5) #Keep the same time in days where there is an increase by 5 starting at 0 and ending at 120.
params <- c(beta=0.3,gamma=1/7) #The original beta and gamma value.
params2 <- c(beta=0.9,gamma=1/7) #Increase the beta and keep the gamma value the same.
params3 <- c(beta=0.3,gamma=1/5) #Increase the gamma and keep the beta value the same.
params4 <- c(beta=0.9,gamma=1/5) #Increase both gamma and beta.
xstart <- c(S=9999/10000,I=1/10000,R=0) #Keep the initial conditions the same.

out <- as.data.frame(lsoda(xstart,times,sir.model.closed,params)) 
out2 <- as.data.frame(lsoda(xstart,times,sir.model.closed,params2))
out3 <- as.data.frame(lsoda(xstart,times,sir.model.closed,params3))
out4 <- as.data.frame(lsoda(xstart,times,sir.model.closed,params4))
#Store the result of each in a dataframe to be used in the plots. 

Ex1 <-ggplot() + 
  geom_line(size = 1.2, data=out, aes(x=time, y=I, colour = "Original")) +
  geom_line(size = 1.2, data=out2, aes(x=time, y=I, colour = "Beta Increased" )) +
  geom_line(size = 1.2, data=out3, aes(x=time, y=I, colour = "Gamma Increased" )) +
  geom_line(size = 1.2, data=out4, aes(x=time, y=I, colour = "Both Increased" )) +
  ylim(0, 0.60) +
  xlim (0,120) +
  labs(title="Exploration of β and γ Parameters in the SIR Model", x = "Time (days)",
       y = "I") +
  theme_bw() +
  theme(text=element_text(size = 12)) +
  theme(legend.justification=c(0.99,0.99), 
        legend.position=c(0.99,0.99), 
        legend.title = element_text(colour="BLACK", size=12),
        text=element_text(size = 12)) +
  theme(legend.title=element_blank())
Ex1
#Created a graph with 4 seperate lines to show the effect of an increase in Beta/Gamma/Both.

#Exercise 2. Explore the dynamics of the system for one set of β and γ at different initial conditions. What happens if there is pre-existing immunity in the population?----

times <- seq(0,120,by=5) #Keep the same time in days where there is an increase by 5 starting at 0 and ending at 120.
params <- c(beta=0.3,gamma=1/7) #The original beta and gamma value.
xchange <- c(S=8999/10000,I=1/10000,R=1000/10000) #Change the initial conditions by increasing R to represent pre-existing immunity in the population.

outchange <- as.data.frame(lsoda(xchange,times,sir.model.closed,params)) 
#Store into a dataframe to be used during plotting.

Ex2 <-ggplot() + 
  geom_line(size = 1.2, data=out, aes(x=time, y=I, colour = "Original")) +
  geom_line(size = 1.2, data=outchange, aes(x=time, y=I, colour = "R Increased" )) +
  ylim(0, 0.60) +
  xlim (0,120) +
  labs(title="The Effects of Pre-existing Immunity in a Population (n=10000)", x = "Time (days)", y = "I") +
  theme_bw() +
  theme(text=element_text(size = 12)) +
  theme(legend.justification=c(0.99,0.99), legend.position=c(0.99,0.99), legend.title = element_text(colour="BLACK", size=12),text=element_text(size = 12)) +
  theme(legend.title=element_blank())
Ex2
#Created a graph with 2 lines: Original initial values and Increased R Values to show the effects of pre-existing immunity in a population. 

#Exercise 3. Modify the code given to study the dynamics of a demographically open SIR model----
openmodel <- function (t, x, parameters) { 
  S <- x[1] #Created a compartment for Susceptible
  I <- x[2] #Created a compartment for Infected
  R <- x[3] #Created a compartment for Recovered
  with( 
    as.list(parameters), 
    {
      dS <-  - (inf.prob*act.rate)*S*I/(S+I+R)  - ds.rate*S + b.rate*(S+I+R)
      dI <- (inf.prob*act.rate)*S*I/(S+I+R) - rec.rate*I - di.rate*I
      dR <- rec.rate*I - dr.rate*R
      dx <- c(dS,dI,dR) 
      list(dx) 
    }
  )
}
#Use an open model as demonstrated in Lab 2. There is the addition of a birth rate and death rate (one for each compartment) to the original SIR model which was closed.

time1 <- seq(0,365, by=1) #Set time1 to 365 days
time2 <- seq(0, 3650, by=1) #Set time2 to 3650 days
time3 <- seq(0, 36500, by=1) #Set time3 to 36500 days
parameters <- c(inf.prob = 0.2,act.rate = 1, rec.rate = 1/20, 
            b.rate = 1/95, ds.rate = 1/100, di.rate = 1/80, 
            dr.rate=1/100)
initial <- c(S=1000, I=1, R=0) #initial conditions

output1 <- as.data.frame(lsoda(initial,time1,openmodel,parameters)) 
output2 <- as.data.frame(lsoda(initial,time2,openmodel,parameters))
output3 <- as.data.frame(lsoda(initial,time3,openmodel,parameters))
#Save the results into a dataframe for each seperate time.

model1 <-ggplot() + 
  geom_line(size = 1.2, data=output1, aes(x=time, y=I, colour = "Infected")) +
  geom_line(size = 1.2, data=output1, aes(x=time, y=S, colour = "Susceptible" )) +
  geom_line(size = 1.2, data=output1, aes(x=time, y=R, colour = "Recovered" )) +
  ylim(0, 1002) +
  xlim (0, 365) +
  labs(title="Open SIR Model: 1 Year", x = "Time (days)",
       y = "Number of individuals") +
  theme_bw() +
  theme(text=element_text(size = 12)) +
  theme(legend.justification=c(0.99,0.99), 
        legend.position=c(0.99,0.99), 
        legend.title = element_text(colour="BLACK", size=12),
        text=element_text(size = 12)) +
  theme(legend.title=element_blank())

model1
#Created a graph with t = 365 days with 3 lines for SIR.

model2 <-ggplot() + 
  geom_line(size = 1.2, data=output2, aes(x=time, y=I, colour = "Infected")) +
  geom_line(size = 1.2, data=output2, aes(x=time, y=S, colour = "Susceptible" )) +
  geom_line(size = 1.2, data=output2, aes(x=time, y=R, colour = "Recovered" )) +
  ylim(0, 1002) +
  xlim (0, 3650) +
  labs(title="Open SIR Model: 10 Years", x = "Time (days)",
       y = "Number of individuals") +
  theme_bw() +
  theme(text=element_text(size = 12)) +
  theme(legend.justification=c(0.99,0.99), 
        legend.position=c(0.99,0.99), 
        legend.title = element_text(colour="BLACK", size=12),
        text=element_text(size = 12)) +
  theme(legend.title=element_blank())

model2
#Created a graph with t = 3650 days with 3 lines for SIR.

model3 <-ggplot() + 
  geom_line(size = 1.2, data=output3, aes(x=time, y=I, colour = "Infected")) +
  geom_line(size = 1.2, data=output3, aes(x=time, y=S, colour = "Susceptible" )) +
  geom_line(size = 1.2, data=output3, aes(x=time, y=R, colour = "Recovered" )) +
  ylim(0, 1002) +
  xlim (0, 36500) +
  labs(title="Open SIR Model: 100 Years", x = "Time (days)",
       y = "Number of individuals") +
  theme_bw() +
  theme(text=element_text(size = 12)) +
  theme(legend.justification=c(0.99,0.99), 
        legend.position=c(0.99,0.99), 
        legend.title = element_text(colour="BLACK", size=12),
        text=element_text(size = 12)) +
  theme(legend.title=element_blank())

model3
#Created a graph with t = 36500 days with 3 lines for SIR.

#Exercise 4. Modify the codes given to study the dynamics of an SEIR model.----
#For this question, use the following parameter values for the latent period: 0.0001 days, 3.5 days, 7 days, and 14 days.

SEIRmodel <- function (time, y, parameters2)
{
  S = y[1]#Created a compartment for Susceptible
  E = y[2]#Created a compartment for Exposed
  I = y[3]#Created a compartment for Infected
  R = y[4]#Created a compartment for Recovered
  with ( 
    as.list (parameters2),
    {
      dS <- - (inf.prob*act.rate)*S*I/(S+I+E+R)  - ds.rate*S + b.rate*(S+I+E+R)
      dE <- (inf.prob*act.rate)*S*I/(S+I+E+R) - (1/latent_period)*E - de.rate*E
      dI <- (1/latent_period)*E - rec.rate*I - di.rate*I
      dR <- rec.rate*I - dr.rate*R
      dres <- c(dS, dE, dI, dR)
      list (dres)
    }
  )
}
#Set up the open SEIR Model with an E compartment, birth rate, death rate, and latent period

par1 <- c(inf.prob = 0.2,act.rate = 1, rec.rate = 1/20, 
                b.rate = 1/95, ds.rate = 1/100, de.rate = 1/90,di.rate = 1/80, 
                dr.rate = 1/100, latent_period = 0.0001) 
par2 <- c(inf.prob = 0.2,act.rate = 1, rec.rate = 1/20, 
          b.rate = 1/95, ds.rate = 1/100, de.rate = 1/90,di.rate = 1/80, 
          dr.rate = 1/100, latent_period = 3.5) 
par3 <- c(inf.prob = 0.2,act.rate = 1, rec.rate = 1/20, 
          b.rate = 1/95, ds.rate = 1/100, de.rate = 1/90,di.rate = 1/80, 
          dr.rate = 1/100, latent_period = 7) 
par4 <- c(inf.prob = 0.2,act.rate = 1, rec.rate = 1/20, 
          b.rate = 1/95, ds.rate = 1/100, de.rate = 1/90,di.rate = 1/80, 
          dr.rate = 1/100, latent_period = 14) 
#Set up the parameters with the 4 different latent periods. Added a death rate for the E compartment.

initial_value <- c(S=1000, E=0, I=1, R=0) #Initial conditions
time <- seq (0, 365, by=1) #Set up time for 365 days.
#Set up the initial values and time. 

ex4_1 <- as.data.frame(lsoda(initial_value, time, SEIRmodel, par1))
ex4_2 <- as.data.frame(lsoda(initial_value, time, SEIRmodel, par2))
ex4_3 <- as.data.frame(lsoda(initial_value, time, SEIRmodel, par3))
ex4_4 <- as.data.frame(lsoda(initial_value, time, SEIRmodel, par4))
#Set up the dataframes to be plotted on the graphs to show the effects of the latent period.

plotEx4_1 <-ggplot() + 
  geom_line(size = 1.2, data=ex4_1, aes(x=time, y=S, colour = "Susceptible")) +
  geom_line(size = 1.2, data=ex4_1, aes(x=time, y=E, colour = "Exposed" )) +
  geom_line(size = 1.2, data=ex4_1, aes(x=time, y=I, colour = "Infected" )) +
  geom_line(size = 1.2, data=ex4_1, aes(x=time, y=R, colour = "Recovered" )) +
  ylim(0, 1200) +
  xlim (0,365) +
  labs(title="Exploration of the Effects of Latent Period: 0.0001 Days", x = "Time (days)",
       y = "Population") +
  theme_bw() +
  theme(text=element_text(size = 12)) +
  theme(legend.justification=c(0.99,0.99), 
        legend.position=c(0.99,0.99), 
        legend.title = element_text(colour="BLACK", size=12),
        text=element_text(size = 12)) +
  theme(legend.title=element_blank())
#Plot the first latent period with the different parameters (S E I R)

plotEx4_2 <-ggplot() + 
  geom_line(size = 1.2, data=ex4_2, aes(x=time, y=S, colour = "Susceptible")) +
  geom_line(size = 1.2, data=ex4_2, aes(x=time, y=E, colour = "Exposed" )) +
  geom_line(size = 1.2, data=ex4_2, aes(x=time, y=I, colour = "Infected" )) +
  geom_line(size = 1.2, data=ex4_2, aes(x=time, y=R, colour = "Recovered" )) +
  ylim(0, 1200) +
  xlim (0,365) +
  labs(title="Exploration of the Effects of Latent Period: 3.5 Days", x = "Time (days)",
       y = "Population") +
  theme_bw() +
  theme(text=element_text(size = 12)) +
  theme(legend.justification=c(0.99,0.99), 
        legend.position=c(0.99,0.99), 
        legend.title = element_text(colour="BLACK", size=12),
        text=element_text(size = 12)) +
  theme(legend.title=element_blank())
#Plot the second latent period with the different parameters (S E I R)

plotEx4_3 <-ggplot() + 
  geom_line(size = 1.2, data=ex4_3, aes(x=time, y=S, colour = "Susceptible")) +
  geom_line(size = 1.2, data=ex4_3, aes(x=time, y=E, colour = "Exposed" )) +
  geom_line(size = 1.2, data=ex4_3, aes(x=time, y=I, colour = "Infected" )) +
  geom_line(size = 1.2, data=ex4_3, aes(x=time, y=R, colour = "Recovered" )) +
  ylim(0, 1200) +
  xlim (0,365) +
  labs(title="Exploration of the Effects of Latent Period: 7 Days", x = "Time (days)",
       y = "Population") +
  theme_bw() +
  theme(text=element_text(size = 12)) +
  theme(legend.justification=c(0.99,0.99), 
        legend.position=c(0.99,0.99), 
        legend.title = element_text(colour="BLACK", size=12),
        text=element_text(size = 12)) +
  theme(legend.title=element_blank())
#Plot the third latent period with the different parameters (S E I R)

plotEx4_4 <-ggplot() + 
  geom_line(size = 1.2, data=ex4_4, aes(x=time, y=S, colour = "Susceptible")) +
  geom_line(size = 1.2, data=ex4_4, aes(x=time, y=E, colour = "Exposed" )) +
  geom_line(size = 1.2, data=ex4_4, aes(x=time, y=I, colour = "Infected" )) +
  geom_line(size = 1.2, data=ex4_4, aes(x=time, y=R, colour = "Recovered" )) +
  ylim(0, 1200) +
  xlim (0,365) +
  labs(title="Exploration of the Effects of Latent Period: 14 Days", x = "Time (days)",
       y = "Population") +
  theme_bw() +
  theme(text=element_text(size = 12)) +
  theme(legend.justification=c(0.99,0.99), 
        legend.position=c(0.99,0.99), 
        legend.title = element_text(colour="BLACK", size=12),
        text=element_text(size = 12)) +
  theme(legend.title=element_blank())
#Plot the fourth latent period with the different parameters (S E I R)

plotEx4_1
plotEx4_2
plotEx4_3
plotEx4_4
#All 4 plots for Exercise 4 to show the effects of different latent periods on the different parameters.
