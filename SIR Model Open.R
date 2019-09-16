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
