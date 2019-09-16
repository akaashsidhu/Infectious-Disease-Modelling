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
