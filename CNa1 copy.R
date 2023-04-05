

library(deSolve)
install.packages("pracma")
library(pracma)
require(plot3D)
library(phaseR)


HH.time <- function(t,state, parameters){
  with(as.list(c(state, parameters)), {
    
    alpha.n <- function(V) (0.01*(V+55))/(1-exp(-0.1*(V+55)))
    beta.n <- function(V) 0.125*exp(-0.0125*(V+65))
    
    alpha.m <- function(V) (0.1*(V+40)) / (1-exp(-0.1*(V+40)))
    beta.m <- function(V) 4*exp(-0.0556*(V+65))
    
    alpha.h <- function(V) 0.07*exp(-0.05*(V+65))
    beta.h <- function(V) 1/(1+exp(-0.1*(V+35)))
   
    tau.n <- function(V) 1/(alpha.n(V) + beta.n(V))
    tau.m <- function(V) 1/(alpha.m(V) + beta.m(V))
    tau.h <- function(V) 1/(alpha.h(V) + beta.h(V))
    
    n.inf <- function(V) alpha.n(V)/(alpha.n(V) + beta.n(V))
    m.inf <- function(V) alpha.m(V)/(alpha.m(V) + beta.m(V))
    h.inf <- function(V) alpha.h(V)/(alpha.h(V) + beta.h(V))
    
    
    im <-  (g.na*(M^3)*H*(V-e.na)) + (g.k*(N^4)*(V-e.k)) + (g.l*(V-e.l))
    
    ie.a <- ifelse(t>=5, 0, -50 )
    
    dV <- (ie.a - im )/cm
    dN <- (n.inf(V) - N)/tau.n(V)
    dM <- (m.inf(V) - M)/tau.m(V)
    dH <- (h.inf(V) - H)/tau.h(V)
  
   
    
    return(list(c(dV, dN, dM, dH)))
  })
}

parameters <- c(cm =10, e.na = 50, e.k = -77, e.l = -54, 
                g.na = 1200, g.k = 360, g.l =3, ie.a = 0 )
state <- c(V = -70, N=0.3177, M = 0.0529, H = 0.5961 )
times <- seq(0, 40, by = 0.1)
out <- ode(y = state, times = times, func = HH.time, parms = parameters)
plot(out, xlab = "time /ms", ylab=c("mV", "", "", ""))



find_peaks <- function (x, m = 3){
  # x: vector of sequence
  # m: number of points on either side of peak that must be smaller than peak
  # return: index of peaks in sequence
  x <- ifelse(x>-54.4, x, -80)
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}




ISI <- function(spike.train, res){
  # spike.train: vector of spike times
  # res: resolution of times (i.e. how many points per millisecond)
  # return: 
  #     - inter spike intervals (ms) for each spike
  #     - average inter spike interval
  #     - average firing rate 
  isi <- integer(length(peaks)-1)
  if (length(spike.train) >2){
    isi <- diff(spike.train)/res
    ave.ISI <- mean(isi)
    ave.FR <- 1000/ave.ISI #1000ms in a second
  }else {
    isi <- NA
    ave.ISI <- NA
    ave.FR <- 1
  }
  
  return(list(isi, ave.ISI, ave.FR ))
}




fr <- c()
for (i in seq(0, 500,10)){
HH.time <- function(t,state, parameters){
  # t: sequence from 0 to T with a certain resolution
  # state: initial values for V, M, H, and N
  # parameters: set values for differential equations i.e. cm, e.na etc. 
  # return: values of V, M, H, and N over time t as matrix
  with(as.list(c(state, parameters)), {
    
    
    alpha.n <- function(V) (0.01*(V+55))/(1-exp(-0.1*(V+55)))
    beta.n <- function(V) 0.125*exp(-0.0125*(V+65))
    
    alpha.m <- function(V) (0.1*(V+40)) / (1-exp(-0.1*(V+40)))
    beta.m <- function(V) 4*exp(-0.0556*(V+65))
    
    alpha.h <- function(V) 0.07*exp(-0.05*(V+65))
    beta.h <- function(V) 1/(1+exp(-0.1*(V+35)))
    
    tau.n <- function(V) 1/(alpha.n(V) + beta.n(V))
    tau.m <- function(V) 1/(alpha.m(V) + beta.m(V))
    tau.h <- function(V) 1/(alpha.h(V) + beta.h(V))
    
    n.inf <- function(V) alpha.n(V)/(alpha.n(V) + beta.n(V))
    m.inf <- function(V) alpha.m(V)/(alpha.m(V) + beta.m(V))
    h.inf <- function(V) alpha.h(V)/(alpha.h(V) + beta.h(V))
    
    
    im <-  (g.na*(M^3)*H*(V-e.na)) + (g.k*(N^4)*(V-e.k)) + (g.l*(V-e.l))
    
    dV <- (ie.a - im )/cm
    dN <- (n.inf(V) - N)/tau.n(V)
    dM <- (m.inf(V) - M)/tau.m(V)
    dH <- (h.inf(V) - H)/tau.h(V)
    
    return(list(c(dV, dN, dM, dH)))
  })
}

parameters <- c(cm =10, e.na = 50, e.k = -77, e.l = -54, 
                g.na = 1200, g.k = 360, g.l =3, ie.a = i )
state <- c(V = -65, N=0.3177, M = 0.0529, H = 0.5961 )
times <- seq(0, 1000, by = 0.1)
out <- ode(y = state, times = times, func = HH.time, parms = parameters)

spike.train <- find_peaks(out[,2])
isi.results <- ISI(spike.train, 10) # to be set as 1/('by" value of "times" variable in seq function)
ave.fr <- isi.results[[3]]
fr <- c(fr, ave.fr)
par(mfrow=c(1,1))
}


par(mfrow=c(1,1))
par(mar=c(4,4,4,4))
plot(fr,pch=19, xlab=expression(paste("Ie/A (nA/mm"^"2", ")", sep="")), 
     ylab="average firing rate (Hz)", xaxt = "n")
axis(side=1, at=c(0, 10, 20, 30, 40, 50), labels=c("0", "100", "200", "300", "400", "500") )
abline(v=7.5, col="lightblue", lty="longdash")







## coupled integrate fire


integrate.fire <- function(time, h, parameters, v.one.init=-68, v.two.init=-70 ){
  with(as.list(parameters),{ 
  
  V1 <- integer(length(seq(0, time, h)))
  V2 <- integer(length(seq(0, time, h)))
  Ps1 <- integer(length(seq(0, time, h)))
  Ps2 <- integer(length(seq(0, time, h)))
  Z1 <- integer(length(seq(0, time, h)))
  Z2 <- integer(length(seq(0, time, h)))
  
  set.seed(4123)
  V1[1] <- v.one.init #runif(min=v.reset,max=v.th, 1)
  V2[1] <- v.two.init #runif(min=v.reset,max=v.th, 1)
  Ps1[1] <- 0.1
  Ps2[1] <- 0.1
  Z1[1] <- 1
  Z2[1] <- 1
  
  fr1 <- c()
  fr2 <-c()
  
  for (t in 1:(length(seq(0, time, h))-1)){
    dPs <- function(Ps, Z) ((exp(1)*pmax*Z - Ps) /tau.s)
    dZ <- function(Z) -Z/tau.s
    dV <- function(V, Ps) (e.l - V - rm.gs*Ps*(V - e.s) + bigrm.ie)/tau.m
    
    V1[t+1] <- V1[t] + h*dV(V1[t], Ps1[t])
    V2[t+1] <- V2[t] + h*dV(V2[t], Ps2[t])
    
    Ps1[t+1] <- Ps1[t] + h*dPs(Ps1[t], Z1[t])
    Ps2[t+1] <- Ps2[t] + h*dPs(Ps2[t], Z2[t])
    
    Z1[t+1] <- Z1[t] + h*dZ(Z1[t])
    Z2[t+1] <- Z2[t] + h*dZ(Z2[t])

    if (V1[t+1] > v.th){
      Z2[t+1] <- 1
      V1[t+1] <- v.reset
      fr2 <- c(fr2, t+1)
      
    }
    
    if (V2[t+1] > v.th){
      Z1[t+1] <- 1
      V2[t+1] <- v.reset
      fr1 <- c(fr1, t+1)
    }
  }
  return(list(V1, V2, Ps1, Ps2, Z1, Z2, fr1, fr2))
  })
}



parameters <- c(e.s = 0, e.l = -70, 
                #e.s : excitatory or inhibitory effect on neuron of a successful input AP
                #e.l : resting membrane potential of neuron
                v.th = -54, v.reset=-80, 
                #v.th : threshold voltage to initiate AP
                #v.reset : voltage reset when AP is fired
                tau.m = 20,tau.s = 10,
                #tau.m : 
                #tau.s : 
                rm.gs = 0.15, bigrm.ie = 18, 
                #rm.gs : 
                #bigrm.ie : 
                pmax =0.5
                #pmax : 
) 


#example excitatory vs inhibitory

time <- 400
coupled.1 <- integrate.fire(time, 0.1, parameters )

coupled.1[[1]][find_peaks(coupled.1[[1]], m=3)] <- 1
coupled.1[[2]][find_peaks(coupled.1[[2]], m=3)] <- 0

#par(mfrow=c(2,1), mar=c(1,4,2,2))
par(mfrow=c(1,1), mar=c(4,4,4,2))
plot(coupled.1[[1]], xlab="Time /ms", ylab="Voltage /mV", type="l", xaxt ="n", col="blue")
axis(side=1, at=seq(0, time*10, 500), labels=seq(0, time, 50))
segments(x0= -300, x1=-50, y0=coupled.1[[1]][1], y1=coupled.1[[1]][1],col="blue")
text(x=-500, y=coupled.1[[1]][1], srt=0, adj = 0, labels = as.character(round(coupled.1[[1]][1])), xpd = TRUE, col="blue") 
title(main="Coupled inhibitory neurons")
#plot(coupled.1[[2]], xlab="Time /ms", ylab="Voltage /mV", type="l", xaxt ="n", col="orange")
lines(coupled.1[[2]], col="orange")
segments(x0= -300, x1=-50, y0=coupled.1[[2]][1], y1=coupled.1[[2]][1],col="orange")
text(x=-500, y=coupled.1[[2]][1], srt=0, adj = 0, labels = as.character(round(coupled.1[[2]][1])), xpd = TRUE, col="orange") 

seq(-75, -53, 2)


#investigate initial voltage

investigate.initial.voltage <- function(time, h, Es){
  vs <- seq(-75, -53, 2)
  firing.rate.V1 <- matrix(0, nrow= length(vs), ncol=length(vs) )
  firing.rate.V2 <- matrix(0, nrow=length(vs), ncol=length(vs) )
  synched <- matrix(0, nrow=length(vs), ncol=length(vs) )
  for (v1 in 1:length(vs)){
    for (v2 in 1:length(vs)){
      parameters <- c(e.s = Es, e.l = -70, 
                      #e.s : excitatory or inhibitory effect on neuron of a successful input AP
                      #e.l : resting membrane potential of neuron
                      v.th = -54, v.reset=-80, 
                      #v.th : threshold voltage to initiate AP
                      #v.reset : voltage reset when AP is fired
                      tau.m = 20,tau.s = 10,
                      #tau.m : 
                      #tau.s : 
                      rm.gs = 0.15, bigrm.ie = 18, 
                      #rm.gs : 
                      #bigrm.ie : 
                      pmax =0.5
                      #pmax : 
      ) 
      coupled <- integrate.fire(time, h, parameters, v.one.init=vs[v1], v.two.init = vs[v2])
      print(c(vs[v1], vs[v2]))
      
      V1.peaks <- find_peaks(coupled[[1]])
      V2.peaks <- find_peaks(coupled[[2]])
      
      V1.isi <- ISI(V1.peaks, 1/h)
      V2.isi <- ISI(V2.peaks, 1/h)
      
      firing.rate.V1[v1,v2] <- V1.isi[[3]]
      firing.rate.V2[v1,v2] <- V2.isi[[3]]
      
      if (abs(V1.peaks[length(V1.peaks)] - V2.peaks[length(V2.peaks)]) <5){
        synched[v1,v2] <- TRUE
      } else{
        synched[v1,v2] <- FALSE
      }    
    }
  }
  return(list(synched, firing.rate.V1,firing.rate.V2))
}



init.volt.excit <- investigate.initial.voltage(4000, 0.1, -80)


vs <- seq(-75, -53, 2)
par(mar=c(4,4,4,2))
image(init.volt.excit[[1]], axes=FALSE, col = grey(seq(0, 1, length = 2)))
axis(1, at=seq(0,1,length.out = 12), label=as.character(vs) ) 
axis(2, at=seq(0,1,length.out = 12), label=as.character(vs) )
title(xlab="V1 starting voltage (mV)", ylab="V2 starting voltage (mV)", 
      main="Synchrony of firing in inhibitory \n coupled neurons after 4s")




#investigate tau of synapses


investigate.tau <- function(time, h, Es){
  vs <- seq(1,60,4)
  firing.rate.V1 <- matrix(0, nrow=length(vs), ncol=length(vs))
  firing.rate.V2 <- matrix(0, nrow=length(vs), ncol=length(vs))
  synched <- matrix(0, nrow=length(vs), ncol=length(vs))
  
  for (i in 1:length(vs)){
    for (j in 1:length(vs)){
      parameters <- c(e.s = Es, e.l = -70, 
                      #e.s : excitatory or inhibitory effect on neuron of a successful input AP
                      #e.l : resting membrane potential of neuron
                      v.th = -54, v.reset=-80, 
                      #v.th : threshold voltage to initiate AP
                      #v.reset : voltage reset when AP is fired
                      tau.m = vs[i],tau.s = vs[j],
                      #tau.m : 
                      #tau.s : 
                      rm.gs = 0.15, bigrm.ie = 18, 
                      #rm.gs : 
                      #bigrm.ie : 
                      pmax =0.5
                      #pmax : 
      ) 
      coupled <- integrate.fire(time, h, parameters)
      
      print(c(vs[i], vs[j]))
      
      V1.peaks <- find_peaks(coupled[[1]])
      V2.peaks <- find_peaks(coupled[[2]])
      
      V1.isi <- ISI(V1.peaks, 1/h)
      V2.isi <- ISI(V2.peaks, 1/h)

      firing.rate.V1[i,j] <- V1.isi[[3]]
      firing.rate.V2[i,j] <- V2.isi[[3]]
      
      if (abs(V1.peaks[length(V1.peaks)] - V2.peaks[length(V2.peaks)]) <5){
        synched[i,j] <- TRUE
      } else{
        synched[i,j] <- FALSE
      }    
    }
  }
  return(list(firing.rate.V1, firing.rate.V2, synched))
}


tau.investigation.inhib <- investigate.tau(4000, 0.1, -80)
tau.investigation.excit <- investigate.tau(4000, 0.1, 0)


#edit tau.investigation for choice of excit or inhib
vs <- seq(1,60,4)
par(mar=c(4,4,4,2))
image(tau.investigation[[3]], axes = FALSE, col = grey(seq(0, 1, length = 2)))
axis(1, at=seq(0,1,length.out = 15), label=as.character(vs) ) 
axis(2, at=seq(0,1,length.out = 15), label=as.character(vs) )
title(xlab="Tau_m (ms)", ylab="Tau_s (ms)", 
      main="Synchrony of firing in inhibitory \n coupled neurons after 4s")

persp(z=tau.investigation[[2]], theta=25, phi=5, scale=T, 
      ticktype = "detailed", 
      xlab="tau_m", 
      ylab="tau_s", 
      zlab="Firing rate", 
      main="Firing rate as altering tau_m and tau_s \n with inhibitory synapses")




#investigate strength of synapses

investigate.strength <- function(time, h, Es){
  vs <- seq(0, 1.5, 0.05)
  firing.rate.V1 <- integer(length(vs))
  firing.rate.V2 <- integer(length(vs))
  synched <- integer(length(vs))
  for (s in 1:length(vs)){
    parameters <- c(e.s = Es, e.l = -70, 
                    #e.s : excitatory or inhibitory effect on neuron of a successful input AP
                    #e.l : resting membrane potential of neuron
                    v.th = -54, v.reset=-80, 
                    #v.th : threshold voltage to initiate AP
                    #v.reset : voltage reset when AP is fired
                    tau.m = 20,tau.s = 10,
                    #tau.m : 
                    #tau.s : 
                    rm.gs = vs[s], bigrm.ie = 18, 
                    #rm.gs : 
                    #bigrm.ie : 
                    pmax =0.5
                    #pmax : 
    ) 
    coupled <- integrate.fire(time, h, parameters, v.one.init=-65, v.two.init = -67)
    print(vs[s])
    
    V1.peaks <- find_peaks(coupled[[1]])
    V2.peaks <- find_peaks(coupled[[2]])
    
    V1.isi <- ISI(V1.peaks, 1/h)
    V2.isi <- ISI(V2.peaks, 1/h)
    
    firing.rate.V1[s] <- V1.isi[[3]]
    firing.rate.V2[s] <- V2.isi[[3]]
    
    if (abs(V1.peaks[length(V1.peaks)] - V2.peaks[length(V2.peaks)]) < 2){
      synched[s] <- TRUE
    } else{
      synched[s] <- FALSE
    }    
    
  }
  return(list(synched, firing.rate.V1,firing.rate.V2))
}


strength.investigation.excit <- investigate.strength(4000, 0.1, 0)
strength.investigation.inhib <- investigate.strength(4000, 0.1, -80)


strengths <-  seq(0, 1.5, 0.05)
plot(strengths, strength.investigation.excit[[3]], 
     xlab=expression(paste("Strength of synapse", "(",r[m]*g[s], "value)")), 
     ylab="firing rate /Hz", 
     main="Firing rate of V1 and V2 \n with changing synapse strength \n in excitatory system", 
     pch=19, col="blue")
points(x= strengths,y=strength.investigation.excit[[2]], col="orange", pch=20)






## Networks

Mee <- 1.25
Mie <- 1
Mii <- -1
Mei <- -1
gamma.e <- -10
gamma.i <- 10
t.e <- 10

lambda.neg <- integer(200)
lambda.pos <- integer(200)
for (t.i in 1:200){
  lambda.neg[t.i] <- 0.5*((Mee-1)/t.e + (Mii-1)/t.i) - sqrt( as.complex(( (Mee-1)/t.e - (Mii-1)/t.i )^2 + (4*Mei*Mie)/(t.e*t.i) ) )
  lambda.pos[t.i] <- 0.5*((Mee-1)/t.e + (Mii-1)/t.i) + sqrt( as.complex(( (Mee-1)/t.e - (Mii-1)/t.i )^2 + (4*Mei*Mie)/(t.e*t.i) ) )
}
lambdas <- list(lambda.neg,lambda.pos)



simulation.network.ode <- function(t,state, parameters){
  with(as.list(c(state, parameters)), {
    dVe <- (-Ve + ifelse((Mee*Ve + Mei*Vi - gamma.e)>0, max(Mee*Ve + Mei*Vi - gamma.e), 0))/tau.e
    dVi <- (-Vi + ifelse((Mii*Vi + Mie*Ve - gamma.i)>0, max(Mii*Vi + Mie*Ve - gamma.i), 0))/tau.i
    return(list(c(dVe, dVi)))
  })
}

parameters <- c(Mee = 1.25, Mei = -1, gamma.e = -10, tau.e = 10,
                Mii = -1, Mie = 1, gamma.i = 10, tau.i = 85)
state <- c(Ve = 35, Vi = 15)
times <- seq(0, 1000, by = 0.1)
out <- ode(y = state, times = times, func = simulation.network.ode, parms = parameters)

par(mar=c(4,4,4,4), mfrow=c(1,1))
plot(out[,c(1,2)], xlab="Time /ms", xaxt='n' ,ylab="Firing rate /Hz", 
     main=expression(paste("Stabilising neuron firing rates with ", 
                           tau[I], " = 55 ms")),
     type="l", ylim=c(15,80), col="blue" )
axis(1, at=seq(0,1000,100), labels=seq(0,1000,100))
lines(out[,1],out[,3], col="orange")

par(mar=c(4,4,4,4), mfrow=c(1,1))
plot(out[,c(1,2)], xlab="Time /ms", xaxt='n' ,ylab="Firing rate /Hz", 
     main=expression(paste("Oscillating neuron firing rates with ", 
                           tau[I], " = 85 ms")),
     type="l", ylim=c(0,100), col="blue" )
axis(1, at=seq(0,1000,100), labels=seq(0,1000,100))
lines(out[,1],out[,3], col="orange")



simulation.network.phaseR <- function(t,y, params){
  Ve <- y[1]
  Vi <- y[2]
  Mee = params[1]
  Mei = params[2]
  gamma.e = params[3]
  tau.e = params[4]
  Mii = params[5]
  Mie = params[6]
  gamma.i = params[7]
  tau.i = params[8]
  dVe <- (-Ve + ifelse((Mee*Ve + Mei*Vi - gamma.e)>0, max(Mee*Ve + Mei*Vi - gamma.e), 0))/tau.e
  dVi <- (-Vi + ifelse((Mii*Vi + Mie*Ve - gamma.i)>0, max(Mii*Vi + Mie*Ve - gamma.i), 0))/tau.i
  dy <- numeric(2)
  dy[1] <- dVe
  dy[2] <-dVi 
  list(dy)
}

params <- c(1.25,-1,-10, 10, -1, 1, 10, 79)

par(pty="m")
sim.network.flowField <- flowField(simulation.network.phaseR,
                                 xlim = c(0,100),
                                 ylim = c(0, 50),
                                 parameters = params,
                                 arrow.type = "proportional",
                                 add = FALSE, 
                                 xlab=expression(paste(v[E], " /Hz")), 
                                 ylab=expression(paste(v[I], " /Hz")), 
                                 main=expression(paste("Phase plane analysis of stable system with ", tau[i], "=79ms")))

nullclines(simulation.network.phaseR,
           xlim = c(0,100),
           ylim = c(0, 100),
           parameters = params, 
           add=TRUE,
           legend=c(expression(paste(v[E])), expression(paste(v[I])) ) )



trajectory(simulation.network.phaseR, y0 = c(59.75,24.75), 
           tlim=c(0,10000),  tstep = 0.1,
           parameters = params)



kill<- 300


# Hopfield

Hopfield <- function(patterns, new.pattern, kill=0){
  weights <- t(patterns)%*%patterns
  diag(weights) <- NA
  k <- 0
  if (kill > 0){
    weights[sample(dim(weights)[1], kill), sample(dim(weights)[2], kill)] <- NA
    weights <- weights + t(weights)
  }
  repeat {
    i <- sample(dim(patterns)[2], 1, replace=T)
    prev.nodes <- new.pattern
    new.pattern[i] <- ifelse( sum(weights[,i]*new.pattern, na.rm = T) >= 0, 1, -1) 
    new.nodes <- new.pattern
    
    if (all.equal(new.nodes, prev.nodes) == T){
      k <- k+1
    }
    if (k >= 2){
      #print("2 iterations with no changes")
      break
    }
  
  #m <- matrix(new.pattern, sqrt(length(new.pattern)), sqrt(length(new.pattern)))
  #par(mar=c(0, 0, 0, 0))
  #image(m, axes = FALSE, col = grey(seq(0, 1, length = 2)))
  }
  return(list(new.pattern, weights))
}


mutate.generate <- function(pattern, mistakes = 3, N){
  patterns <- matrix(pattern, nrow=N, ncol=length(pattern))
  for (n in 1:N){
    for (i in sample(1:length(pattern), mistakes)){
      patterns[n, i] <- pattern[i]*-1
    }
  }
  return(patterns)
}



# par(mar=c(0, 0, 0, 0), pty="s")
# m <- matrix(patterns[1,], sqrt(length(patterns[1,])), sqrt(length(patterns[1,])))
# image(m, axes = FALSE, col = grey(seq(0, 1, length = 2)))
# 
# new.pattern <- mutate.generate(patterns[2,], mistakes = I^2/50, N=1)
# new.pattern.image <- matrix(new.pattern, sqrt(length(new.pattern)), sqrt(length(new.pattern)))
# image(new.pattern.image, axes = FALSE, col = grey(seq(0, 1, length = 2)))

check.hopf <- function(patterns, hopf.pattern, strict = FALSE){
  differences <- integer(dim(patterns)[1]) 
  overlap <- integer(dim(patterns)[1]) 
  if (strict == TRUE ){
    max.diff <- 0 
  } else {
    x <- as.integer(sum(patterns[1,]==1)/(sqrt(dim(patterns)[2])*4))
    max.diff <- ifelse(x==0, 1, x)
  }
  for (pattern in 1:(dim(patterns)[1])){
    differences[pattern] <- sum(patterns[pattern,] != hopf.pattern)
    overlap[pattern] <- as.numeric( ( patterns[pattern,]%*%as.vector(hopf.pattern) ) /dim(patterns)[2])
  }
  norm.diff <- differences/(dim(patterns)[2])
  best.match <- which.min(differences)
  learned <- ifelse(min(differences) <= max.diff, T, F)
  
  return(list(max(overlap),learned, best.match ))
}




I <- 20
num.patterns <- 10
mistakes <- I^2/2
new.mistakes <- as.integer(mistakes/30)

patterns <- mutate.generate(rep(-1, I^2), mistakes= mistakes , N= num.patterns )

new <- sample(1:num.patterns, 1)
new.pattern <- mutate.generate(patterns[new,], mistakes = ifelse(new.mistakes==0, 1, new.mistakes) , N=1)
#new.pattern <- patterns[new, ]
image(matrix(patterns[new,], I, I), axes = FALSE, col = grey(seq(0, 1, length = 2)))
image(matrix(new.pattern, I, I), axes = FALSE, col = grey(seq(0, 1, length = 2)))


hopf.temp <- Hopfield(patterns, new.pattern)
hopf.pattern <- hopf.temp[[1]]
image(matrix(hopf.pattern, I, I), axes = FALSE, col = grey(seq(0, 1, length = 2)))

test.check <- check.hopf(patterns, hopf.pattern, strict = F)

test.check

length(patterns[new,])



storage.capacity <- function(I, K, S, single.pattern.set=F, kill){
  # I: sqrt of number of nodes in network
  # S: sparseness of patterns
  # K: number of repetitions 
  # p: sequence of number of patters 
  if (single.pattern.set==T){
    p <- c(I^2*0.1)
  } else{
    p <- seq(1, I^2/2, by=I/2)
  }
  learned <- matrix(NA, nrow=length(p), ncol=K)
  overlap <- matrix(NA, nrow=length(p), ncol=K)
  match <- matrix(NA, nrow=length(p), ncol=K)
  for (num.patterns in 1:length(p) ){
    patterns <- mutate.generate(rep(-1, I^2), mistakes=S, N=p[num.patterns] )
    new <- sample(1:p[num.patterns], 1)
    new.mistakes <- as.integer(S/ (3*I))
    
    #new.pattern <- mutate.generate(patterns[new,] , 
    #                               mistakes = ifelse(new.mistakes==0, 1, new.mistakes),
    #                               N=1)
    new.pattern <- patterns[new, ]
    for (k in 1:K){
      #print(c( k, p[num.patterns], S))
      hopf.temp <- Hopfield(patterns, new.pattern, kill)
      hopf.pattern <- as.vector(hopf.temp[[1]])
      weights <- hopf.temp[[2]]
      
      l.temp <- check.hopf(patterns, hopf.pattern, strict=T )
      overlap[num.patterns, k] <- l.temp[[1]]
      learned[num.patterns, k] <- l.temp[[2]]
      
      if ( isTRUE(all.equal(l.temp[[3]], new)) == T) {
        match[num.patterns, k] <- T
      } else{
        match[num.patterns,k ] <- F
      }
    }
  }
  return(list(overlap, learned, match, weights))
}

test.storage.cap <- storage.capacity(20, 100, 200)


par(pty="s")
plot(apply(test.storage.cap[[2]], 1, sum), xaxt="n", 
     xlab="Number of patters/Size of network", 
     ylab="Percentage correct match", 
     main="Storage capacity of Hopfield network \n at sparcity = 50%", pch=20)
axis(1, at=seq(0,20,length.out=10), label=round(
  seq(1, 20^2/2, length.out=10)/400, 2) )
abline(v=7.5, col="lightblue")

storage.capacity.sparseness <- function(I){
  s <- seq(1, I^2, by=4*I)
  p <- seq(1, I^2/2, by=I/2)
  overlap <- matrix(NA, nrow=length(s), ncol= length(p) )
  learned <- matrix(NA, nrow=length(s), ncol=length(p) )
  match <- matrix(NA, nrow=length(s), ncol= length(p) )
  for (ones in 1:length(s) ){
    storage.res <- storage.capacity(I, K=100, S=s[ones])
    overlap[ones,] <- apply(storage.res[[1]], 1, sum)
    learned[ones, ] <- apply(storage.res[[2]], 1, sum)
    match[ones,] <- apply(storage.res[[3]], 1, sum)
    
  }
  return(list(overlap, learned, match ))
}

sparce.pattern.capacity <- storage.capacity.sparseness(20)


sparce.pattern.capacity[[2]]
persp(z=sparce.pattern.capacity[[2]], theta=130, phi=20, scale=T, 
      ticktype = "simple", 
      xlab="Sparseness", 
      ylab="Number of patters", 
      zlab="Percentage correct", 
      main="Percentage correct while altering sparseness \n and number of patterns learned")

plot(sparce.pattern.capacity[[2]][3,], type='l',
     xaxt="n", col="purple",
     xlab="Number of patters/Size of network", 
     ylab="Percentage correct match", 
     main="Storage capacity of Hopfield network \n while varying sparseness", pch=20)
axis(1, at=seq(0,20,length.out=10), label=round(
  seq(1, 20^2/2, length.out=10)/400, 2) )
lines(sparce.pattern.capacity[[2]][4,], col="green")
lines(sparce.pattern.capacity[[2]][5,], col="orange")
lines(sparce.pattern.capacity[[2]][2,], col="cyan")

legend(x=13, y=99, legend="Sparseness \n Cyan: 20% \n Purple: 40% \n Green: 60% \n Orange: 80%", bty = "n")



robustness <- function(I){
  r <- seq(0, I^2, I)
  p <- seq(1, I^2/2, by=I/2)
  overlap <- matrix(NA, length(r), length(p))
  learned <- matrix(NA, length(r), length(p))
  match <- matrix(NA, length(r), length(p))
  #weights <- list()

  for (k in 1:length(r) ){
    print(r[k])
    storage.res <- storage.capacity(I, K=100, S=I^2/2, single.pattern.set = F, kill=r[k])
    overlap[k,] <- apply(storage.res[[1]], 1, sum)
    learned[k,] <- apply(storage.res[[2]], 1, sum)
    match[k,] <- apply(storage.res[[3]], 1, sum)
    #weights[[k]] <- storage.res[[4]]
    }
  return(list(overlap,learned, match,weights ))
}

storage.res <- storage.capacity(I, K=100, S=I^2/2, single.pattern.set = T, kill=390)
apply(storage.res[[2]], 1, sum)

test.single.robust <- robustness(20)

test.single.robust[[2]]

par(pty="s")
plot(test.single.robust[[2]], xaxt="n", 
     xlab="Percentage weights lost", 
     ylab="Percentage correct match", 
     main="Robustness of Hopfield network \n with N/I = 0.1 and sparcity = 50%", pch=20)
axis(1, at=1:21, label=seq(0, I^2, I)/400*100 )
abline(v=8.5, col="lightblue")


persp(z=test.single.robust[[2]], theta=70, phi=20, scale=T, 
      ticktype = "simple", 
      xlab="Percentage weights lost", 
      ylab="Number of patters", 
      zlab="Percentage correct", 
      main="Percentage correct while altering number of weights lost \n and number of patterns learned")

plot(test.single.robust[[2]][1,], type='l', ylim=c(0,100),
     xaxt="n", col="purple",
     xlab="Number of patters/Size of network", 
     ylab="Percentage correct match", 
     main="Robustness of Hopfield network  while varying \n number of patterns and loss of weights", pch=20)
axis(1, at=seq(0,20,length.out=10), label=round(
  seq(1, 20^2/2, length.out=10)/400, 2) )
lines(test.single.robust[[2]][5,], col="green")
lines(test.single.robust[[2]][10,], col="navy")
lines(test.single.robust[[2]][15,], col="orange")
lines(test.single.robust[[2]][20,], col="cyan")

legend(x=0, y=30, legend="Loss of weights \n Purple: 0% Green: 20% \n Navy: 45% Orange: 70%, Cyan: 100%", bty = "n")

test.single.robust[[2]][15,]




### new method

improved.Hopfield <- function(patterns, new.pattern, L){
  weights <- t(patterns)%*%patterns
  diag(weights) <- 0
  
  t <- patterns
  t[t==-1] <- 0
  for (l in 1:L){
    diag(weights) <- 0
    activations <- patterns %*% weights
    outputs <- sigmoid(activations)
    errors <- t - outputs
    gradients <- t(patterns) %*% errors
    gradients <- gradients + t(gradients)
    
    weights <- weights + 0.01*(gradients - 0.1*weights) 
  }
}








######################


# nullclines <- function(){
#   Ve.xnull <- integer(101)
#   Vi.xnull <- integer(101)
#   Ve.ynull <- integer(101)
#   Vi.ynull <- integer(101)
#   Mii <- -1
#   Mee <- 1.25
#   Mei <- -1
#   Mie <- 1
#   gamma.i <- 10
#   gamma.e <- -10
# 
#   #x/Ve nullcline
#   for (e in 0:100){
#     Vi.xnull[e+1] <- e*((1-Mee)/Mei) + gamma.e/Mei
#   }
#   for (i in 0:100){
#     Ve.xnull[i+1] <- i*(Mei/(1-Mee)) - gamma.e/(1-Mee)
#   }
# 
#   #ynullcline
#   for (e in 0:100){
#     Vi.ynull[e+1] <- (e*Mie)/(1-Mii) - gamma.i/(1-Mii)
#   }
#   for (i in 0:100){
#     Ve.ynull[i+1] <- i*((1-Mii)/Mie) + gamma.i/Mie
#   }
# 
#   return(list(Vi.xnull,Ve.xnull, Vi.ynull, Ve.ynull))
# }
# 
# p <- nullclines()
# 
# plot(x = 0:100, y = p[[1]], type="l", col='red', xlim=c(0,100), ylim=c(0, 100))
# lines(x=p[[2]], y=0:100, col='red')
# lines(x=0:100, y=p[[3]], col='blue')
# lines(x=p[[4]], y=0:100, col='blue')


