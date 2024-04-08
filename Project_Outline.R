## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

library(pacman)
library(EpiModel)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(reshape2)


## ----3 Vaccination of HPV-----------------------------------------------------
Vaccination <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    # 1.Track the population sizes
        m.num = s.m.num + i.m.num + v.m.num
        f.num = s.f.num + i.f.num + v.f.num
        num = m.num + f.num

    
    # 2. Define Contact Rate based on Sex.Define contact rate for females based on the contact rate for males.
        ce.f <- (ce.m * m.num) / f.num
      # We are isolating to only heterosexual contacts.
      
    
   # 3.Define lambda
        #For Males
        lambda.m <- (tau) * (ce.m) * (i.f.num / f.num)

        #For Females
        lambda.f <- (tau) * (ce.f) * (i.m.num / m.num)

  # 4. Define additional parameters
    mu <- 1/lifespan
    delta <- 1/dur.inf
    
    # 5.Write six differential equations
    
# MALES
  dSm <- (-lambda.m * s.m.num) + ((1 - omega.m) * (mu * m.num)) + (delta * i.m.num)  - (mu * s.m.num)

  dIm <- (lambda.m * v.m.num * (1-psi))  + (lambda.m * s.m.num) - (delta * i.m.num)  - (mu * i.m.num)
  
  dVm <- (-lambda.m * v.m.num * (1-psi)) + (omega.m * mu * m.num) - (mu * v.m.num)

# FEMALES
  dSf <- (-lambda.f * s.f.num) + ((1 - omega.f)* (mu * f.num)) + (delta * i.f.num)  - (mu * s.f.num)

  dIf <- (lambda.f * v.f.num * (1-psi))  + (lambda.f * s.f.num) - (delta * i.f.num)  - (mu * i.f.num)
  
  dVf <- (-lambda.f * v.f.num * (1-psi)) + (omega.f * mu * f.num) - (mu * v.f.num)
    
    # 6.Outputs
    list(c(dSm, dIm,dVm,
            dSf, dIf, dVf,
            si.m.flow = (lambda.m * s.m.num),
            si.f.flow = (lambda.f * s.f.num),
            vi.m.flow = (lambda.m * v.m.num * (1-psi)),
            vi.f.flow = (lambda.f * v.f.num * (1-psi))
       ))
  })
}


param <- param.dcm(omega.m = seq(0.1, 0.4, by = 0.1), omega.f = 0.40, lifespan = (11*365),
                   psi = 0.90, dur.inf = (1.5*365), ce.m = (1.8/365), tau = 0.74)

init <- init.dcm(s.m.num = 3800000, i.m.num = 1300000, v.m.num = 304000,
                 s.f.num = 3000000, i.f.num = 2000000, v.f.num = 240000,
                 si.m.flow = 0, si.f.flow = 0,
                 vi.m.flow = 0,vi.f.flow = 0)

control <- control.dcm(nsteps = (11*365), new.mod = Vaccination)

mod <- dcm(param, init, control)
mod

test <- as.data.frame(mod)
View(test)


## -----------------------------------------------------------------------------
mod <- mutate_epi(mod, m.num = s.m.num + i.m.num + v.m.num, f.num = s.f.num + i.f.num + v.f.num)
mod <- mutate_epi(mod, m.prev = i.m.num / m.num, f.prev = i.f.num / f.num)


par(mfrow = c(2,2))
plot(mod, y = "m.prev", main = "Proportional Prevalence for Males\nomega.m = 0.1-0.4", legend = "full", ylim = c(0,0.5))
plot(mod, y = "f.prev", main = "Proportional Prevalence for Females\nomega.m = 0.1-0.4", legend = "full", ylim = c(0,0.5))
plot(mod, y = "si.m.flow", main = "Incidence in Unvaccinated Males\nomega.m = 0.1-0.4", legend = "full")
plot(mod, y = "si.f.flow", main = "Incidence in Unvaccinated Females\nomega.m = 0.1-0.4", legend = "full")


## -----------------------------------------------------------------------------
param2 <- param.dcm(omega.m = 0.1, omega.f = seq(0.40,0.80, by = 0.2), lifespan = (11*365),
                   psi = 0.90, dur.inf = (1.5*365), ce.m = (1.8/365), tau = 0.74)

init2 <- init.dcm(s.m.num = 3800000, i.m.num = 1300000, v.m.num = 304000,
                 s.f.num = 3000000, i.f.num = 2000000, v.f.num = 240000,
                 si.m.flow = 0, si.f.flow = 0,
                 vi.m.flow = 0,vi.f.flow = 0)

control2 <- control.dcm(nsteps = (11*365), new.mod = Vaccination)

mod2 <- dcm(param2, init2, control2)
mod2

mod2 <- mutate_epi(mod2, m.num = s.m.num + i.m.num + v.m.num, f.num = s.f.num + i.f.num + v.f.num)
mod2 <- mutate_epi(mod2, m.prev = i.m.num / m.num, f.prev = i.f.num / f.num)


par(mfrow = c(2,2))
plot(mod2, y = "m.prev", main = "Proportional Prevalence for Males\nomega.f = 0.4,0.6,0.8", legend = "full", ylim = c(0,0.5))
plot(mod2, y = "f.prev", main = "Proportional Prevalence for Females\nomega.f = 0.4,0.6,0.8", legend = "full", ylim = c(0,0.5))
plot(mod2, y = "si.m.flow", main = "Incidence in Unvaccinated Males\nomega.f = 0.4,0.6,0.8", legend = "full")
plot(mod2, y = "si.f.flow", main = "Incidence in Unvaccinated Females\nomega.f = 0.4,0.6,0.8", legend = "full")


## -----------------------------------------------------------------------------
m.num = s.m.num + i.m.sy.num+ i.m.ay.num + v.m.num
f.num = s.f.num + i.f.sy.num+ i.f.ay.num + v.f.num
num = m.num + f.num

# males
dSm <- (-lambda.m.uv.sy * s.m.num) + (-lambda.m.uv.ay * s.m.num) + (1-(0.5*omega*chi*mu*num)) + delta(i.m.sy.num) + delta(i.m.ay.num)

dIm <- (lambda.m.v.sy * v.m.num) + (lambda.m.v.ay * v.m.num) + (lambda.m.uv.sy * s.m.num) + (lambda.m.uv.ay * s.m.num) - delta(i.m.sy.num) - delta(i.m.ay.num) 

dVm <- (-lambda.m.v.sy * v.m.num) + (-lambda.m.v.ay * v.m.num) + (0.5*omega*chi*mu*num)

# females
dSf <- (-lambda.f.uv.sy * s.f.num) + (-lambda.f.uv.ay * s.f.num) + (1-(0.5*omega*chi*mu*num)) + delta(i.f.sy.num) + delta(i.f.ay.num)

dIf <- (lambda.f.v.sy * v.f.num) + (lambda.f.v.ay * v.f.num) + (lambda.f.uv.sy * s.f.num) + (lambda.f.uv.ay * s.f.num) - delta(i.f.sy.num) - delta(i.f.ay.num) 

dVf <- (-lambda.f.v.sy * v.f.num) + (-lambda.f.v.ay * v.f.num) + (0.5*omega*chi*mu*num)



#defining contact rates

# We are altering the contact rate between symptomatic and asymptomatic by nu.
# we are keeping the transmission probability between sympto and asympt the same aka tau
ce.f.ay <- (ce.m.ay * m.num) / f.num
ce.f.sy <- ce.f.ay * nu
ce.m.sy <- ce.m.ay * nu
#defining lambdas

#For Males
lambda.m.uv.ay <- (tau.uv) * (ce.m.ay) * (i.m.ay.num / m.uv.num)
lambda.m.uv.sy <- (tau.uv) * (ce.m.sy) * (i.m.sy.num / m.uv.num)

lambda.m.v.ay <- (tau.v) * (ce.m.ay) * (i.m.ay.num / m.v.num)
lambda.m.v.sy <- (tau.v) * (ce.m.sy) * (i.m.sy.num / m.v.num)


#For Females
lambda.f.uv.ay <- (tau.uv) * (ce.f.ay) * (i.f.ay.num / f.uv.num)
lambda.f.uv.sy <- (tau.uv) * (ce.f.sy) * (i.f.sy.num / f.uv.num)

lambda.f.v.ay <- (tau.v) * (ce.f.ay) * (i.f.ay.num / f.v.num)
lambda.f.v.sy <- (tau.v) * (ce.f.sy) * (i.f.sy.num / f.v.num)

#Output
list(c(dSm,dIm,dVm,
       dSf,dIf,dVf,
       si.m.flow = (lambda.m.uv.sy * s.m.num + lambda.m.uv.ay * s.m.num),
       si.f.flow = (lambda.f.uv.sy * s.f.num + lambda.f.uv.ay * s.f.num)
       ))

# param <- param.dcm(omega = ?, chi = ?, mu = ?, delta = ?,
#                    ,ce.m.ay = ?, nu = ?, tau.uv = ?, tau.v = ?)
# 
# init <- init.dcm(s.m.num = ?, i.m.sy.num = ?, i.m.ay.num = ?, v.m.num = ?, 
#                  s.f.num = ?, i.f.sy.num = ?, i.f.ay.num = ?,v.f.num = ?,)
# 
# control <- control.dcm(nsteps = ?, new.mod = AoN)

