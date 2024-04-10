

#1. Track Population Sizes
m.num = s.m.ue + s.m.e.num + i.m.s.a.num + i.m.s.s.num + v.m.ue + v.m.e.num + i.m.v.a.num + i.m.v.s.num + s.m.e.num
f.num = s.f.ue + s.f.e.num + i.f.s.a.num + i.f.s.s.num + v.f.ue + v.f.e.num + i.f.v.a.num + i.f.v.s.num + s.f.e.num
num = f.num + m.num

#2. Define Contact Rate
ce.f <- (ce.m * m.num) / f.num

#3. Define Lambda
  # For Susceptible Males
lambda.m.s.a <- (tau.s) * (ce.m) * (i.f.num / f.num)

  # For Susceptible Females
lambda.f.s.a <-  (tau.s) * (ce.f) * (i.m.num / m.num)
  
  # For Vaccinated Males
lambda.m.v.a <- (tau.v) * (ce.m) * (i.f.num / f.num)

# For Vaccinated Females
lambda.f.v.a <-  (tau.v) * (ce.f) * (i.m.num / m.num)


#4. Define Additional parameters
  mu <- 1/lifespand
  delta <- 1/dur.inf

#5. Define Differential Equations
  
#zeta is proportion of asympomatic that become symptomatic at each timestep, kappa is the proportion of susceptible exposed that become vaccinated exposed at each timestep
# MALES
dS.m.ue <- (1-(0.5 * omega * chi * mu * num)) - ((1/8) * s.m.ue)

dS.m.e <- ((1/8) * s.m.ue) - (lambda.m.s.a * s.m.e.num) + (delta.s.a * i.m.s.a.num) + (delta.s.s * i.m.s.s.num) - (kappa.m * s.m.e.num)

dI.m.s.a <- (lambda.m.s.a * s.m.e.num) - (zeta.s * i.m.s.a.num) - (delta.s.a * i.m.s.a.num)

dI.m.s.s <- (zeta.s * i.m.s.a.num) - (delta.s.s * i.m.s.s.num)


dV.m.ue <- (0.5 * omega * chi * mu * num) - ((1/8) * v.m.ue)

dV.m.e <- ((1/8) * v.m.ue) - (lambda.m.v.a * v.m.e.num) + (delta.v.a * i.m.v.a.num) + (delta.v.s * i.m.v.s.num) + (kappa.m * s.m.e.num)

dI.m.v.a <- (lambda.m.v.a * v.m.e.num) - (zeta.v * i.m.v.a.num) - (delta.v.a * i.m.v.a.num)

dI.m.v.s <- (zeta.v * i.m.v.a.num) - (delta.v.s * i.m.v.s.num)



# FEMALES
dS.f.ue <- (1-(0.5 * omega * chi * mu * num)) - ((1/8) * s.f.ue)

dS.f.e <- ((1/8) * s.f.ue) - (lambda.f.s.a * s.f.e.num) + (delta.s.a * i.f.s.a.num) + (delta.s.s * i.f.s.s.num) - (kappa.f * s.f.e.num)

dI.f.s.a <- (lambda.f.s.a * s.f.e.num) - (zeta.s * i.f.s.a.num) - (delta.s.a * i.f.s.a.num)

dI.f.s.s <- (zeta.s * i.f.s.a.num) - (delta.s.s * i.f.s.s.num)


dV.f.ue <- (0.5 * omega * chi * mu * num) - ((1/8) * v.f.ue)

dV.f.e <- ((1/8) * v.f.ue) - (lambda.f.v.a * v.f.e.num) + (delta.v.a * i.f.v.a.num) + (delta.v.s * i.f.v.s.num) + (kappa.f * s.f.e.num)

dI.f.v.a <- (lambda.f.v.a * v.f.e.num) - (zeta.v * i.f.v.a.num) - (delta.v.a * i.f.v.a.num)

dI.f.v.s <- (zeta.v * i.f.v.a.num) - (delta.v.s * i.f.v.s.num)

# 6.Outputs
list(c(dS.m.ue, dS.m.e, dI.m.s.a, dI.m.s.s, 
       dV.m.ue, dV.m.e, dI.m.v.a, dI.m.v.s,
       si.m.flow = (lambda.m.s.a * s.m.e.num),
       vi.m.flow = (lambda.m.v.a * v.m.e.num),
       dS.f.ue, dS.f.e, dI.f.s.a, dI.f.s.s, 
       dV.f.ue, dV.f.e, dI.f.v.a, dI.f.v.s,
       si.f.flow = (lambda.f.s.a * s.f.e.num),
       vi.f.flow = (lambda.f.v.a * v.f.e.num)))
       

       
       


