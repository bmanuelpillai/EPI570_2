
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

param <- param.dcm(omega = ?, chi = ?, mu = ?, delta = ?,
                   ,ce.m.ay = ?, nu = ?, tau.uv = ?, tau.v = ?)

init <- init.dcm(s.m.num = ?, i.m.sy.num = ?, i.m.ay.num = ?, v.m.num = ?, 
                 s.f.num = ?, i.f.sy.num = ?, i.f.ay.num = ?,v.f.num = ?,)

control <- control.dcm(nsteps = ?, new.mod = AoN)

