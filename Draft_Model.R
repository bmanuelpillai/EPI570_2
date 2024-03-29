# males
dSm <- (-lambda.m.uv.sy * s.m.num) + (-lambda.m.uv.ay * s.m.num) + (1-(omega*chi*mu*num)) + delta(i.m.sy.num) + delta(i.m.ay.num)

dIm <- (lambda.m.v.sy * v.m.num) + (lambda.m.v.ay * v.m.num) + (lambda.m.uv.sy * s.m.num) + (lambda.m.uv.ay * s.m.num) - delta(i.m.sy.num) - delta(i.m.ay.num) 

dVm <- (-lambda.m.v.sy * v.m.num) + (-lambda.m.v.ay * v.m.num) + + (omega*chi*mu*num)

# females
dSf <- (-lambda.f.uv.sy * s.f.num) + (-lambda.f.uv.ay * s.f.num) + (1-(omega*chi*mu*num)) + delta(i.f.sy.num) + delta(i.f.ay.num)

dIf <- (lambda.f.v.sy * v.f.num) + (lambda.f.v.ay * v.f.num) + (lambda.f.uv.sy * s.f.num) + (lambda.f.uv.ay * s.f.num) - delta(i.f.sy.num) - delta(i.f.ay.num) 

dVf <- (-lambda.f.v.sy * v.f.num) + (-lambda.f.v.ay * v.f.num) + (omega*chi*mu*num)



#defining contact rates

# We have to think about how symptomatic status impacts force of infection. Changing contact rate would be difficult
ce.f <- (c.m * m.num) / f.num

#defining lambdas

#For Males
lambda.m.uv.sy <- (tau.uv.sy) * (ce.m) * (i.m.sy.num / m.uv.num)
lambda.m.uv.ay <- (tau.uv.ay) * (ce.m) * (i.m.ay.num / m.uv.num)

lambda.m.v.sy <- (tau.v.sy) * (ce.m) * (i.m.sy.num / m.v.num)
lambda.m.v.ay <- (tau.v.ay) * (ce.m) * (i.m.ay.num / m.v.num)


#For Females
lambda.f.uv.sy <- (tau.uv.sy) * (ce.f) * (i.f.sy.num / f.uv.num)
lambda.f.uv.ay <- (tau.uv.ay) * (ce.f) * (i.f.ay.num / f.uv.num)

lambda.f.v.sy <- (tau.v.ay) * (ce.f) * (i.f.sy.num / f.v.num)
lambda.f.v.ay <- (tau.v.sy) * (ce.f) * (i.f.ay.num / f.v.num)

