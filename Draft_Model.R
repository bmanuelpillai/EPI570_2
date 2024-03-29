dSm <- (-lambda.m.uv.sy * s.m.num) + (-lambda.m.uv.ay * s.m.num) + (1-(omega*chi*mu*num)) + delta(i.m.sy.num) + delta(i.m.ay.num)

dIm <- (lambda.m.v.sy * v.m.num) + (lambda.m.v.ay * v.m.num) + (lambda.m.uv.sy * s.m.num) - delta(i.m.sy.num) - delta(i.m.ay.num) 

dVm <- (-lambda.m.v.sy * v.m.num) + (-lambda.m.v.ay * v.m.num) + + (omega*chi*mu*num)


dSf <- (-lambda.f.uv.sy * s.f.num) + (-lambda.f.uv.ay * s.f.num) + (1-(omega*chi*mu*num)) + delta(i.f.sy.num) + delta(i.f.ay.num)

dIf <- (lambda.f.v.sy * v.f.num) + (lambda.f.v.ay * v.f.num) + (lambda.f.uv.sy * s.f.num) - delta(i.f.sy.num) - delta(i.f.ay.num) 

dVf <- (-lambda.f.v.sy * v.f.num) + (-lambda.f.v.ay * v.f.num) + (omega*chi*mu*num)
