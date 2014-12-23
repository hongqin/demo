rm(list=ls())
source("/Users/hongqin/lib/R/lifespan.r")

##two parameter Gompertz model 
N=1000
I =0.005
G = 0.10
t= seq(1, 100, by=1)
s = G.s(c(I,G,0), t)
plot( s ~ t)
mu = I * exp(G*t)
plot ( s*mu ~ t )
pmf = s * mu  # prob mass function
pmf = pmf / sum(pmf)
summary ( pmf * N)
hist( round( pmf * N) )
lifespanT = rep( t, round(pmf*N))
hist(lifespanT, br=20)
summary(lifespanT)


##### log likelihood function, simple gompertz mortality model
#s = exp( (I/G) *(1 - exp(G* my.data)) )  ;
#m = I * exp( G * my.data );   
llh.gompertz.single.run <- function( IG, lifespan, trace=0 ) {
  #print(IG)
  I = IG[1]; G = IG[2]; 
  my.data = lifespan[!is.na(lifespan)];
  log_s = (I/G) *(1 - exp(G* my.data))
  #if( I< 0 ) { I = 1E-10 }
  log_m = log(I) +  G * my.data ; 
  my.lh = sum(log_s)  + sum(log_m);
  if (trace) {
    print (IG ); #trace the convergence
  }
  ret = - my.lh # because optim seems to maximize 
}

lifespan = sample(lifespanT, 30)
shat = calculate.s( lifespan )
plot(shat$s ~ shat$t )

ret1a = optim ( c(I,G)*2, fn=llh.gompertz.single.run, lifespan=lifespan, lower=c(1E-10, 1E-5), method="L-BFGS-B" );
ret1a$par / c(I,G)

ret1a = optim ( c(I, G), fn=llh.gompertz.single.run, lifespan=lifespan, lower=c(1E-10, 1E-5), method="L-BFGS-B" );
ret1a$par / c(I,G)

ret1b = optim ( c(I, G)*c(0.01, 1.5), fn=llh.gompertz.single.run, lifespan=lifespan, lower=c(1E-10, 1E-5), method="L-BFGS-B" );
ret1b$par / c(I,G)

ret1c = optim ( c(I,G)*c(0.01,1.5), fn=llh.gompertz.single.run, lifespan=lifespan, trace=1  );
ret1c$par / c(I, G)

ret1d = optim ( c(0.05, 0.15), fn=llh.gompertz.single.run, lifespan=lifespan, trace=1  );
ret1d$par / c(I,G)

ret1e = optim ( c(0.05, 0.1), fn=llh.gompertz.single.run, lifespan=lifespan  );
ret1e$par / ret1d$par

#################
## 3  parameter Gompertz model 
N=1000
I =0.005; G = 0.10; M = I
t= seq(1, 100, by=1)
s = GM.s(c(I,G,0), t)
plot( s ~ t)
mu = I * exp(G*t) + M
plot ( s*mu ~ t )
pmf = s * mu  # prob mass function
pmf = pmf / sum(pmf)
summary ( pmf * N)
hist( round( pmf * N) )
lifespanT = rep( t, round(pmf*N))
hist(lifespanT, br=20)
summary(lifespanT)


lifespan = sample(lifespanT, 30)
shat = calculate.s( lifespan )
plot(shat$s ~ shat$t )

ret1a = optim ( c(I,G, M)*2, fn=llh.GM.single.run, lifespan=lifespan, lower=c(1E-10, 1E-5, 0), method="L-BFGS-B" );
ret1b = optim ( c(I, G)*2, fn=llh.G.single.run, lifespan=lifespan, lower=c(1E-10, 1E-5), method="L-BFGS-B" );
ret1a$par / c(I,G, M)
ret1b$par / c(I,G)
#ignore larger M can inflate R

ret1a = optim ( par = c(I, G, M*2), fn=llh.GM.single.run, lifespan=lifespan, lower=c(1E-10, 1E-5, 1E-15), method="L-BFGS-B" );
ret1b = optim ( par = c(I*0.1, G), fn=llh.G.single.run, lifespan=lifespan, lower=c(1E-10, 1E-5), method="L-BFGS-B" );
ret1a$par / c(I,G, M)
ret1b$par / c(I,G)

