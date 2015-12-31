library(foreach)
library(doMC)
registerDoMC(2)
x = foreach(i=1:(3*2)) %dopar% sqrt(i)
x

x2 = foreach(i=1:(3*2), .combine='cbind') %dopar% sqrt(i)
x2
