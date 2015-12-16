#example from GetoptLong

# Rscript test-getoptlong.R --number 4 --verbose T
# Rscript test-getoptlong.R -n 4 --verbose T
library(GetoptLong)

cutoff = 0.05
GetoptLong(c(
  "number|n=i", "Number of items, integer, mandatory option",
  "cutoff|c=f", "cutoff to filter results, optional, default (0.05)",
  "verbose|v", "print messages"
))

print(c("input:", 'number'=number, 'cutoff'=cutoff, 'verbose'=verbose))


