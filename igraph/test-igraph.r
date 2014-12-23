source("http://bioconductor.org/biocLite.R")
biocLite("GO.db")
library(GO.db)
library(igraph)
BP <- toTable(GOBPPARENTS)
CC <- toTable(GOCCPARENTS)
MF <- toTable(GOMFPARENTS)
g <- graph.data.frame( rbind(BP,CC,MF) )


terms <- toTable(GOTERM)[,2:5]
terms <- terms[ !duplicated(terms[,1]), ]
g <- graph.data.frame( rbind(BP,CC,MF), vertices=terms )

g.d = degree(g)

#continental US demo
#list.files(path="continental_states")
#states = read.table( "continental_states/states.tab", colClass=c("character","character",NA))

require(RCurl)
URL = "https://docs.google.com/spreadsheet/pub?key=0ArLJZixvTlU7dDI5c3dnTzdRX1dndzBORUk1UC0wYlE&single=true&gid=0&output=csv"
states = read.csv(textConnection(getURL(URL)), colClass = c("character", "character"))

g = graph.data.frame(states, directed=F)
g.degree = degree(g)
g.degree [g.degree == max(g.degree)] #TN and MO have 8 bordering states

g.shortestpath.m = shortest.paths(g)
str(g.shortestpath.m)
sorted.names = sort( rownames(g.shortestpath.m) )
gsm = g.shortestpath.m[, sorted.names]
gsm = gsm[sorted.names, ]


