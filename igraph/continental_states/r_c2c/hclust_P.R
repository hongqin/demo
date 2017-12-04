 library(cluster);

 layout _ matrix(c (0.1, 0.5, 0.2, 0.8,
                 0.5, 0.9,   0.2, 0.8,
                 ),
                 nrow=2, ncol=4, byrow=T);
 close.screen(all+T);
 split.screen(figs=layout);

 cnum  _ 3;  	# the number of categories, i.e., the matrix dimension

# zt _ read.table("_Z._c2c_out3.10K.dat", header=T);
 zm0 _ matrix(scan("Z.dat", n=cnum*cnum, skip=1),cnum, cnum, byrow= TRUE)  

 cat.tb _ read.table("zlabels.dat", as.is=F, sep ="\t" );

 xx _ seq(1:cnum);
 yy  _ seq(1:cnum);

 pm.up _  	pnorm(zm0, mean=0, sd=1, lower.tail = F);
 pm.low _  	pnorm(zm0, mean=0, sd=1, lower.tail = T);

# convert p-value matrix to distance matrix-class and hclust()
 dx.up <-as.dist(pm.up);
 attr(dx.up, "Labels") _ cat.tb$V1 ;
 hc.up <- hclust(dx.up, "ave");
 screen(1);
 plot(hc.up, main="Interaction Significance Tree",
	xlab="average-linkage",ylab="p value difference" );


 dx.low <-as.dist(pm.low);
 attr(dx.low, "Labels") _ cat.tb$V1 ; 
 hc.low <- hclust(dx.low, "ave");     
 screen(2);
 plot(hc.low, main="Separation Significance Tree",
    xlab="average-linkage",ylab="p value difference" );        



q('no')

