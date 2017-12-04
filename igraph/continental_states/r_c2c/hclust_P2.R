# here I experiment with different clustering methods.
 library(cluster);

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

 screen(1);
 hc.up <- hclust(dx.up, "single");
 plot(hc.up, main="Interaction Significance Tree, single-linkage",
	xlab="average-linkage",ylab="p value difference" );

 screen(2);
 hc.up <- hclust(dx.up, "ave");
 plot(hc.up, main="Interaction Significance Tree, average linkage",
	xlab="average-linkage",ylab="p value difference" );

 screen(3);
 hc.up <- hclust(dx.up, "complete");
 plot(hc.up, main="Interaction Significance Tree, complete linkage",
	xlab="average-linkage",ylab="p value difference" );


 dx.low <-as.dist(pm.low);
 attr(dx.low, "Labels") _ cat.tb$V1 ; 

 screen(4);
 hc.low <- hclust(dx.low, "single");     
 plot(hc.low, main="Separation Significance Tree, single linkage",
    xlab="average-linkage",ylab="p value difference" );        

 screen(4);
 hc.low <- hclust(dx.low, "ave");     
 plot(hc.low, main="Separation Significance Tree, average linkage",
    xlab="average-linkage",ylab="p value difference" );        

 screen(4);
 hc.low <- hclust(dx.low, "complete");     
 plot(hc.low, main="Separation Significance Tree, complete linkage",
    xlab="average-linkage",ylab="p value difference" );        



q('no')

