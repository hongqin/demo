 layout _ matrix(c (0, 0.8, 0, 1,
                 0.8, 1,   0, 1,
                 ),
                 nrow=2, ncol=4, byrow=T);
 close.screen(all+T);
 split.screen(figs=layout);


 newColHSVV _ function(h1, h2, n) {
    nside _ (n/2) - 1;
    s _ 1;
    ret _ c( hsv(h1, s, seq(1, 1/nside, length=nside), gamma=1),
             "#000000",
             hsv(h2, s, seq(1/nside, 1, length=nside), gamma=1) );
 };

 cnum  _ 3;  	# the number of categories, i.e., the matrix dimension

 zm _ matrix(scan("Z.dat", n=cnum*cnum, skip=1),cnum, cnum, byrow= TRUE)    

 zcutoff _ 10;
 zm[zm>zcutoff] _ zcutoff;		# level the highest outliers
 zm[zm< -zcutoff] _ -zcutoff;		# level the highest outliers

 zmax _ max(zm);
 xx _ seq(1:cnum);		# this should be done automatically
 yy  _ seq(1:cnum);

 breaks _ seq(-zmax, zmax, 1);
 colmap _ newColHSVV(1/3, 0, length(breaks));   # green-red gradient
 zlim _ c(-zmax-0.5, zmax);

# do image
 screen(1);
 image(xx, yy, zm[1:cnum, 1:cnum], zlim, col = colmap, axes=FALSE, 
   main="Z-score profile of interaction among US continental states",
   xlab="isotemporal categories", ylab="isotemporal categories");  

# do labels
 cat.tb _ read.table("zlabels.dat", as.is=TRUE );
 xlab _ paste( cat.tb$V1, sep=" ");
 axis(1, at=xx, labels=xlab);
 axis(2, at=yy, labels=xlab);

# now generate the scale-panel for Z-score
 midpnt _ seq( -zmax, zmax, 1 )
 midpnt _ matrix(midpnt, nrow=1, ncol=length(midpnt));
 screen(2);
 image(1:1, as.vector(midpnt), midpnt, zlim,  axes=FALSE,
       col=colmap, ylab="Z-score", xlab="");
 ylab _ round( as.vector(midpnt), digits=1);
 axis(2, at=as.vector(midpnt), labels=ylab);

# at this step, the image for Z-scores is done

 quit("no");
