#' @importFrom colorspace hex2RGB mixcolor polarLUV

testcol <- vector()
for(i in 0:199){
  testcol <- c(testcol, hex(mixcolor(i/199,polarLUV(70,50,30),polarLUV(70,50,120))))
}

testcol2 <- vector()
for(i in 0:199){
  testcol2 <- c(testcol2, hex(mixcolor(i/199,polarLUV(70,50,300),polarLUV(70,50,210))))
}

testcols <- vector()
for(i in 0:199){
  testcols <- rbind(testcols, hex(mixcolor(i/199,hex2RGB(testcol),hex2RGB(testcol2))))
}


rescale<-function(x,newrange) {
  if(nargs() > 1 && is.numeric(x) && is.numeric(newrange)) {
    if(newrange[1] > newrange[2]) {
      newmin<-newrange[2]
      newrange[2]<-newrange[1]
      newrange[1]<-newmin
    }
    xrange<-range(x)
    if(xrange[1] == xrange[2]) stop("can't rescale a constant vector!")
    mfac<-(newrange[2]-newrange[1])/(xrange[2]-xrange[1])
    return(newrange[1]+(x-xrange[1])*mfac)
  }
  else {
    cat("Usage: rescale(x,newrange)\n")
    cat("\twhere x is a numeric object and newrange is the min and max of the new range\n")
  }
}

#' Generate diverging colors in HCL colour space.
#'
#' A function to generate colors in Hue-Chroma-Luminance colour scheme for mapping bioregions.
#'
#' @param x An object of class \code{\link[vegan]{metaMDS}}
#'
#' @rdname hexcols
#'
#' @keywords phyloregion
#' @importFrom colorspace hex
#' @export
#' @return A range of discrete colors differentiating between bioregions in
#' terms of their shared relationships.
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' library(vegan)
#' data(dune)
#' c1 <- metaMDS(dune)
#' hexcols(c1)
#' plot(c1$points, pch = 21, cex = 7, bg = hexcols(c1), las = 1)
#'
hexcols <- function(x){
  k=x$nobj
  plotcols <- rep(0, k)
  names(plotcols) <- as.character(1:k)
  rans <- c(max(x$points[,1])-min(x$points[,1]),
            max(x$points[,2])-min(x$points[,2]))

  stan <- rans/max(rans)
  reord <- cbind(rescale(x$points[,1], c(stan[1],0)),
                 rescale(x$points[,2],c(stan[2],0)))  # rescale function in r stuff
  for(i in rownames(reord))
  {
    plotcols[i] <- testcols[ceiling(reord[i,1]*199)+1,ceiling(reord[i,2]*199)+1]
  }
  plotcols
}
