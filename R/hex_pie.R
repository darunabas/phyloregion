# assing colors to pie sliced based on omega matrix
pie.col<-function(cols,d.mat){ # create pie slice colors according to number of phyloregions each hex cell belong
    new.col<-unlist(mapply(rep,cols,d.mat))
    if(length(new.col)<12){
        tt<-data.frame(table(new.col))
        new.col<-sort(c(new.col,as.character(tt[which(tt[,2]==max(tt[,2])),1])))
    }
    if(length(new.col)>12){
        tt<-data.frame(table(new.col))
        t.r<-as.character(tt[which(tt[,2]==max(tt[,2])),1])[1]
        new.col <- sort(new.col[-(which(new.col==t.r)[1])])
    }
    return(new.col)
}

#' Visualize biogeographic regions using hexagonal pie charts
#'
#' @description Similar to plot_pie function, this function takes a SpatialPolygonsDataframe with an omega matrix and plots biogreographic patterns using hexagona pies for a more aesthetic appeal.
#'
#' @param poly.df an object of SpatialPolygonsDataframe with an omega matrix
#' @keywords visualization and mapping
#' @importFrom graphics polygon
#' @importFrom rgeos gCentroid
#' @return A hexagonal pies plot on a map in geographic space
#' @author Piyal Karunarathne \email{piyalkarumail@@yahoo.com}
#'
#' @examples
#' require(sp)
#' data(africa)
#' wm<-africa$polys
#' omeg<-africa$omega
#' hx <- spsample(wm,type="hexagonal",n=nrow(omeg))
#' HexPols <- HexPoints2SpatialPolygons(hx)
#' nex.hex <- as(HexPols,"SpatialPolygonsDataFrame")
#' nex.hex@@data<-data.frame(omeg[1:length(nex.hex),])
#' hex.pie(nex.hex)
#'
#' @export
hex.pie<-function(poly.df){ # uses a SpatialPolygonsDataFrame with an omega matrix
    d.mat<-round(proportions(as.matrix(poly.df@data),margin=1)*12,0)
    cols=hue(12)
    plot(poly.df,border=F)
    for(i in seq_along(poly.df)){
        poly<-poly.df[i,]
        new.col<- pie.col(cols,d.mat[i,])
        xy<- data.frame(poly@polygons[[1]]@Polygons[[1]]@coords)
        x <- xy$x[-1] + (xy$x[-nrow(xy)] - xy$x[-1]) / 2
        y <- xy$y[-1] + (xy$y[-nrow(xy)] - xy$y[-1]) / 2
        m <- data.frame(cbind(x,y))
        coords<-data.frame(suppressWarnings(mapply(rbind,xy,m)))
        mids<-gCentroid(poly)@coords
        df0<-suppressWarnings(cbind(mids,data.frame(coords[-14,])))
        for(i in 1:13){
            polygon(c(df0[i,1],df0[i,3],df0[i+1,3]),c(df0[i,2],df0[i,4],df0[i+1,4]),col=new.col[i],border=F)
        }
    }
}

