data(africa)
tree <- africa$phylo
x <- africa$comm

subphy <- match_phylo_comm(tree, x)$phy
submat <- match_phylo_comm(tree, x)$com

pbc <- phylobeta(submat, subphy)
y <- phyloregion(pbc[[1]], shp=africa$polys)

plot_NMDS(y, cex=6)
text_NMDS(y, cex=2)
plot(y, cex=1, palette="NMDS")
plot(y, cex=1)

y_sf <- sf::st_as_sf(y$shp)
ggplot() +
  geom_sf(mapping = aes(geometry = geometry, fill = cluster), data = y_sf) +
  scale_fill_manual(data = y_sf, aes(color = COLOURS))

y$shp


x <- pbc[[1]]
k = 10
method = "average"
shp = africa$polys
function (x, k = 10, method = "average", shp = NULL, grid.membership = TRUE, ...) 
{
  Q <- as.dist(x)
  P1 <- hclust(Q, method = method)
  g <- cutree(P1, k)
  dx <- data.frame(grids = names(g), cluster = unname(g))
  x <- as.matrix(x)
  if(grid.membership == TRUE){
    n.groups <- k
    comm.groups <- dx[, 2]
    names(comm.groups) <- names(g)
    Gs <- lapply(1:n.groups, function(x){
      which(comm.groups == x)
    })
    names(Gs) <- paste("G", 1:n.groups, sep = "")
    dist.matrix <- Q
    dist.matrix <- as.matrix(dist.matrix)
    PGall <- lapply(Gs, function(x){
      dist.matrix[x, x]
    })
    
    
    afilliation_by_grp <- lapply(PGall, function(x){
      afilliation_by_grp <- matrix(NA, nrow(x), 2, dimnames = list(rownames(x), c("membership", "group")))
      for (z in 1:nrow(x)) {
        dis <- as.data.frame(x[z,])[-z,]
        afilliation_by_grp[z, 1] <- mean(dis)
      }
      return(afilliation_by_grp)
    })
    list_afilliation_by_grp <- vector(mode = "list", length = length(afilliation_by_grp))
    for(l in 1:length(afilliation_by_grp)){
      afilliation_by_grp_pad <- scales::rescale(afilliation_by_grp[[l]], c(0, 1))
      afilliation_by_grp_pad[, 2] <- l
      list_afilliation_by_grp[[l]] <- afilliation_by_grp_pad
    }
    matrix_afilliation <- do.call(rbind, list_afilliation_by_grp)
    df_membership_grid <- data.frame(membership = matrix_afilliation[, 1], grids = rownames(matrix_afilliation))
  }
  colnames(x) <- rownames(x)
  region.mat <- matrix(NA, k, k, dimnames = list(1:k, 1:k))
  for (i in 1:k) {
    for (j in 1:k) {
      region.mat[i, j] <- mean(x[names(g)[g == i], names(g)[g == 
                                                              j]])
    }
  }
  region.dist <- as.dist(region.mat)
  region.mat <- as.matrix(region.dist)
  evol_distinct <- colSums(region.mat)/(nrow(region.mat) - 
                                          1)
  evol_distinct <- data.frame(ED = evol_distinct)
  evol_distinct <- cbind(cluster = rownames(evol_distinct), 
                         data.frame(evol_distinct, row.names = NULL))
  if (length(shp) == 0) {
    r <- list(membership = dx, k = k, evol_distinct = evol_distinct, 
              region.dist = region.dist, region.df = dx)
    class(r) <- c("phyloregion")
    r
  }
  else {
    m <- sp::merge(shp, cbind(dx, membership.grid = df_membership_grid[, -2]), by = "grids")
    if (!inherits(m, "SpatialPolygons")) {
      stop("Invalid geometry, may only be applied to polygons")
    }
    m <- m[!is.na(m@data$cluster), ]
    region <- raster::aggregate(m, by = "cluster")
    m1 <- sp::merge(region, evol_distinct, by = "cluster")
    proj4string(m1) <- sp::proj4string(shp)
    c1 <- vegan::metaMDS(region.dist, trace = 0)
    c2 <- vegan::metaMDS(Q, trace = 0)
    
    v <- data.frame(colorspace::hex2RGB(hexcols(c1))@coords)
    v$r <- v$R * 255
    v$g <- v$G * 255
    v$b <- v$B * 255
    
    v2 <- data.frame(colorspace::hex2RGB(rep(hexcols(c1), table(dx$cluster)))@coords)
    v2$r <- v2$R * 255
    v2$g <- v2$G * 255
    v2$b <- v2$B * 255
    v$COLOURS <- rgb(v$r, v$g, v$b, maxColorValue = 255)
    v2$COLOURS <- rgb(v2$r, v2$g, v2$b, maxColorValue = 255)
    v2$cluster <- rownames(v2)
    y <- Reduce(function(x, y) merge(x, y, by = "cluster", 
                                     all = TRUE), list(region, v, m1))
    index <- match(dx$cluster, y$cluster)
    z <- cbind(dx, ED = y$ED[index], COLOURS = y$COLOURS[index], )
    r <- list(membership = dx, k = k, shp = y, region.dist = region.dist, 
              region.df = z, NMDS = c1)
    class(r) <- "phyloregion"
    r
  }
}

