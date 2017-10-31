genNet <- function(N=N, age=age, n.try=100) {
  require(igraph) # TODO: add check points for invalid inputs
  if (N < 100) warning("Considering increase N!")
  # Utils functions -------------------------------------------------------
  `%notin%` <- function(x, y) is.na(match(x, y, nomatch=NA_integer_))
  # pick age to make contacts according the prob. from contact matrix
  pickAge <- function(x) sample(Mnmr, ncont[x], 1, M[, Grp[x]])
  # pick id that is free and in respected age groups
  pickId <- function(x) {
    pool <- freeI[freeI %in% which(Grp==apick[x])]
    if (length(pool)==1) return(pool)
      else return(sample(pool, npick[x]))
  }
  makeLinks <- function(i, j) c(sapply(j, function(x) c(i, x)))
  # templates and storages ------------------------------------------------
  g <- make_empty_graph(n = N, directed = FALSE)
  freeI <- which(ncont > 0) # some might have no contact
  # Main loop -------------------------------------------------------------
  for (i in 1:N) {
    cat('\r', "Node:", i); flush.console()
    if (ncont[i] == 0) next # am I free?
    freeI <- freeI[freeI %notin% i] # exclude myself: no loop back
    # Check if enough peoples to contact ----------------------------------
    freeA <- rle(sort(Grp[freeI])) #freenodes by age-grp
    afree <- freeA$values  # age-grp that have freenodes
    nfree <- freeA$lengths  # of freenodes by ages
    isEnough <- noWay <- FALSE
    i.try    <- 0  # trying to find enough matched age partners
    while(!isEnough) {
      i.try <- i.try + 1
      if ( i.try > n.try ) isEnough <- noWay <- TRUE # force stop
      conti <- rle(sort(pickAge(i)))  # pick age of my contacts
      npick <- conti$lengths  # of contact by ages
      apick <- conti$values  # age-grp to pick
      if (!all(apick %in% afree)) next  # if no free age, try again
      ok <- sapply(seq_along(npick), function(x) npick[x] <= nfree[which(afree==apick[x])]) # if there are, check number enough?
      if (all(ok)) isEnough <- TRUE
    }
    if (noWay) next  # leave this guy alone
    pickedId <- unlist(sapply(seq_along(npick), function(x) pickId(x)))
    if (any(which(ncont==0) %in% pickedId) ) break
    g <- add_edges(g, makeLinks(i, pickedId) ) # add weights here if needed
    ncont[i] <- 0  # update me, no more friend
    ncont[pickedId] <- ncont[pickedId] - 1  # update friends
    freeI <- which(ncont > 0) # update free id
  }
  g <- set_vertex_attr(g, 'age', value=age)
  return(g)
}
genNet <- compiler::cmpfun(genNet)
plotNet <- function(g, vs=log(degree(g)+2)*2, vc=vertex_attr(g, "age"), ly=layout_with_fr, ...) plot(g, vertex.size=vs, vertex.label=NA, edge.arrow.mode=0, layout=ly, vertex.color=vc, edge.color="gray90",...) 
# calculate ro from sims outputs
# -------------------------------------------------------------------------
R0net <- function(out = sim, net=g, All=FALSE) {
    if (class(out$SIR) != "smidSIR") stop("Incompatible")
    n  <- igraph::delete.vertices(net, which(out$aux == 0))
    n  <- igraph::as.directed(n, mode = c("arbitrary"))
    dg <- igraph::degree(n,  mode='out')
    if (All==TRUE) r0 <- dg
      else r0 <- mean(dg)
    return(r0)
}
