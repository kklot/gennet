# =========================================================================
# Generating network age-specific
# -------------------------------------------------------------------------
# - a target population age-structure
# - contact distribution (Mossong et al Plos Medicine)
# - contact matrix (Mossong et al Plos Medicine)
# -------------------------------------------------------------------------
N    <- 1000
seed <- 123
# -------------------------------------------------------------------------
# assign age: depends on the target population age distribution
# -------------------------------------------------------------------------
  SR             <- dget("SR")  # Sierra Leon population
  if (!is.null(seed)) set.seed(seed)
  age           <- sample(SR$age, N, TRUE, SR$prob)
  # add ages 70/80+
  age[age==70] <- sample(70:80, length(age[age==70]), replace=1)
# -------------------------------------------------------------------------
# assign #contacts: depends on both target age and POLYMOD data
# -------------------------------------------------------------------------
  agecont      <- orderBy(dget("POLYMODtab1"), 1) # POLYMOD Tab.1
  agecont$rk   <- rank(agecont$age)  # Rank age-group by contact
  POLYMODbreak <- c(0,4,9,14,19,29,39,49,59,69,80) 
  ageGrp <- cut(age, breaks=POLYMODbreak, include.lowest=1)
  ageGrp <- as.numeric(ageGrp) # mapping age to ageGrp
  agerk  <- sapply(ageGrp, function(x) agecont$rk[x==agecont$nmr])
  age    <- rev(orderBy(cbind(age, agerk), 2)[, 1]) # large to small
  # Sampling from contact distribution
  # TODO: assign all with zero contact to >18 (base on POLYMOD data)
  # -----------------------------------------------------------------------
  distCont <- dget("distCont")
  ncont    <- sort(sample(distCont$freq, N, 1, distCont$prob), TRUE)
# -------------------------------------------------------------------------
# Contact matrix (prob.) POLYMOD data: averaging all countries
# -------------------------------------------------------------------------
  M <- dget("M")  # image(M)
  Mbrk <- c(0,4,9,14,19,24,29,34,39,44,49,54,59,64,69,70)
  Mnmr <- 1:15
  Grp  <- as.numeric(cut(age, breaks=Mbrk, include.lowest=1))
  Grp[is.na(Grp)] <- 15 # 70+ age-grp

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
  # g <- set_vertex_attr(g, 'n', value=ncont)
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
system.time(g <- genNet(N, age))
system.time(g1 <- genNet(N, age))
plotNet <- function(M, vs=log(degree(M)+2)*2, vc=vertex_attr(g, "age"), ly=layout_with_fr, ...) plot(g, vertex.size=vs, vertex.label=NA, edge.arrow.mode=0, layout=ly, vertex.color=vc, edge.color="gray90",...) 
plotNet(g)
# -------------------------------------------------------------------------
# Validate
# -------------------------------------------------------------------------
require(gridExtra) 
require(lattice)
# Compare age-distribution
ks.test(age, vertex_attr(g, "age"))
# Compare contact-distribution
findDg <- function(x, graph=g) {
  # find neibours
  age <- vertex_attr(graph, "age")
  Grp <- as.numeric(cut(age, breaks=Mbrk, include.lowest=1))
  Grp[is.na(Grp)] <- 15
  n <- length(which(Grp==x))
  tmp <- unlist(sapply(which(Grp==x), function(x) degree(graph, x))) 
  return(m = mean(tmp), sd=sd(tmp)) 
}
# Compare contact-distribution
ks.test(degree(g), ncont)
# Plot contact-pattern
findAge <- function(x, graph) {
  # find neibours
  age <- vertex_attr(graph, "age")
  Grp <- as.numeric(cut(age, breaks=Mbrk, include.lowest=1))
  Grp[is.na(Grp)] <- 15 # 70+ age-grp
  n <- length(which(Grp==x))
  tmp <- unlist(sapply(which(Grp==x), function(x) neighbors(graph, x))) 
  tab1 <- as.data.frame(table(Grp[tmp])) # find their age Grp and count
  tab1[, 1] <- as.numeric(levels(tab1[,1]))
  tab1[, 2] <- tab1[, 2]/n
  noC <- which(agen %notin% tab1[,1]) # add extra zeros
  tmp <- data.frame(Var1 = noC, Freq = rep(0, length(noC)))
  # combine sort and take only frequency
  return(orderBy(rbind(tab1, tmp), 1)[, 2]) 
}
v1 <- sapply(1:15, findAge, graph=g)
M1 <- as.matrix(read.csv("./data/DE.csv", header=0))

plot1 <- levelplot(log1p(v1), colorkey=FALSE, par.settings = list(axis.line = list(col = "white")), scales=list(x=list(rot=90, cex=.5, label=ages, at=1:15), y=list(cex=.5, at=1:15, label=ages)), xlab=list(cex=.7, label="Subject's age"), ylab=list(cex=.7, label="Contacts' age"), main=list(label="Simulated Matrix", cex=.7), sub=list(font=1, label="n = 1000 (Using Sierra Leon age's distribution)", cex=.7), region = TRUE)
plot2 <- levelplot(M, colorkey=FALSE, par.settings = list(axis.line = list(col = "white")), scales=list(x=list(rot=90, cex=.5, label=ages, at=1:15), y=list(cex=.5, at=1:15, label=ages)), xlab=list(cex=.7, label="Subject's age"), ylab=list(cex=.7, label="Contacts' age"), main=list(label="Data Matrix", cex=.7), sub=list(font=1, label="n = 7290, 8 European countires", cex=.7), region = TRUE)
grid.arrange(plot1,plot2, ncol=2)
savePDF("matrixCompare")
savePs("matrixCompare")
