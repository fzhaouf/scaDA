sample.loc <- function(object,group.1,method=c("balanced","unbalanced")){
  # for one vs. all other case, use stratified sampling to get balanced sample size
  if(method=="balanced"){
    group.1.loc <- which(object@colData==group.1)
    num.group.1 <- length(group.1.loc)

    celltypes <- names(sort(table(object@colData), decreasing = TRUE))
    cellnums <- sort(table(object@colData), decreasing = TRUE)
    mat <- data.frame(celltype=celltypes,cellnum=as.numeric(cellnums))

    mat2 <- mat[mat$celltype!=group.1,]
    mat2$prop <- as.numeric(mat2$cellnum)/sum(as.numeric(mat2$cellnum))
    ncell <- nrow(mat2)

    group.2.loc <- NULL
    for (i in 1:ncell){
      cell <- mat2$celltype[i]
      cell.loc <- which(object@colData==cell)
      size <- floor(num.group.1 * mat2$prop[i])
      cell.loc.samp <- sample(cell.loc, size, replace = FALSE)
      group.2.loc <- c(group.2.loc, cell.loc.samp)
    }
  }

  if(method=="unbalanced"){
    group.1.loc <- which(object@colData == group.1)
    group.2.loc <- which(object@colData != group.1)
  }

  return(list(group.1.loc=group.1.loc, group.2.loc=group.2.loc))
}


#loglikelihood function for ZINB
zinb.loglink <- function(counts,p,u,k){
  counts <- as.numeric(counts)
  dens <- numeric(length(counts))
  for (i in 1:length(counts)){
    if (counts[i]==0){
      dens[i] <- p+(1-p)*(1/(1+k*u))^(1/k)
    } else {
      g <- gamma(counts[i]+1/k)/(gamma(1/k)*gamma(counts[i]+1))
      dens[i] <- (1-p)*g*(1/(1+k*u))^(1/k)*((k*u)/(1+k*u))^counts[i]
    }
    dens <- unlist(dens)
  }
  loglink <- sum(log(dens))
  return(loglink)
}

posterior.phi <- function(counts,p,u,k,m,tao){
  counts <- as.numeric(counts)
  zros <- counts[counts==0]
  nzros <- counts[counts!=0]
  n1 <- length(zros)
  n2 <- length(nzros)
  alpha <- 1/k
  temp1 <- 1/(1+k*u)
  get.phi <- n1*log(p+(1-p)*(temp1^alpha))+sum(log(gamma(nzros+alpha)))-n2*log(gamma(alpha))+
    n2*alpha*log(temp1)+sum(nzros*(log(k*u)-log(1+k*u)))-((log(k)-m)^2)/(2*(tao^2))-
    log(k)-log(tao)
  return(get.phi)
}

#### two kind of tasks: normal sample - human brain; disease sample - AD dataset

