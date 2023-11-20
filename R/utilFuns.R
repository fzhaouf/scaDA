sample.loc = function(object,celltype2,method=c("balanced","unbalanced")){

  if(method=="balanced"){
    celltype2.use = celltype2
    celltype2.loc=which(object@colData==celltype2.use)
    num.celltype2 = length(celltype2.loc)

    celltypes=names(sort(table(object@colData), decreasing = TRUE))
    cellnums=sort(table(object@colData), decreasing = TRUE)
    mat=data.frame(celltype=celltypes,cellnum=as.numeric(cellnums))
    mat

    mat2 = mat[mat$celltype!=celltype2,]
    mat2$prop = as.numeric(mat2$cellnum)/sum(as.numeric(mat2$cellnum))
    ncell = nrow(mat2)

    celltype1.loc = NULL
    for (i in 1:ncell){
      cell = mat2$celltype[i]
      cell.loc=which(object@colData==cell)
      size = floor(num.celltype2 * mat2$prop[i])
      cell.loc.samp = sample(cell.loc, size, replace = FALSE)
      celltype1.loc = c(celltype1.loc, cell.loc.samp)
    }
  }

  if(method=="unbalanced"){
    celltype2.use = celltype2
    celltype2.loc=which(object@colData == celltype2.use)
    celltype1.loc = which(object@colData != celltype2.use)
  }

  return(list(celltype1.loc=celltype1.loc,celltype2.loc=celltype2.loc))
}


#loglikelihood function for ZINB
zinb.loglink = function(counts,p,u,k){
  counts = as.numeric(counts)
  dens = numeric(length(counts))
  for (i in 1:length(counts)){
    if (counts[i]==0){
      dens[i]=p+(1-p)*(1/(1+k*u))^(1/k)
    } else {
      g = gamma(counts[i]+1/k)/(gamma(1/k)*gamma(counts[i]+1))
      dens[i] = (1-p)*g*(1/(1+k*u))^(1/k)*((k*u)/(1+k*u))^counts[i]
    }
    dens = unlist(dens)
  }
  loglink = sum(log(dens))
  return(loglink)
}

posterior.phi = function(counts,p,u,k,m,tao){
  counts = as.numeric(counts)
  zros = counts[counts==0]
  nzros = counts[counts!=0]
  n1 = length(zros)
  n2 = length(nzros)
  alpha = 1/k
  temp1 = 1/(1+k*u)
  get.phi= n1*log(p+(1-p)*(temp1^alpha))+sum(log(gamma(nzros+alpha)))-n2*log(gamma(alpha))+
    n2*alpha*log(temp1)+sum(nzros*(log(k*u)-log(1+k*u)))-((log(k)-m)^2)/(2*(tao^2))-
    log(k)-log(tao)
  return(get.phi)
}


# penalized likelihood function for estimate hyperparameters m and tao
#loglikelihood function for lognormal distribution
penalized_likelihood_lognormal <- function(params, log_data) {
  mu = params[1]
  sigma = params[2]
  n = length(log_data)
  if(sigma <= 0) return(Inf)  # Prevent non-finite values
  ll = -n/2 * log(2 * pi) - n * log(sigma) - sum((log_data - mu)^2) / (2 * sigma^2)
  return(-ll)  # Negative because most optimizers minimize
}


