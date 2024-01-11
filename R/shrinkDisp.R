#' shrinking dispersion
#'
#' @param object scaDAdataset object
#'
#' @return scaDAdataset object containing shrinked dispersion estimates
#' @export
#' @importFrom stats median IQR optimise
#' @importFrom progress progress_bar
#'
#' @examples
shrinkDisp <- function(object){ #output phi is dispersion!
  message("start shrink dispersion parameter")

  count <- object@count
  group.1.loc <- object@params$g1  # assuming g1 is case
  group.2.loc <- object@params$g2
  dat <- count[,c(group.1.loc, group.2.loc)]
  npeak <- dim(dat)[1]
  nsam <- dim(dat)[2]
  nsam1 <- length(group.1.loc)
  nsam2 <- length(group.2.loc)
  poolCol <- c(1:nsam)
  cond1Col <- c(1:nsam1)
  cond2Col <- c((nsam1+1):nsam)
  total_iterations <- npeak*3
  pb <- progress_bar$new(
    format = "  [:bar] :percent :elapsedfull",
    total = total_iterations, clear = FALSE, width = 60
  )
  suppressWarnings({
  #shrink phi for pooled data
  est_params_pooled <- data.frame(object@params$param_pooled)
  phi_shrink_pooled <- NULL
  m <- median(log(1/est_params_pooled$phi))
  max.val <- max(1/est_params_pooled$phi)
  sigma2.mar <- (IQR(log(1/est_params_pooled$phi), na.rm=TRUE) / 1.349)^2
  tao <- sqrt(max(sigma2.mar*0.5,1e-2))
  for (i in 1:npeak){
    counts <- dat[i,poolCol]
    prev <- est_params_pooled[i,]$p0
    nb_mu <- est_params_pooled[i,]$mu
    phimax <- optimise(posterior.phi,c(0.01,max.val),tol=0.0001, maximum = TRUE, counts=counts,p=prev,u=nb_mu,m=m,tao=tao)
    phimax <- phimax$maximum
    phi_shrink_pooled <- c(phi_shrink_pooled,phimax)
    pb$tick()
  }
  est_params_pooled$phi <- phi_shrink_pooled
  #shrink phi for cond1 data
  est_params_cell1 <- data.frame(object@params$param_g1)
  phi_shrink_cond1 <- NULL
  m <- median(log(1/est_params_cell1$phi))
  max.val <- max(1/est_params_cell1$phi)
  sigma2.mar <- (IQR(log(1/est_params_cell1$phi), na.rm=TRUE) / 1.349)^2
  tao <- sqrt(max(sigma2.mar*0.5,1e-2))
  for (i in 1:npeak){
    counts <- dat[i,cond1Col]
    prev <- est_params_cell1[i,]$p0
    nb_mu <- est_params_cell1[i,]$mu
    phimax <- optimise(posterior.phi,c(0.01,max.val),tol=0.00001,maximum = TRUE, counts=counts,p=prev,u=nb_mu,m=m,tao=tao)
    phimax <- phimax$maximum
    phi_shrink_cond1 <- c(phi_shrink_cond1,phimax)
    pb$tick()
  }
  est_params_cell1$phi <- phi_shrink_cond1
  ## condition 2 data
  est_params_cell2 <- data.frame(object@params$param_g2)
  phi_shrink_cond2 <- NULL
  m <- median(log(1/est_params_cell2$phi))
  max.val <- max(1/est_params_cell2$phi)
  sigma2.mar <- (IQR(log(1/est_params_cell2$phi), na.rm=TRUE) / 1.349)^2
  tao <- sqrt(max(sigma2.mar*0.5,1e-2))
  for (i in 1:npeak){
    counts <- dat[i,cond2Col]
    prev <- est_params_cell2[i,]$p0
    nb_mu <- est_params_cell2[i,]$mu
    phimax <- suppressWarnings({optimise(posterior.phi,c(0.01,max.val),tol=0.00001,maximum = TRUE, counts=counts,p=prev,u=nb_mu,m=m,tao=tao)})
    phimax <- phimax$maximum
    phi_shrink_cond2 <- c(phi_shrink_cond2,phimax)
    pb$tick()
  }
  est_params_cell2$phi <- phi_shrink_cond2
  object@params <- list(g1 = group.1.loc,
                        g2 = group.2.loc,
                        param_pooled = est_params_pooled,
                        param_g1 = est_params_cell1,
                        param_g2 = est_params_cell2)
  return(object)
  })
}
