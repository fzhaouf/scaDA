#' optimization of mu and phi parameters
#'
#' @param object scaDAdataset object
#'
#' @return scaDAdataset object
#' @return
#' @export
#' @importFrom stats optimise pchisq p.adjust
#' @importFrom progress progress_bar
#' @examples
optParams <- function(object){ ## output phi is dispersion
  message("start optimize parameter estimates")
  count = object@count
  celltype1.loc = object@params$c1  # assuming c1 is reference
  celltype2.loc = object@params$c2
  # subset dat
  dat = count[,c(celltype1.loc, celltype2.loc)]
  npeak = dim(dat)[1]
  nsam = dim(dat)[2]
  nsam1 = length(celltype1.loc)
  nsam2 = length(celltype2.loc)
  poolCol = c(1:nsam)
  cond1Col = c(1:nsam1)
  cond2Col = c((nsam1+1):nsam)

  tol = 1e-3
  nitr = 100

  total_iterations <- npeak*3
  pb <- progress::progress_bar$new(total = total_iterations, clear = FALSE)

  #opt mu and prev for cond1
  est_params_cell1 = data.frame(object@params$param_c1)
  mu_opt_c1 = NULL
  prev_opt_c1 = NULL
  for (i in 1:npeak){
    counts = dat[i,cond1Col]
    prev = est_params_cell1[i,]$p0
    nb_mu = est_params_cell1[i,]$mu
    max.mu = max(est_params_cell1$mu)
    nb_phi = est_params_cell1[i,]$phi
    pb$tick()
    for (k in 1:nitr){
      prev0 = prev
      nb_mu0=nb_mu
      nb_mu = optimise(zinb.loglink,c(0.01,max.mu),tol=0.0001,maximum = TRUE, counts=counts,p=prev0,k=nb_phi)$maximum
      prev = optimise(zinb.loglink,c(0.01,1),tol=0.0001,maximum = TRUE, counts=counts,u=nb_mu,k=nb_phi)$maximum
      if(abs(nb_mu0-nb_mu)/abs(nb_mu0)<tol){
        if(abs(prev0-prev)/abs(prev0)<tol){
          break
        }
      }
    }
    mu_opt_c1 = c(mu_opt_c1,nb_mu)
    prev_opt_c1 = c(prev_opt_c1,prev)
  }
  est_params_cell1$mu = mu_opt_c1
  est_params_cell1$p0 = prev_opt_c1

  #opt mu and prev for cond2 data
  est_params_cell2 = data.frame(object@params$param_c2)
  mu_opt_c2 = NULL
  prev_opt_c2 = NULL
  for (i in 1:npeak){
    counts = dat[i,cond2Col]
    prev = est_params_cell2[i,]$p0
    nb_mu = est_params_cell2[i,]$mu
    max.mu = max(est_params_cell2$mu)
    nb_phi = est_params_cell2[i,]$phi
    pb$tick()
    for (k in 1:nitr){
      prev0 = prev
      nb_mu0=nb_mu
      nb_mu = optimise(zinb.loglink,c(0.01,max.mu),tol=0.0001,maximum = TRUE, counts=counts,p=prev0,k=nb_phi)$maximum
      prev = optimise(zinb.loglink,c(0.01,1),tol=0.0001,maximum = TRUE, counts=counts,u=nb_mu,k=nb_phi)$maximum
      if(abs(nb_mu0-nb_mu)/abs(nb_mu0)<tol){
        if(abs(prev0-prev)/abs(prev0)<tol){
          break
        }
      }
    }
    mu_opt_c2 = c(mu_opt_c2,nb_mu)
    prev_opt_c2 = c(prev_opt_c2,prev)
  }
  est_params_cell2$mu = mu_opt_c2
  est_params_cell2$p0 = prev_opt_c2
  est_params_pooled = object@params$param_pooled

  ## testing using shrinked phi and optimized mu and prev
  pval_zinb_shrink_opt = NULL
  tstats = NULL
  for (i in 1:npeak){
    counts = dat[i,poolCol]
    prev = est_params_pooled[i,]$p0
    nb_mu = est_params_pooled[i,]$mu
    nb_phi_shrink = est_params_pooled[i,]$phi
    logL_null = zinb.loglink(counts=counts,p=prev,u=nb_mu,k=nb_phi_shrink)

    counts = dat[i,cond1Col]
    prev = est_params_cell1[i,]$p0
    nb_mu = est_params_cell1[i,]$mu
    nb_phi_shrink = est_params_cell1[i,]$phi
    logL_alter_1 = zinb.loglink(counts=counts,p=prev,u=nb_mu,k=nb_phi_shrink)

    counts = dat[i,cond2Col]
    prev = est_params_cell2[i,]$p0
    nb_mu = est_params_cell2[i,]$mu
    nb_phi_shrink = est_params_cell2[i,]$phi
    logL_alter_2 = zinb.loglink(counts=counts,p=prev,u=nb_mu,k=nb_phi_shrink)
    logL_alter = logL_alter_1+logL_alter_2

    test.stats = -2*(logL_null - logL_alter)

    pvl = pchisq(test.stats, df=3, lower.tail = FALSE)
    pval_zinb_shrink_opt = c(pval_zinb_shrink_opt,pvl)
    tstats = c(tstats,test.stats)
    pb$tick()
  }
  result = data.frame(tstats=tstats, pval=pval_zinb_shrink_opt)

  result$FDR = p.adjust(result$pval,method='fdr')
  # calculate fold change
  m1 = (1-est_params_cell1$p0)*est_params_cell1$mu
  m2 = (1-est_params_cell2$p0)*est_params_cell2$mu
  foch = m2/m1
  result$log2fc = log(foch,2)

  object@params = list(c1 = celltype1.loc,c2 = celltype2.loc, result = result,
                       param_pooled = est_params_pooled,param_c1 = est_params_cell1,param_c2 = est_params_cell2)

  return(object)
}


