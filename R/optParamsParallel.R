#' optimization of mu and phi parameters
#'
#' @param object scaDAdataset object
#'
#' @return scaDAdataset object
#' @export
#' @import parallel
#' @importFrom stats optimise pchisq p.adjust
#' @importFrom progress progress_bar
optParamsParallel <- function(object) {
  message("start optimize parameter estimates")

  count <- object@count
  group.1.loc <- object@params$g1
  group.2.loc <- object@params$g2

  dat <- count[,c(group.1.loc, group.2.loc)]
  npeak <- dim(dat)[1]
  nsam <- dim(dat)[2]
  nsam1 <- length(group.1.loc)
  nsam2 <- length(group.2.loc)
  poolCol <- c(1:nsam)
  cond1Col <- c(1:nsam1)
  cond2Col <- c((nsam1+1):nsam)

  # Initial parameter estimates
  est_params_cell1 <- data.frame(object@params$param_g1)
  est_params_cell2 <- data.frame(object@params$param_g2)

  tol <- 1e-3
  nitr <- 10

  # Setting up parallel backend
  no_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(no_cores)
  parallel::clusterExport(cl, c("dat", "cond1Col", "cond2Col", "nitr", "tol", "est_params_cell1", "est_params_cell2", "zinb.loglink"),envir = environment())
  # Parallelize the loop for condition 1
  # Parallelize optimization for condition 1
  results_c1 <- parallel::parLapply(cl, 1:npeak, function(i) {
    counts <- dat[i, cond1Col]
    prev <- est_params_cell1[i,]$p0
    nb_mu <- est_params_cell1[i,]$mu
    max.mu <- max(est_params_cell1$mu)
    nb_phi <- est_params_cell1[i,]$phi
    for (k in 1:nitr) {
      prev0 <- prev
      nb_mu0 <- nb_mu
      nb_mu <- optimise(zinb.loglink, c(0.01, max.mu), tol = 1e-4, maximum = TRUE, counts = counts, p = prev0, k = nb_phi)$maximum
      prev <- optimise(zinb.loglink, c(0.01, 1), tol = 1e-4, maximum = TRUE, counts = counts, u = nb_mu, k = nb_phi)$maximum
      if (abs(nb_mu0 - nb_mu) / abs(nb_mu0) < tol && abs(prev0 - prev) / abs(prev0) < tol) {
        break
      }
    }
    return(list(mu = nb_mu, prev = prev))
  })
  # Process results for condition 1
  mu_opt_c1 <- sapply(results_c1, function(x) x$mu)
  prev_opt_c1 <- sapply(results_c1, function(x) x$prev)

  # Parallelize optimization for condition 2
  results_c2 <- parallel::parLapply(cl, 1:npeak, function(i) {
    counts <- dat[i, cond2Col]
    prev <- est_params_cell2[i,]$p0
    nb_mu <- est_params_cell2[i,]$mu
    max.mu <- max(est_params_cell2$mu)
    nb_phi <- est_params_cell2[i,]$phi
    for (k in 1:nitr) {
      prev0 <- prev
      nb_mu0 <- nb_mu
      nb_mu <- optimise(zinb.loglink, c(0.01, max.mu), tol = 1e-4, maximum = TRUE, counts = counts, p = prev0, k = nb_phi)$maximum
      prev <- optimise(zinb.loglink, c(0.01, 1), tol = 1e-4, maximum = TRUE, counts = counts, u = nb_mu, k = nb_phi)$maximum
      if (abs(nb_mu0 - nb_mu) / abs(nb_mu0) < tol && abs(prev0 - prev) / abs(prev0) < tol) {
        break
      }
    }
    return(list(mu = nb_mu, prev = prev))
  })

  # Process results for condition 2
  mu_opt_c2 <- sapply(results_c2, function(x) x$mu)
  prev_opt_c2 <- sapply(results_c2, function(x) x$prev)

  parallel::stopCluster(cl)

  # updates parameter estimates
  est_params_cell1$mu <- mu_opt_c1
  est_params_cell1$p0 <- prev_opt_c1
  est_params_cell2$mu <- mu_opt_c2
  est_params_cell2$p0 <- prev_opt_c2
  est_params_pooled <- object@params$param_pooled

  ## testing using shrinked phi and optimized mu and prev
  pval_zinb_shrink_opt <- NULL
  tstats <- NULL
  for (i in 1:npeak){
    counts <- dat[i,poolCol]
    prev <- est_params_pooled[i,]$p0
    nb_mu <- est_params_pooled[i,]$mu
    nb_phi_shrink <- est_params_pooled[i,]$phi
    logL_null <- zinb.loglink(counts=counts,p=prev,u=nb_mu,k=nb_phi_shrink)

    counts <- dat[i,cond1Col]
    prev <- est_params_cell1[i,]$p0
    nb_mu <- est_params_cell1[i,]$mu
    nb_phi_shrink <- est_params_cell1[i,]$phi
    logL_alter_1 <- zinb.loglink(counts=counts,p=prev,u=nb_mu,k=nb_phi_shrink)

    counts <- dat[i,cond2Col]
    prev <- est_params_cell2[i,]$p0
    nb_mu <- est_params_cell2[i,]$mu
    nb_phi_shrink <- est_params_cell2[i,]$phi
    logL_alter_2 <- zinb.loglink(counts=counts,p=prev,u=nb_mu,k=nb_phi_shrink)
    logL_alter <- logL_alter_1+logL_alter_2

    test.stats <- -2*(logL_null - logL_alter)

    pvl <- pchisq(test.stats, df=3, lower.tail = FALSE)
    pval_zinb_shrink_opt <- c(pval_zinb_shrink_opt,pvl)
    tstats <- c(tstats,test.stats)
  }
  result <- data.frame(tstats=tstats, pval=pval_zinb_shrink_opt)

  result$FDR <- p.adjust(result$pval,method='fdr')
  # calculate fold change
  m1 <- (1-est_params_cell1$p0)*est_params_cell1$mu
  m2 <- (1-est_params_cell2$p0)*est_params_cell2$mu
  foch <- m2/m1
  result$log2fc <- log(foch,2)
  # update params estimates
  object@params <- list(g1 = group.1.loc,
                        g2 = group.2.loc,
                        param_g1 = est_params_cell1,
                        param_g2 = est_params_cell2)
  # save DA test results to object@result slot
  object@result <- result
  return(object)
}
