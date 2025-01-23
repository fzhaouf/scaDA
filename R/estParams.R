#' estimating ZINB parameters
#'
#' @param object scaDAdataset object
#' @param group.1 cell type 1
#' @param group.2 cell type 2
#'
#' @return a scaDAdataset object contains parameter estimates in params slot
#' @export
#' @importFrom stats median plogis pchisq
#' @importFrom progress progress_bar
estParams <- function(object, group.1=NULL, group.2=NULL){
  message("start initial parameter estimates")

  count <- object@count
  peak_names <- rownames(count)  # Preserve peak names
  # normlaization factor
  sfs <- apply(count,2,mean); sfs=sfs/median(sfs)

  if(is.null(group.2)){
    cells_loc <- sample.loc(object=object, group.1 = group.1, method="balanced")
    group.1.loc <- cells_loc$group.1.loc  # assuming group 1 is case
    group.2.loc <- cells_loc$group.2.loc  # assuming group 2 is reference
  } else { # if both group.1 and group.2 are specified
    group.1.loc <- which(object@colData==group.1)  # assuming c1 is case
    group.2.loc <- which(object@colData==group.2)  # assuming c2 is reference
  }
  # subset dat
  dat <- count[,c(group.1.loc, group.2.loc)]
  npeak <- dim(dat)[1]
  nsam <- dim(dat)[2]
  nsam1 <- length(group.1.loc)
  nsam2 <- length(group.2.loc)
  poolCol <- c(1:nsam)
  cond1Col <- c(1:nsam1)
  cond2Col <- c((nsam1+1):nsam)

  counts_pooled <- dat[,poolCol]
  counts_cell1 <- dat[,cond1Col]
  counts_cell2 <- dat[,cond2Col]
  # seq depth adj factors
  sfs_pooled <- sfs[c(group.1.loc, group.2.loc)]
  sfs_cell1 <- sfs[group.1.loc]
  sfs_cell2 <- sfs[group.2.loc]
  # df store estimates
  est_params_pooled <- matrix(NA,dim(counts_pooled)[1],ncol=3)
  est_params_cell1 <- matrix(NA,dim(counts_cell1)[1],ncol=3)
  est_params_cell2 <- matrix(NA,dim(counts_cell2)[1],ncol=3)
  colnames(est_params_pooled) <- c("mu","phi","p0")
  rownames(est_params_pooled) <- peak_names
  colnames(est_params_cell1) <- c("mu","phi","p0")
  rownames(est_params_cell1) <- peak_names
  colnames(est_params_cell2) <- c("mu","phi","p0")
  rownames(est_params_cell2) <- peak_names

  pval_zinb3p <- NULL
  tstats <- NULL
  pb <- progress_bar$new(
    format = "  [:bar] :percent :elapsedfull",
    total = npeak, clear = FALSE, width = 60
  )
  suppressWarnings({
  for (i in 1:npeak){
    counts.peak <- as.numeric(counts_pooled[i, ])
    dat.whole.peak <- data.frame(counts = counts.peak)

    ctrl <- pscl::zeroinfl.control(method = "L-BFGS-B")
    ctrl$reltol <- NULL
    ctrl$factr <- 1e-3/.Machine$double.eps
    #fit zero infl ng to whole sample
    m1 <- pscl::zeroinfl(counts ~ 1+offset(log(sfs_pooled)) | 1+offset(log(sfs_pooled)),
                   data = dat.whole.peak, dist = "negbin" ,control = ctrl)
    mu <- exp(m1$coefficients$count)
    if (m1$theta>150){
      theta <- 150
    } else {
      theta <- m1$theta
    }
    p0 <- plogis(m1$coefficients$zero)
    est_params_pooled[i,] <- c(mu,theta,p0)
    logL_null <- zinb.loglink(counts=counts.peak,p=p0,u=mu,k=1/theta)
    #fit zinb to cond1
    counts.peak <- as.numeric(counts_cell1[i, ])
    dat.cell1.peak <- data.frame(counts = counts.peak)
    m2 <- pscl::zeroinfl(counts ~ 1+offset(log(sfs_cell1)) | 1+offset(log(sfs_cell1)),
                   data = dat.cell1.peak, dist = "negbin" ,control = ctrl)
    mu <- exp(m2$coefficients$count)
    if (m2$theta>150){
      theta <- 150
    } else {
      theta <- m2$theta
    }
    p0 <- plogis(m2$coefficients$zero)
    est_params_cell1[i,] <- c(mu,theta,p0)
    logL_alter_1 <- zinb.loglink(counts=counts.peak,p=p0,u=mu,k=1/theta)
    #fit zinb to cond2
    counts.peak <- as.numeric(counts_cell2[i, ])
    dat.cell2.peak <- data.frame(counts = counts.peak)
    m3 <- pscl::zeroinfl(counts ~ 1 + offset(log(sfs_cell2)) | 1+offset(log(sfs_cell2)),
                   data = dat.cell2.peak, dist = "negbin" ,control = ctrl)
    mu <- exp(m3$coefficients$count)
    if (m3$theta>150){
      theta <- 150
    } else {
      theta <- m3$theta
    }
    p0 <- plogis(m3$coefficients$zero)
    est_params_cell2[i,] <- c(mu,theta,p0)
    logL_alter_2 <- zinb.loglink(counts=counts.peak,p=p0,u=mu,k=1/theta)
    logL_alter <- logL_alter_1+logL_alter_2

    #LRT test follow x 3 degress of freedom
    test.stats <- -2*(logL_null - logL_alter)
    pvl <- pchisq(test.stats, df=3, lower.tail = FALSE)
    pval_zinb3p <- c(pval_zinb3p,pvl)
    tstats <- c(tstats, test.stats)
    pb$tick()
  }

  est_params_pooled <- data.frame(est_params_pooled)
  est_params_cell1 <- data.frame(est_params_cell1)
  est_params_cell2 <- data.frame(est_params_cell2)

  object@params <- list(g1 = group.1.loc,
                        g2 = group.2.loc,
                        param_pooled = est_params_pooled,
                        param_g1 = est_params_cell1,
                        param_g2 = est_params_cell2)

  })

  return(object)

}


