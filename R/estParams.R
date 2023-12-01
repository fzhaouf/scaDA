#' estimating ZINB parameters
#'
#' @param object scaDAdataset object
#' @param celltype2 cell type of interest
#'
#' @return a scaDAdataset object contains parameter estimates in params slot
#' @export
#' @importFrom stats median plogis pchisq
#' @importFrom progress progress_bar
#'
#' @examples
estParams <- function(object, celltype2){
  message("start initial parameter estiamte")

  count = object@count
  # normlaization factor
  sfs=apply(count,2,mean); sfs=sfs/median(sfs)
  cells_loc = sample.loc(object=object, celltype2 = celltype2, method="balanced")
  celltype1.loc = cells_loc$celltype1.loc  # assuming c1 is reference
  celltype2.loc = cells_loc$celltype2.loc  # assuming c2 is case

  # subset dat
  dat = count[,c(celltype1.loc, celltype2.loc)]
  npeak = dim(dat)[1]
  nsam = dim(dat)[2]
  nsam1 = length(celltype1.loc)
  nsam2 = length(celltype2.loc)
  poolCol = c(1:nsam)
  cond1Col = c(1:nsam1)
  cond2Col = c((nsam1+1):nsam)

  counts_pooled = dat[,poolCol]
  counts_cell1= dat[,cond1Col]
  counts_cell2 = dat[,cond2Col]

  sfs_pooled=sfs[c(celltype1.loc, celltype2.loc)]
  sfs_cell1=sfs[celltype1.loc]
  sfs_cell2=sfs[celltype2.loc]

  est_params_pooled = matrix(NA,dim(counts_pooled)[1],ncol=3)
  est_params_cell1 = matrix(NA,dim(counts_cell1)[1],ncol=3)
  est_params_cell2 = matrix(NA,dim(counts_cell2)[1],ncol=3)

  colnames(est_params_pooled) = c("mu","phi","p0")
  colnames(est_params_cell1) = c("mu","phi","p0")
  colnames(est_params_cell2) = c("mu","phi","p0")

  pval_zinb3p = NULL
  tstats = NULL

  pb <- progress_bar$new(
    format = "  [:bar] :percent :elapsedfull",
    total = npeak, clear = FALSE, width = 60
  )

  suppressWarnings({

  for (i in 1:npeak){

    counts.peak = as.numeric(counts_pooled[i, ])
    dat.whole.peak = data.frame(counts = counts.peak)

    ctrl=pscl::zeroinfl.control(method = "L-BFGS-B")
    ctrl$reltol=NULL
    ctrl$factr=1e-3/.Machine$double.eps

    #fit zero infl ng to whole sample
    m1 <- pscl::zeroinfl(counts ~ 1+offset(log(sfs_pooled)) | 1+offset(log(sfs_pooled)),
                   data = dat.whole.peak, dist = "negbin" ,control = ctrl)
    summary(m1)

    mu=exp(m1$coefficients$count)
    if (m1$theta>150){
      theta = 150
    } else {
      theta = m1$theta
    }
    p0=plogis(m1$coefficients$zero)
    est_params_pooled[i,]=c(mu,theta,p0)
    logL_null = zinb.loglink(counts=counts.peak,p=p0,u=mu,k=1/theta)
    # logL_null = logLik(m1)[1]

    #fit zinb to cond1
    counts.peak = as.numeric(counts_cell1[i, ])
    dat.cell1.peak = data.frame(counts = counts.peak)
    m2 <- pscl::zeroinfl(counts ~ 1+offset(log(sfs_cell1)) | 1+offset(log(sfs_cell1)),
                   data = dat.cell1.peak, dist = "negbin" ,control = ctrl)

    summary(m2)
    mu=exp(m2$coefficients$count)
    if (m2$theta>150){
      theta = 150
    } else {
      theta = m2$theta
    }
    p0=plogis(m2$coefficients$zero)
    est_params_cell1[i,]=c(mu,theta,p0)
    logL_alter_1 = zinb.loglink(counts=counts.peak,p=p0,u=mu,k=1/theta)
    # logL_alter_1 = logLik(m2)[1]

    #fit zinb to cond2
    counts.peak = as.numeric(counts_cell2[i, ])
    dat.cell2.peak = data.frame(counts = counts.peak)
    m3 <- pscl::zeroinfl(counts ~ 1 + offset(log(sfs_cell2)) | 1+offset(log(sfs_cell2)),
                   data = dat.cell2.peak, dist = "negbin" ,control = ctrl)

    summary(m3)
    mu=exp(m3$coefficients$count)
    if (m3$theta>150){
      theta = 150
    } else {
      theta = m3$theta
    }
    p0=plogis(m3$coefficients$zero)
    est_params_cell2[i,]=c(mu,theta,p0)
    logL_alter_2 = zinb.loglink(counts=counts.peak,p=p0,u=mu,k=1/theta)
    # logL_alter_2 = logLik(m3)[1]

    logL_alter=logL_alter_1+logL_alter_2

    #LRT test follow x 3 degress of freedom : null: p1=p2, mu1=mu2 phi1=phi2 vs p1!=p2, mu1!=mu2, phi1!=phi2
    test.stats = -2*(logL_null - logL_alter)
    pvl = pchisq(test.stats, df=3, lower.tail = FALSE)
    pval_zinb3p = c(pval_zinb3p,pvl)
    tstats = c(tstats, test.stats)

    pb$tick()

  }

  est_params_pooled = data.frame(est_params_pooled)
  est_params_cell1 = data.frame(est_params_cell1)
  est_params_cell2 = data.frame(est_params_cell2)

  object@params = list(c1 = celltype1.loc,c2 = celltype2.loc,
                       param_pooled = est_params_pooled,param_c1 = est_params_cell1,param_c2 = est_params_cell2)
  return(object)
  })
}




