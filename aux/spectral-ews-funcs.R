# spectral-ews-funcs.R
#  Functions for calculating spectral EWS. From Miller et al. 2017.
#
# Author: Paige Miller

wavelet_filter <- function (y1, x1 = seq_along(y1), h = 55, l = 115){ #
  t1 <- data.frame(time=x1, incid=y1)
  
  ## Continuous bias corrected wavelet transform (liu et al 2007)
  wt.t1 <- biwavelet::wt(t1, do.sig=FALSE)
  h=max(which(wt.t1$scale<h))
  l=min(which(wt.t1$scale>l))
  if (!all(is.finite(c(h, l)))) {
    warning(paste("Selected wavelet bands are outside of range of computed scales.",
                  "Power is always zero by definition."))
    return(rep(0, nrow(t1)))
  }
  Power=wt.t1$power.corr[l:h, , drop = FALSE]
  stat <- apply(X=Power, FUN=sum, MARGIN=2)
  return(stat)
}

wavelet_median <- function (y1, x1 = seq_along(y1)){
  t1 <- data.frame(time=x1, incid=y1)
  
  ## Continuous bias corrected wavelet transform (liu et al 2007)
  wt.t1 <- biwavelet::wt(t1, do.sig=FALSE)
  
  Power=wt.t1$power.corr
  probfun <- function(x) cumsum(x) / sum(x)
  ppower <- apply(Power, 2, probfun)
  islower50 <- ppower < 0.5
  get_med_ind <- function(x) match(FALSE, x)
  inds <- apply(islower50, 2, get_med_ind)
  median <- wt.t1$scale[inds]
  median
}