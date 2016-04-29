##' @import data.table
##' @import magrittr
##' @importFrom parallel mclapply
NULL

globalVariables(".")

##' block.glm
##' 
##' Regression models for genomic data often assume there is
##' independence between neighbouring genomic elements when, in
##' reality, there is spatial dependence.  This function implements a
##' block bootstrap method for estimating correct variances of
##' paramter estimates.
##'
##' Note that this function uses `mclapply` to parallelise the
##' bootstrapping.  Please set `mc.cores` to something sensible, eg
##' \code{options(mc.cores=10)}
##' if you have 10 cores.
##'
##' @param f.lhs character vector, left hand side of a formula, the
##'     model(s) to be fit will be defined by
##'     `paste(f.lhs, f.rhs, sep=" ~ ")`
##' @param f.rhs character string, right hand side of a formula
##' @param data data.table containing the columns referred to in f.lhs
##'     and f.rhs
##' @param order.by if not `NULL`, the name of a column in `data` on
##'     which it should be sorted
##' @param strat.by if not `NULL`, the name of a column in `data` on
##'     which it should be stratified before block sampling.  Eg, if
##'     you are considering genomic data, you should stratify by
##'     chromosome as there should be no spatial correlation between
##'     chromosomes
##' @param block.size size of blocks of contiguous observations that
##'     will be sampled for bootstrap estimation of variance of
##'     parameter estimates
##' @param B number of bootstrap estimates
##' @param ... other arguments passed to `glm()` (eg
##'     `family="binomial"`)
##' @return data.table giving the estimated effect ("beta") of each
##'     item in f.rhs on each item in f.lhs, together with block
##'     bootstrap estimates of confidence interval (beta.025,
##'     beta.975) and standard error (se.beta) and the number of
##'     bootstraps on which those estimates are based.
##' @author Chris Wallace and Oliver Burren
##' @export
block.glm<-function(f.lhs,
                    f.rhs,
                      data,
                      order.by=NULL,
                      strat.by=NULL,
                      block.size=20,
                      B=200,
                      ...){

#  message(f)

    ## compute valid starting points for blocks
    if(!is.null(order.by))
        data<-data[order(data[[order.by]]),]
    if(!is.null(strat.by)) {
        dsplit <- split(1:nrow(data),as.character(data[[strat.by]]))
        dsplit <- dsplit[ sapply(dsplit,length)>=block.size ]
    } else {
        dsplit <- list(1:nrow(data))
    }
    dsplit.idx<- lapply(dsplit,function(i){
        ret<-i[1:(length(i)-block.size)]
        names(ret)<-NULL
        ret
    }) %>% unlist()

  ##   effect size estimates  
    f1 <- paste(f.lhs,f.rhs,sep=" ~ ")
    fun <- function(data,...) {
        tmp <- lapply(f1, function(f) {
            glm(as.formula(f), data=data) %>% coef()
        })  %>% do.call("rbind",.)
        rownames(tmp) <- f.lhs
        as.data.table(tmp)
    }
    effects <- fun(data, ...)


  ## we want to select n/block.size blocks
  sample.no<-ceiling(nrow(data)/block.size)
  ## bootstrap B
  RESULTS<-mclapply(1:B,function(z){
    cat(".")
    samp<-sample(dsplit.idx,sample.no)
    idx<-lapply(samp,function(x){
      seq(from=x,to=x+block.size,by=1)
    }) %>% unlist()
    fun(data[idx,],...)
  })

    byvars <- c("y","x")
    RESULTS %<>% lapply(., melt, varnames=byvars) %>%
    do.call("rbind",.)
    results <- as.data.table(RESULTS)
    setkeyv(results,byvars)
    results <- results[,.(beta.med=median(value),
                          beta.025=quantile(value,0.025),
                          beta.975=quantile(value,0.975),
                          se.beta=sd(value),
                          nboot=.N),
                       by=byvars]
  
  beta <- melt(effects,varnames=byvars)
  beta <- as.data.table(beta)
  setnames(beta,"value","beta")
  results <- merge(results,beta,by=byvars)
  return(results)  
}
