#Extract and impose class to SEIR fit.
extract_SEIR <- function(x){
    ext <- rstan::extract(x)
    class(ext) <- append(class(ext), 'SEIR_fit', 0)
    ext
}

#Get CI for compartment
._ci_model_ <- function(obj, n, group, .f=mean, .quantile){
    .quantile <- ifelse(.quantile>1, .quantile/100, .quantile)
    time <- seq_len(dim(obj$SEIR)[3])
    # browser()
    m <- sapply(time, function(i) sapply(group, function(j) .f(obj$SEIR[,n,i,j])))
    l <- sapply(time, function(i) sapply(group, function(j) quantile(obj$SEIR[,n,i,j], min(.quantile))))
    u <- sapply(time, function(i) sapply(group, function(j) quantile(obj$SEIR[,n,i,j], max(.quantile))))
    
    structure(list(mean=m, lower=l, upper=u), class='SEIR_obj')
}

ci_SEIR <- function(obj, what = c('S', 'E', 'I', 'R', 'incidence', 'reported'), group, .f=mean, .quantile = c(2.5, 97.5)){
    n <- which(c('S', 'E', 'I', 'R', 'incidence', 'reported') == match.arg(what))
    if (missing(group)) group <- seq_len(dim(obj$SEIR)[4])
    ._ci_model_(obj, n, group, .f, .quantile)
}

ci_SEIcIscR <- function(obj, what = c('S', 'E', 'Ic', 'Isc', 'R', 'incidence', 'subclinical'), group, .f=mean, .quantile = c(2.5, 97.5)){
    n <- which(c('S', 'E', 'Ic', 'Isc', 'R', 'incidence', 'subclinical') == match.arg(what))
    if (missing(group)) group <- seq_len(dim(obj$SEIR)[4])
    ._ci_model_(obj, n, group, .f, .quantile)
}

#Plotting method for each compartment
ggplot.SEIR_obj <- function(obj, group, ci=TRUE){
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x=1:length(obj$mean[group,])))+
    ggplot2::geom_line(ggplot2::aes(y=obj$mean[group,]))
    
    if (!ci) return(p)
    p + ggplot2::geom_ribbon(ggplot2::aes(ymin=obj$lower[group,], ymax=obj$upper[group,]), alpha=.2)
}

#Plotting method for all
ggplot.SEIR_fit <- function(obj, what, group, ci=TRUE,.quantile = c(2.5, 97.5), .t_range = c(0, Inf), .f=mean){
    require(ggplot2)
    fun <- if (length(obj) == 6) ci_SEIR else ci_SEIcIscR
    tab <- lapply(what, function(w){
        this <- as.data.frame.list(fun(obj, w, group,.f=.f, .quantile=.quantile))
        this$t <- seq_len(nrow(this))
        this$compartment <- w
        this
    })
    tab <- do.call(rbind, tab)
    # browser()
    # print(head(tab))
    tab <- subset(tab, t >= min(.t_range) & t <= max(.t_range))
    p <- ggplot(data=tab, aes(x=t))+geom_line(aes(y=mean, color=compartment))+ylab('cases')

    if (!ci) return(p)
    p + geom_ribbon(aes(ymin=lower,ymax=upper, fill=compartment), alpha=.2)
}



