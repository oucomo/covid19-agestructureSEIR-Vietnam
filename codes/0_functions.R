#Extract and impose class to SEIR fit.
extract_SEIR <- function(x){
    ext <- rstan::extract(x)
    class(ext) <- append(class(ext), 'SEIR_fit', 0)
    ext
}

#Get CI for compartment
ci_SEIR <- function(obj, what = c('S', 'E', 'I', 'R', 'incidence', 'reported'), group, .quantile = c(2.5, 97.5)){
    n <- which(c('S', 'E', 'I', 'R', 'incidence', 'reported') == match.arg(what))
    .quantile <- ifelse(.quantile>1, .quantile/100, .quantile)
    if (missing(group)) group <- seq_len(dim(obj$SEIR)[4])
    time <- seq_len(dim(obj$SEIR)[3])
    browser()
    m <- sapply(time, function(i) sapply(group, function(j) mean(obj$SEIR[,n,i,j])))
    l <- sapply(time, function(i) sapply(group, function(j) quantile(obj$SEIR[,n,i,j], min(.quantile))))
    u <- sapply(time, function(i) sapply(group, function(j) quantile(obj$SEIR[,n,i,j], max(.quantile))))
        
    structure(list(mean=m, lower=l, upper=u), class='SEIR_obj')
}

#Plotting method for each compartment
ggplot.SEIR_obj <- function(obj, group, ci=TRUE){
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x=1:length(obj$mean[group,])))+
    ggplot2::geom_line(ggplot2::aes(y=obj$mean[group,]))
    
    if (!ci) return(p)
    p + ggplot2::geom_ribbon(ggplot2::aes(ymin=obj$lower[group,], ymax=obj$upper[group,]), alpha=.5)
}

#Plotting method for all
ggplot.SEIR_fit <- function(obj, what, group, ci=TRUE){
    require(ggplot2)
    tab <- lapply(what, function(w){
        this <- as.data.frame.list(ci_SEIR(obj, w, group))
        this$t <- seq_len(nrow(this))
        this$compartment <- w
        this
    })
    tab <- do.call(rbind, tab)
    browser()
    print(head(tab))

    p <- ggplot(data=tab, aes(x=t, color=compartment, fill=compartment))+geom_line(aes(y=mean))

    if (!ci) return(p)
    p + geom_ribbon(aes(ymin=lower,ymax=upper))
}



