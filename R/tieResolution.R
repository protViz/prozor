
inter_x <- function(x, y){
    i1 <- which(x > 0)
    i2 <- which(y > 0)
    int <- intersect(i1, i2)
    return(sum(x[int]))
}

inter_y <- function(x, y){
    i1 <- which(x > 0)
    i2 <- which(y > 0)
    int <- intersect(i1,i2)
    return(sum(y[int]))
}

if (FALSE) {
    debug(prozor:::sdiff_x)
    prozor:::sdiff_x(c(1,1,0,0,1),c(0,0,1,1,1))
    debug(prozor:::sdiff_y)
    prozor:::sdiff_y(c(1,1,0,0,1),c(0,0,1,1,1))
}

sdiff_x <- function(x, y){
    i1 <- which(x > 0)
    i2 <- which(y > 0)
    si1 <- setdiff(i1,i2)
    return(sum(x[si1]))
}

sdiff_y <- function(x, y){
    i1 <- which(x > 0)
    i2 <- which(y > 0)
    si2 <- setdiff(i2,i1)
    return(sum(y[si2]))
}


euc <- function(x, y) {
    d <- sqrt( sum( (x - y) ^ 2))
    return(d)
}


mdist <- function(x, .func = euc, diag = FALSE, upper = FALSE){
    res <- matrix(nrow = nrow(x), ncol = nrow(x))
    for (j in seq_len(nrow(x)) ) {
        if (diag) {
            for (i in j:nrow(x)) {
                res[i, j] <- .func(x[i,],x[j,])
            }
        } else {
            if (j < nrow(x)) {
                for (i in ((j + 1):nrow(x))) {
                    res[i, j] <- .func(x[i,],x[j,])
                }
            }
        }
    }
    return(if (upper) { t(res) } else { res })
}


if (FALSE) {
    dd <- data.frame(x = 1:10, y = 10:1, z = 2:11, w = rep(5,10))
    xx <- mdist(dd, diag = FALSE, upper = FALSE)
    xx <- as.dist(xx)
    cS <- mean( colSums(mm) )

    ii_x <- mdist(t( as.matrix(mm) ), .func = prozor:::inter_x)
    ii_y <- mdist(t( as.matrix(mm) ), .func = prozor:::inter_y)
    all(na.omit(ii_x == ii_y))
    ib <- which(ii_x == cS, arr.ind = TRUE)
    if(nrow(ib) > 1) {
         mmerg <- mm[,unique(as.integer(ib))]
         paste(colnames((mmerg)), collapse= ";")
         md <- mm[,-unique(as.integer(ib))]
    } else {

    }

    xx_x <- mdist(t(as.matrix(md)), .func = prozor:::sdiff_x)
    which(xx_x == cS, arr.ind = TRUE)
    xx_y <- mdist(t(as.matrix(mm)), .func = prozor:::sdiff_y)
}


