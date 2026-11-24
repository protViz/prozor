
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


