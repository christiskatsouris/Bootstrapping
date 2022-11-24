> boot::tsboot 

function (tseries, statistic, R, l = NULL, sim = "model", endcorr = TRUE, 
    n.sim = NROW(tseries), orig.t = TRUE, ran.gen = function(tser, 
        n.sim, args) tser, ran.args = NULL, norm = TRUE, ..., 
    parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 
        1L), cl = NULL) 
{
    if (missing(parallel)) 
        parallel <- getOption("boot.parallel", "no")
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    if (parallel != "no" && ncpus > 1L) {
        if (parallel == "multicore") 
            have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") 
            have_snow <- TRUE
        if (!have_mc && !have_snow) 
            ncpus <- 1L
        loadNamespace("parallel")
    }
    statistic
    tscl <- class(tseries)
    R <- floor(R)
    if (R <= 0) 
        stop("'R' must be positive")
    call <- match.call()
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    t0 <- if (orig.t) 
        statistic(tseries, ...)
    else NULL
    ts.orig <- if (!isMatrix(tseries)) 
        as.matrix(tseries)
    else tseries
    n <- nrow(ts.orig)
    if (missing(n.sim)) 
        n.sim <- n
    class(ts.orig) <- tscl
    if ((sim == "model") || (sim == "scramble")) 
        l <- NULL
    else if ((is.null(l) || (l <= 0) || (l > n))) 
        stop("invalid value of 'l'")
    fn <- if (sim == "scramble") {
        rm(ts.orig)
        function(r) statistic(scramble(tseries, norm), ...)
    }
    else if (sim == "model") {
        rm(ts.orig)
        ran.gen
        ran.args
        function(r) statistic(ran.gen(tseries, n.sim, ran.args), 
            ...)
    }
    else if (sim %in% c("fixed", "geom")) {
        if (sim == "geom") 
            endcorr <- TRUE
        i.a <- ts.array(n, n.sim, R, l, sim, endcorr)
        ran.gen
        ran.args
        function(r) {
            ends <- if (sim == "geom") 
                cbind(i.a$starts[r, ], i.a$lengths[r, ])
            else cbind(i.a$starts[r, ], i.a$lengths)
            inds <- apply(ends, 1L, make.ends, n)
            inds <- if (is.list(inds)) 
                matrix(unlist(inds)[1L:n.sim], n.sim, 1L)
            else matrix(inds, n.sim, 1L)
            statistic(ran.gen(ts.orig[inds, ], n.sim, ran.args), 
                ...)
        }
    }
    else stop("unrecognized value of 'sim'")
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
        if (have_mc) {
            parallel::mclapply(seq_len(R), fn, mc.cores = ncpus)
        }
        else if (have_snow) {
            list(...)
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", 
                  ncpus))
                if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
                  parallel::clusterSetRNGStream(cl)
                res <- parallel::parLapply(cl, seq_len(R), fn)
                parallel::stopCluster(cl)
                res
            }
            else parallel::parLapply(cl, seq_len(R), fn)
        }
    }
    else lapply(seq_len(R), fn)
    t <- matrix(, R, length(res[[1L]]))
    for (r in seq_len(R)) t[r, ] <- res[[r]]
    ts.return(t0 = t0, t = t, R = R, tseries = tseries, seed = seed, 
        stat = statistic, sim = sim, endcorr = endcorr, n.sim = n.sim, 
        l = l, ran.gen = ran.gen, ran.args = ran.args, call = call, 
        norm = norm)
}

