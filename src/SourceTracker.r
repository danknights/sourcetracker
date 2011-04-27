# File: SourceTracker.r
# Author: Dan Knights
# Contact: danknights@gmail.com
# License: GPL
# Copyright: Copyright 2011, Dan Knights


# Function "sourcetracker"
#
# Gets prior counts of taxa in source environments for use by SourceTracker.
#
# Params:
# train - a sample x feature observation count matrix
# envs - a factor or character vector indicating the sample environments
# maxdepth - if not NULL, all samples with > maxdepth sequences are rarified
#
# Value: 
# Returns an object of class "sourcetracker". This is a list containing:
# sources - matrix of the total counts of each taxon in each environment
# train - a copy of the input training data
# envs - a copy of the input environment vector
"sourcetracker" <- function(train, envs, maxdepth=NULL){
    envs <- factor(envs)
    train.envs <- levels(envs)
    
    # rarefy samples above maxdepth if requested
    if(!is.null(maxdepth))    train <- rarefy(train, maxdepth)
    
    # get source environment counts, smoothed by alpha * mean source env sum
    # sources is nenvs X ntaxa
    sources <- t(sapply(split(data.frame(train), envs), colSums))
    envsums <- rowSums(sources)
    sources <- rbind(sources, rep(0,ncol(train)))
    #sources <- sweep(sources, 1, alpha1 * envsums, '+')
    #sources <- rbind(sources, rep(alpha2*mean(envsums),ncol(train)))
    rownames(sources) <- c(train.envs,"Unknown")
    colnames(sources) <- colnames(train)
    
    ret <- list(sources=sources, train=train, envs=envs)
    class(ret) <- "sourcetracker"
    return(invisible(ret))
}


# Function "predict.sourcetracker"
#
# S3-level function to estimate source proportions using a sourcetracker object.
# Returns an object of class "sourcetracker.fit"
#
# Params:
# stobj - output from function "sourcetracker"
# test - test data, a matrix of sample x feature counts
#     if test is NULL, performs leave-one-out predictions of the training
#   samples
# burnin - number of "burn-in" passes for Gibbs sampling
# nrestarts - number of times to restart the Gibbs sampling process
# n.draws.per.restart - number of Gibbs draws to collect per restart
# delay - number of passes between draws (ignored if n.draws.per.restart is 1)
# alpha1 - prior counts of each species in the training environments,
#     relative to the total number of sequences per training environment.
#   Higher values decrease the trust in the training data, and make the 
#   source environment distributions over taxa smoother
# alpha2 - prior counts of each species in the Unknown environment,
#   relative to the number of sequences in a given test sample.
#   Higher values make the Unknown environment smoother and less prone to
#   over-fitting a given training sample. 
# beta - prior counts of test sequences in each environment,
#   relative to the number of sequences in a given test sample.
#   Higher values cause a smoother distribution over source environments.
# maxdepth - if not NULL, all test samples with > maxdepth sequences are rarified
# verbosity - if > 0, print progress updates while running
#
# Value: A list containing:
# draws - an array of dimension (ndraws X nenvironments X nsamples),
#     containing all draws from gibbs sampling
# proportions - the mean proportions over all Gibbs draws for each sample
# proportions - standard deviation of the mean proportions for each sample
# train.envs - the names of the source environments
# samplenames - the names of the test samples
"predict.sourcetracker" <- function(stobj, test=NULL, 
            burnin=25, nrestarts=10, ndraws.per.restart=1, delay=10,
            alpha1=0.0001, alpha2=0.1, beta=0.0001, maxdepth=NULL,
            verbosity=1){

    if(!is.null(test)){
        sources <- stobj$sources
        T <- ncol(sources) # number of taxa
        V <- nrow(sources) # number of source envs
        N <- nrow(test) # number of sink samples
        
        
        draws <- run.gibbs(sources, test, V, T, N,
                burnin=burnin, nrestarts=nrestarts, 
                ndraws.per.restart=ndraws.per.restart, delay=delay,
                alpha1=alpha1, alpha2=alpha2, beta=beta, maxdepth=maxdepth,
                verbosity=verbosity)
    } else {  # leave-one-out 
    
        train <- stobj$train
        envs <- stobj$envs
        
        T <- ncol(train) # number of taxa
        V <- nrow(stobj$sources) # number of source envs
        N <- nrow(train) # number of sink samples
        ndraws <- nrestarts * ndraws.per.restart # total number of draws
        draws <- array(dim=c(ndraws, V, N))
        for(i in (1:N)){
            stobj.i <- sourcetracker(train[-i,], envs[-i], maxdepth=maxdepth)
            sources <- stobj.i$sources
            if(verbosity >= 1) cat(sprintf('%3d: ',i))
            draws.i <- run.gibbs(sources, train[i,], V, T, 1,
                burnin=burnin, nrestarts=nrestarts, 
                ndraws.per.restart=ndraws.per.restart, delay=delay,
                alpha1=alpha1, alpha2=alpha2, beta=beta, maxdepth=maxdepth,
                verbosity=verbosity)
            draws[,,i] <- drop(draws.i)
        }
    }

    proportions <- matrix(nrow=N, ncol=V)
    proportions_sd <- matrix(nrow=N, ncol=V)
    for(i in 1:N){
        proportions[i,] <- apply(matrix(draws[,,i], ncol=V),2,mean)
        proportions_sd[i,] <- apply(matrix(draws[,,i], ncol=V),2,sd)
    }
    
    res <- list(draws=draws, proportions=proportions,
                proportions_sd=proportions_sd,
                train.envs=rownames(sources), samplenames=rownames(test))
    class(res) <- "sourcetracker.fit"
    return(invisible(res))
}


# Function "plot.sourcetracker.fit"
#
# S3-level function to plot the SourceTracker output.
#
# Params:
# stresult - output from function "predict.sourcetracker"
# labels - Labels for samples; if NULL, uses rownames of data table
# type - One of 'pie', 'bar', or 'dist' (distribution plots). Default is 'pie'.
# gridsize - number of samples to plot per row; if NULL, will be estimated
# ... - Additional graphical parameters
"plot.sourcetracker.fit" <- function(stresult, labels=NULL, 
        type=c('pie','bar','dist')[1], gridsize=NULL, env.colors=NULL, 
        titlesize=NULL, indices=NULL, include.legend=FALSE, ...){
    if(is.null(env.colors)){
        #env.colors <- c('blue',rgb(0,128/255,0),'red',rgb(0,191/255,191/255), rgb(191/255,0,191/255))
        env.colors <- std.env.colors
    }
            
    if(is.null(indices)) indices <- 1:dim(stresult$draws)[3]
    N <- length(indices)
    V <- dim(stresult$draws)[2]

    if(include.legend) N <- N + 1

    if(!is.null(gridsize) && gridsize**2 < N)
        stop(sprintf('Please choose a gridsize of at least %d.',ceiling(sqrt(N))))
        
    if(is.null(labels)) labels <- stresult[['samplenames']]
    if(is.null(gridsize)) gridsize <- ceiling(sqrt(N))
    if(is.null(titlesize)){
        if(gridsize > 1){
            titlesize <- .8/log10(gridsize)
        } else {
            titlesize=1
        }
    } 

    ngridrows <- ceiling(N / gridsize)
    par(mfrow=c(ngridrows,gridsize))
    par(oma=c(1,1,1,1), mar=c(0,0,titlesize,0))

    # legend will occupy one full plot in the upper left
    if(include.legend){
        plot(0,0,xlim=c(0,1), ylim=c(0,1), type='n', axes=FALSE)
        leg.cex <- 0.1
        leg <- legend('topleft',stresult$train.envs, fill=env.colors, bg='white', cex=leg.cex, plot=FALSE)
        maxdim <- max(leg$rect$w, leg$rect$h)

        # resize legend to be just big enough to fill the plot (80%), or to have maximum text size 2
        while(maxdim <.8 && leg.cex <= 2){
            leg.cex <- leg.cex + 0.01
            leg <- legend('topleft',stresult$train.envs, fill=env.colors, bg='white', cex=leg.cex, plot=FALSE)
            maxdim <- max(leg$rect$w, leg$rect$h)
        }
        leg <- legend('topleft',stresult$train.envs, fill=env.colors, bg='white', cex=leg.cex)
    }

    if(type=='pie') plot.sourcetracker.pie(stresult, labels, gridsize, env.colors, titlesize, indices=indices, ...)
    if(type=='bar') plot.sourcetracker.bar(stresult, labels, gridsize, env.colors, titlesize, indices=indices, ...)
    if(type=='dist') plot.sourcetracker.dist(stresult, labels, gridsize, env.colors, titlesize, indices=indices, ...)
}


######### Internal function below ####################


# Internal SourceTracker function to run Gibbs sampling
"run.gibbs" <- function(sources, test, V, T, N,
        burnin=25, nrestarts=100, ndraws.per.restart=1, delay=10,
        alpha1=0.1, alpha2=0.1, beta=0.0001, maxdepth=NULL,
        verbosity=1){

    train.envs <- rownames(sources)
    ndraws <- nrestarts * ndraws.per.restart # total number of draws
    npasses <- burnin + (ndraws.per.restart-1) * delay + 1 # passes per restart

    # draws will hold all draws (ndraws x V x N)
    draws <- array(dim=c(ndraws, V, N))
    
    # rarefy samples above maxdepth if requested
    if(!is.null(maxdepth))    test <- rarefy(test, maxdepth)

    # sink samples must have integer counts
    if(is.null(dim(test))) test <- matrix(test,ncol=T)
    test <- round(test)
    
    # store original prior counts for "Unknown"
    unknown.prior <- sources[V,]
    
    # for each sink sample
    for(i in 1:N){
        sink <- test[i,]
        D <- sum(sink) # sink sample depth
        
        # precalculate denominator for Pr(env in sample)
        p_v_denominator = max(1,(D-1) + V*beta*D)
        
        # get taxon index for each sequence
        tax.cumsum <- cumsum(sink)
        tax.ix <- sapply(1:D,function(x) min(which(x<=tax.cumsum)))

        drawcount <- 1 # keeps running count of draws for this sample
        # for each restart
        for(j in 1:nrestarts){
            
            z <- sample(V,D,replace=TRUE) # random env assignments
            sources[V,] <- unknown.prior # prior counts of taxa in Unknown
            sources[-V,] <- sources[-V,] + alpha1 * D # add relative alpha prior counts
            sources[V,] <- sources[V,] + alpha2 * D # add relative alpha prior counts
            
            # tally counts in envs
            envcounts <- rep(beta * D, V)
            for(ix in 1:D){
                if(z[ix] == V)    sources[V,tax.ix[ix]] <- sources[V,tax.ix[ix]] + 1
                envcounts[z[ix]] <- envcounts[z[ix]] + 1
            }

            for(rep in 1:npasses){
                rand.ix <- sample(D) # random order for traversing sequence

                cnt <- 0
                for(ix in rand.ix){
                    taxon <- tax.ix[ix]
                     # remove this sequence from all counts
                    envcounts[z[ix]] <- envcounts[z[ix]] - 1
                    if(z[ix] == V)    sources[V,taxon] <- sources[V,taxon] - 1


                    # get relative PDF over env assignments
                    p_tv <- sources[,taxon] / rowSums(sources) # Pr(taxon | env)
                    p_v <- envcounts/p_v_denominator# Pr(env in sample)
                    
                    
                    # re-sample this sequence's env assignment
                    z[ix] <- sample(1:V, prob=p_tv * p_v, size=1)

                    if(1==0){
                        cat(sprintf('otu index %d:\n', taxon))
                        cat(sprintf('%.5f\t%.5f\n',p_tv, p_v))
                        cat('Actual counts: ')
                        cat(sprintf('%.5f',sources[,taxon]),sep='\t')
                        
                        cat('\n\n')
                        if(cnt > 10    ) stop()
                        cnt <- cnt + 1
                    }
                    # replace this sequence in all counts
                    envcounts[z[ix]] <- envcounts[z[ix]] + 1
                    if(z[ix] == V)    sources[V,taxon] <- sources[V,taxon] + 1
                    

                }
                # take sample
                if(rep > burnin && (((rep-burnin) %% delay)==1 || delay<=1)){
                        
                    # save current mixing proportions
                    draws[drawcount,,i] <- round((envcounts - beta*D) / D,7)
                    draws[drawcount,,i] <- draws[drawcount,,i] / sum(draws[drawcount,,i])
                     drawcount <- drawcount + 1
                }
            }
            
            sources[-V,] <- sources[-V,] - alpha1 * D # add relative alpha prior counts

        }
        if(verbosity>=1){
            cat(sprintf('%d of %d depth=%4d, ', i, N, D))
            props <- colMeans(matrix(draws[,,i],ncol=V))
            prop_devs <- apply(matrix(draws[,,i], ncol=V), 2, sd)
            cat(sprintf('%.2f (%.2f)', props, prop_devs),sep='\t')
            cat('\n')
        }
    }
    
    return(draws=draws)
}

# Internal SourceTracker function to perform rarefaction analysis
"rarefy" <- function(x,maxdepth){
    if(!is.element(class(x), c('matrix', 'data.frame','array')))
        x <- matrix(x,nrow=1)
    nr <- nrow(x)
    nc <- ncol(x)

    for(i in 1:nrow(x)){
        if(sum(x[i,]) > maxdepth){
            s <- sample(nc, size=maxdepth, prob=x[i,], replace=T)
            x[i,] <- hist(s,breaks=seq(.5,nc+.5,1), plot=FALSE)$counts
        }
    }
    return(x)
}

# Internal SourceTracker function to plot pie charts
"plot.sourcetracker.pie" <- function(stresult, labels, 
        gridsize, env.colors, titlesize, indices, ...){
    V <- length(stresult$train.envs)
    for(i in indices){
        props <- apply(matrix(stresult$draws[,,i], ncol=V),2,mean)
        pie(props, labels=NA, col=env.colors, main=labels[i],cex.main=titlesize, ...)
    }
}

# Internal SourceTracker function to plot bar plots
"plot.sourcetracker.bar" <- function(stresult, labels    , 
        gridsize, env.colors, titlesize, indices, ...){
    V <- length(stresult$train.envs)
    # add extra space at top for title
    new.margins <- par('mar')
    new.margins[3] <- max(.5, new.margins[3] * 1.5)
    new.margins[2] <- .25 * titlesize
    par(mar=new.margins)
    for(i in indices){
        props <- apply(matrix(stresult$draws[,,i], ncol=V),2,mean)
        prop_devs <- apply(matrix(stresult$draws[,,i], ncol=V),2,sd)
        centers <- barplot(props, col=env.colors, main=labels[i],
                cex.main=titlesize, axes=FALSE, axisnames=FALSE,
                ylim=c(0,1.5), ...)
        sourcetracker.error.bars(centers, props, prop_devs)
        for(j in 1:4)  axis(j, at=c(-100,100),labels=FALSE)
    }
}

# Internal SourceTracker function to plot distribution plots
"plot.sourcetracker.dist" <- function(stresult, labels, 
        gridsize, env.colors, titlesize, indices, sortmethod=c('divergence', 'multilevel')[1], ...){
    # stop conditions
    if(dim(stresult$draws)[1] < 2)
        stop('Distribution plots require more than one draw.')
    V <- length(stresult$train.envs)
    for(i in indices){
        x <- stresult$draws[,,i]
        rownames(x) <- 1:nrow(x)
        if(sortmethod=='multilevel'){
            # sort by size of column
            sortby.ix <- sort(colMeans(x),index=T, dec=T)$ix
            x <- matrix(x, ncol=V)
            ix <- sortmatrix(x[,sortby.ix])
            x <- x[ix,]
        } else {
            ix <- cmdscale(jsdmatrix(x),k=1)
            ix <- sort(ix,index=T)$ix
            x <- x[ix,]
        }
        
        centers <- barplot(t(x), beside=FALSE, col=env.colors,
                space=0, border=NA, axes=FALSE, axisnames=FALSE,
                main=labels[i],cex.main=titlesize,
                ylim=c(-.05,1.05), ...)
        bounds <- c(0, min(centers)-.5, 1, max(centers)+.5)
        
        lines(c(bounds[2], bounds[4]), c(bounds[1], bounds[1]), lty=1, lwd=1)
        lines(c(bounds[2], bounds[4]), c(bounds[3], bounds[3]), lty=1, lwd=1)
        lines(c(bounds[2], bounds[2]), c(bounds[1], bounds[3]), lty=1, lwd=1)
        lines(c(bounds[4], bounds[4]), c(bounds[1], bounds[3]), lty=1, lwd=1)
    }
}

# Internal SourceTracker function to plot error bars on a bar chart
"sourcetracker.error.bars" <- function(x,centers,spread,...){
    width = min(.01, .25/length(x))
    xlim <- range(x)
    barw <- diff(xlim) * width * 4
    
    upper <- centers + spread
    lower <- centers - spread

    segments(x, upper, x, lower, lwd=1.5, ...)
    segments(x - barw, upper, x + barw, upper, lwd=1.5, ...)
    segments(x - barw, lower, x + barw, lower, lwd=1.5, ...)
}

# sorts a matrix by first column, breaks ties by the 2nd, 3rd, etc. columns
# returns row indices
"sortmatrix" <- function(x){
    # sort by last column, then 2nd-to-last, etc.
    ix <- 1:nrow(x)
    for(j in ncol(x):1){
        ixj <- sort(x[ix,j], index=T)$ix
        ix <- ix[ixj]
    }
    return(ix)
}

"jsdmatrix" <- function(x){
    d <- matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(i in 1:(nrow(x)-1)){
        for(j in (i+1):nrow(x)){
            d[i,j] <- jsd(x[i,], x[j,])
            d[j,i] <- d[i,j]
        }
    }
    return(d)
}

"jsd" <- function(p,q){
    m <- (p + q)/2
    return((kld(p,m) + kld(q,m))/2)
}

"kld" <- function(p,q){
    nonzero <- p>0 & q>0
    return(sum(p[nonzero] * log(p[nonzero]/q[nonzero])))    
}

"figlegend" <- function(){
}


# global definition of standard env colors
std.env.colors <- c(rgb(63,82,162,m=255), 
                rgb(18,130,68,m=255),
                rgb(235,20,45,m=255),
                rgb(29,189,188,m=255),
                rgb(164,57,149,m=255),
                '#885588','#6B78B4','#CC6666','#663333',
                '#47697E','#5B7444','#79BEDB','#A3C586',
                '#FFCC33','#e93E4A','#B1BDCD','#266A2E',
                '#FCF1D1','#660F57','#272B20','#003366')
