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
# maxdepth - if present, all samples with > maxdepth sequences are rarified
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
# 	if test is NULL, performs leave-one-out predictions of the training
#   samples
# burnin - number of "burn-in" passes for Gibbs sampling
# nrestarts - number of times to restart the Gibbs sampling process
# n.draws.per.restart - number of Gibbs draws to collect per restart
# delay - number of passes between draws (ignored if n.draws.per.restart is 1)
# alpha1 - prior counts of each species in the training environments,
# 	relative to the total number of sequences per training environment.
#   Higher values decrease the trust in the training data, and make the 
#   source environment distributions over taxa smoother
# alpha2 - prior counts of each species in the Unknown environment,
#   relative to the number of sequences in a given test sample.
#   Higher values make the Unknown environment smoother and less prone to
#   over-fitting a given training sample. 
# beta - prior counts of test sequences in each environment,
#   relative to the number of sequences in a given test sample.
#   Higher values cause a smoother distribution over source environments.
# maxdepth - if present, all test samples with > maxdepth sequences are rarified
# verbosity - if > 0, print progress updates while running
#
# Value: A list containing:
# draws - an array of dimension (ndraws X nenvironments X nsamples),
# 	containing all draws from gibbs sampling
# proportions - the mean proportions over all Gibbs draws for each sample
# proportions - standard deviation of the mean proportions for each sample
# train.envs - the names of the source environments
# samplenames - the names of the test samples
"predict.sourcetracker" <- function(stobj, test=NULL, 
            burnin=25, nrestarts=100, ndraws.per.restart=1, delay=10,
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
        for(i in 1:N){
            stobj.i <- sourcetracker(train[-i,], envs[-i], maxdepth=maxdepth)
            sources <- stobj.i$sources
            cat(sprintf('%3d: ',i))
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
		type=c('pie','bar','dist')[1], gridsize=NULL, ...){
    if(type=='pie') plot.sourcetracker.pie(stresult, labels=labels, gridsize=gridsize, ...)
    if(type=='bar') plot.sourcetracker.bar(stresult, labels=labels, gridsize=gridsize, ...)
    if(type=='dist') plot.sourcetracker.dist(stresult, labels=labels, gridsize=gridsize, ...)
}




######### Internal function below ####################


# Internal SourceTracker function to run Gibbs sampling
"run.gibbs" <- function(sources, test, V, T, N,
        burnin=25, nrestarts=100, ndraws.per.restart=1, delay=10,
        alpha1=0.1, alpha2=0.1, beta=0.0001, maxdepth=1000,
		verbosity=1){

    train.envs <- rownames(sources)
    ndraws <- nrestarts * ndraws.per.restart # total number of draws
    npasses <- burnin + (ndraws.per.restart-1) * delay + 1 # passes per restart

    # draws will hold all draws (ndraws x V x N)
    draws <- array(dim=c(ndraws, V, N))
    
    # rarefy samples above maxdepth if requested
    if(!is.null(maxdepth))    test <- rarefy(test, maxdepth)

    # sink samples must have integer counts
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
            envcounts <- rep(beta * D, V)
            sources[V,] <- unknown.prior # prior counts of taxa in Unknown
            sources[-V,] <- sources[-V,] + alpha1 * D # add relative alpha prior counts
            sources[V,] <- sources[V,] + alpha2 * D # add relative alpha prior counts
            
            # tally counts in envs
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

                    if(0==1){
                        cat(sprintf('%.5f',sources[,taxon]),sep='\t'); cat('\n')
                        cat(sprintf('%.5f\t%.5f\n',p_tv, p_v))
                        cat(sprintf('%.2f',p_tv * p_v),sep='\t')
                        cat('\n')
                        #cat(sprintf('env chosen: %d\n',z[ix]))
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
"plot.sourcetracker.pie" <- function(stresult, labels=NULL, 
        gridsize=NULL, env.colors=NULL, ...){
    if(is.null(env.colors)){
        std.colors <- c('blue',rgb(0,128/255,0),'red',rgb(0,191/255,191/255), rgb(191/255,0,191/255))
    }
    
    std.colors <- c('blue',rgb(0,128/255,0),'red',rgb(0,191/255,191/255), rgb(191/255,0,191/255))
    N <- dim(stresult$draws)[3]
    V <- dim(stresult$draws)[2]

    if(!is.null(gridsize) && gridsize**2 < N)
        stop(sprintf('Please choose a gridsize of at least %d.',ceiling(sqrt(N))))

    if(is.null(labels)) labels <- stresult[['samplenames']]
    if(is.null(gridsize)) gridsize <- ceiling(sqrt(N))
    
	ngridrows <- ceiling(N / gridsize)
    par(mfrow=c(ngridrows,gridsize))
    par(oma=c(1,1,1,1), mar=c(0,0,1,0))
    plot.ix <- 1
    group.ix <- 1
    for(i in 1:N){
        props <- apply(matrix(stresult$draws[,,i], ncol=V),2,mean)
        pie(props, labels=NA, col=std.colors, main=labels[i],cex.main=.8/log10(gridsize), ...)
    }
}

# Internal SourceTracker function to plot bar plots
"plot.sourcetracker.bar" <- function(stresult, labels=NULL, 
        gridsize=NULL, env.colors=NULL, ...){
    if(is.null(env.colors)){
        std.colors <- c('blue',rgb(0,128/255,0),'red',rgb(0,191/255,191/255), rgb(191/255,0,191/255))
    }
    
    std.colors <- c('blue',rgb(0,128/255,0),'red',rgb(0,191/255,191/255), rgb(191/255,0,191/255))
    N <- dim(stresult$draws)[3]
    V <- dim(stresult$draws)[2]
    N <- 20

    if(!is.null(gridsize) && gridsize**2 < N)
        stop(sprintf('Please choose a gridsize of at least %d.',ceiling(sqrt(N))))

    if(is.null(labels)) labels <- stresult[['samplenames']]
    if(is.null(gridsize)) gridsize <- ceiling(sqrt(N))

	ngridrows <- ceiling(N / gridsize)
    par(mfrow=c(ngridrows,gridsize))
    par(oma=c(1,1,1,1), mar=c(0,0,1,0))
    plot.ix <- 1
    group.ix <- 1
    for(i in 1:N){
        props <- apply(matrix(stresult$draws[,,i], ncol=V),2,mean)
        prop_devs <- apply(matrix(stresult$draws[,,i], ncol=V),2,sd)
        centers <- barplot(props, col=std.colors, main=labels[i],
                cex.main=.8/log10(gridsize), axes=FALSE, axisnames=FALSE,
                ylim=c(0,1.5), ...)
        sourcetracker.error.bars(centers, props, prop_devs)
        for(j in 1:4)  axis(j, at=c(-100,100),labels=FALSE)

    }
}

# Internal SourceTracker function to plot distribution plots
"plot.sourcetracker.dist" <- function(stresult, labels=NULL, 
        gridsize=NULL, env.colors=NULL, ...){
    if(is.null(env.colors)){
        env.colors <- c('blue',rgb(0,128/255,0),'red',rgb(0,191/255,191/255), rgb(191/255,0,191/255))
    }
 
    
    N <- dim(stresult$draws)[3]
    V <- dim(stresult$draws)[2]

    # stop conditions
    if(dim(stresult$draws)[1] < 2)
        stop('Distribution plots require more than one draw.')
    if(!is.null(gridsize) && gridsize**2 < N)
        stop(sprintf('Please choose a gridsize of at least %d.',ceiling(sqrt(N))))

    if(is.null(labels)) labels <- stresult[['samplenames']]
    if(is.null(gridsize)) gridsize <- ceiling(sqrt(N))
    
	ngridrows <- ceiling(N / gridsize)
    par(mfrow=c(ngridrows,gridsize))
    par(oma=c(1,1,1,1), mar=c(0,0,1,0))
    plot.ix <- 1
    group.ix <- 1
    for(i in 1:N){
        x <- stresult$draws[,,i]
        
        # sort by size of column
        sortby.ix <- sort(colMeans(x),index=T)$ix
        for(j in 1:ncol(x)){
            ix <- sort(x[,sortby.ix[j]],index=T)$ix
            x <- x[ix,]
        }

        barplot(t(x), beside=FALSE, col=env.colors,
                space=-.2, border=NA, axes=FALSE, axisnames=FALSE,
                main=labels[i],cex.main=.8/log10(gridsize),
				ylim=c(-.05,1.05), ...)
    }
}

# Internal SourceTracker function to plot error bars on a bar chart
"sourcetracker.error.bars" <- function(x,centers,spread,...){
    width = min(.010,.25/length(x))
    xlim <- range(x)
    barw <- diff(xlim) * width

    upper <- centers + spread
    lower <- centers - spread

    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
}
