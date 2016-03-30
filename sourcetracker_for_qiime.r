# USAGE
# get help:
# Rscript sourcetracker_for_qiime.r-h 
#
# run sink predictions using QIIME taxon abundance file:
# Rscript sourcetracker_for_qiime.r-t taxa.txt -m map.txt 
#
# run leave-one-out source-sample predictions using QIIME taxon abundance file:
# Rscript sourcetracker_for_qiime.r-t taxa.txt -m map.txt -s 
#
# run sink predictions using QIIME OTU table:
# Rscript sourcetracker_for_qiime.r-i otutable.txt -m map.txt 
#
# run sink predictions using QIIME OTU table with 1000 burnins, 25 random restarts, and rarefaction depth of 100:
# Rscript sourcetracker_for_qiime.r-i otutable.txt -m map.txt -b 1000 -n 25 -r 100< sourcetracker_for_qiime.r 
#
# run sink predictions using QIIME taxon abundance file and and input file listing the sampleids to predict:
# Rscript sourcetracker_for_qiime.r-t taxa.txt -m map.txt -f sampleid_file.txt 
#
# Note: you must add the path to your SourceTracker top-level directory to an environment variable, e.g.:
# 
# echo "" >> ~/.bash_profile; echo "export SOURCETRACKER_PATH=$HOME/path/to/your/sourcetracker/repository/folder" >> ~/.bash_profile; source ~/.bash_profile
# or, if that does not work, you can add it to your .Renviron file:
# 
# echo >> $HOME/.Renviron
# echo "SOURCETRACKER_PATH=$PWD" >> $HOME/.Renviron

# load SourceTracker package
envvars <- as.list(Sys.getenv())
if(is.element('SOURCETRACKER_PATH', names(envvars))){
    sourcefile <- sprintf('%s/src/SourceTracker.r',envvars[['SOURCETRACKER_PATH']])
    source(sourcefile)
} else {
    stop("Please add SOURCETRACKER_PATH environment variable pointing to the SourceTracker top-level directory (containing 'sourcetracker_for_qiime.r')")
}


helpstr <- c(
"-i otu table: QIIME-formatted OTU table (first line is a comment line starting with '#', second line starts with '#OTU ID' followed by sample names). You must supply either this or the taxon table via '-t'.",
"-t taxon table: output from QIIME script summarize_taxa.py. You must supply either this or the otu table via '-i'.",
"-m mapfile: mapping file with an 'Env' column giving the source environments, and a 'SourceSink' column giving 'source' for source samples and 'sink' for sink samples.",
"-n number of restarts of Gibbs sampling (default 10)",
"-b number of burn-in iterations for Gibbs sampling (default 100)",
"-r rarefaction depth, 0 for none (default 1000)",
"--train_rarefaction training data rarefaction depth, 0 for none (default 1000)",
"-f sampleid_file, file containing list of samples to predict. Useful for parallel processing (default None).",
"-o outdir: output directory; default is '.'",
"-s predict source samples using leave-one-out predictions (default: FALSE)",
"--suppress_full_results suppress writing of full per-taxon predictions (default: FALSE)",
"--alpha1 alpha1: Dirichlet hyperparameter for taxa/genes in known environments (default: 1e-3)",
"--alpha2 alpha2: Dirichlet hyperparameter for taxa/genes in unknown environments (default: 1e-1)",
"--beta beta: Dirichlet hyperparameter for mixture of environments (default: 1e-2)",
"-R results file from previous run. If given, no predictions are made, only plotting and output files)",
"--tune_alphas: tune_ntrials Tune alpha values using cross-validation on the training set with this many trials (suggest at least 25); (default: 0, no tuning)",
"--color_ix: comma-separated list of color indices for alphabetical source environments",
"--eval_fit: fit_ntrials Evaluate quality of fit to the data using simulations. Ignored if less than or equal to --tune_alpha ntrials (default: 0)",
"-v: verbose output (default FALSE)")

allowed.args <- list('-i'=NULL,'-t'=NULL,'-m'=NULL,'-n'=10,'-b'=100,'-r'=1000,
                     '--train_rarefaction'=1000,
                     '-o'='.', '-v'=FALSE, '-s'=FALSE, '-f'=NULL,'-R'=NULL,
                     '--alpha1'=1e-3, '--alpha2'=1e-1, '--beta'=1e-2, '--tune_alphas'=0, '--eval_fit'=0,
                     '--color_ix'=NULL, "--suppress_full_results"=FALSE)

# Parse command-line params
# assumes that NO args are positional
# allows flags without argument
"parse.args" <- function(allowed.args,helplist=NULL){
    argv <- commandArgs(trailingOnly=TRUE)
    # print help string if requested
    if(!is.null(helpstr) && sum(argv == '-h')>0){
        cat('',helpstr,'',sep='\n')
        q(runLast=FALSE)
    }
    argpos <- 1
    for(name in names(allowed.args)){
        argpos <- which(argv == name)
        if(length(argpos) > 0){
            # test for flag without argument
            if(argpos == length(argv) || substring(argv[argpos + 1],1,1) == '-')
                allowed.args[[name]] <- TRUE
            else {
                allowed.args[[name]] <- argv[argpos + 1]
            }
        }
    }
    return(allowed.args)
}

# parse arg list
arglist <- parse.args(allowed.args)
if(is.null(arglist[['-m']])) stop('Please supply a mapping file.')

if(is.null(arglist[['-R']])){
    if(    (is.null(arglist[['-i']]) && is.null(arglist[['-t']])) 
        || (!is.null(arglist[['-i']]) && !is.null(arglist[['-t']]))) stop('Please supply a QIIME OTU table or a QIIME taxon summary.')
}
nrestarts <- as.numeric(arglist[['-n']])
burnin <- as.numeric(arglist[['-b']])
rarefaction <- as.numeric(arglist[['-r']])
train.rarefaction <- as.numeric(arglist[['--train_rarefaction']])
sourceonly <- arglist[['-s']]
outdir <- arglist[['-o']]
predictfile <- arglist[['-f']]
resultsfile <- arglist[['-R']]
alpha1 <- as.numeric(arglist[['--alpha1']])
alpha2 <- as.numeric(arglist[['--alpha2']])
beta <- as.numeric(arglist[['--beta']])
tune.alphas.ntrials <- as.numeric(arglist[['--tune_alphas']])
eval.fit.ntrials <- as.numeric(arglist[['--eval_fit']])
if(rarefaction==0) rarefaction <- NULL
env.colors <- NULL
if(!is.null(arglist[['--color_ix']])){
    env.ix <- as.numeric(strsplit(arglist[['--color_ix']],',')[[1]])
    env.colors <- std.env.colors[env.ix]
}

# create output directory
if(!is.null(outdir)) {
    dir.create(outdir,showWarnings=FALSE, recursive=TRUE)
} else outdir <- '.'

# save command that was run
sink(sprintf('%s/command.txt',arglist[['-o']]))
cat(paste(commandArgs(),collapse=' '),'\n',sep='')
sink(NULL)

# load list of samples to predict
predictlist <- NULL
if(!is.null(predictfile)){
    predictlist <- as.character(read.table(predictfile)[,1])
}

# load mapping file
map <- read.table(arglist[['-m']],sep='\t',comment='',head=T,row.names=1,check=FALSE,colClasses='character')

if(sum(colnames(map)=='Env')==0) stop("The mapping file must contain an 'Env' column naming the source environment for each sample.")
if(sourceonly){
    # if no sourcesink column, use all samples for leave-one-out
    if(sum(colnames(map)=='SourceSink')==0){
        map$SourceSink <- factor(rep('source',nrow(map)),levels=c('source','sink','ignore'))
    }
} else {
    if(sum(colnames(map)=='SourceSink')==0) stop("The mapping file must contain a 'SourceSink' column indicating 'source' or 'sink' for each sample.")
    map$SourceSink <- factor(map$SourceSink,levels=c('source','sink','ignore'))
}


if(!is.null(resultsfile)){
    load(resultsfile)
    if(sourceonly){
        filebase <- 'source_predictions'
    } else {
        filebase <- 'sink_predictions'
    }
} else {

    # load otus/taxa
    if(!is.null(arglist[['-i']])){
        otus <- read.table(arglist[['-i']],sep='\t',comment='',head=T,row.names=1,check=FALSE, skip=1)
    } else {
        otus <- read.table(arglist[['-t']],sep='\t',head=T,row.names=1,check=FALSE)
    }
    # drop "Consensus Lineage" column if present
    otus <- as.matrix(t(otus[,!grepl('Consensus|Metadata|taxonomy',colnames(otus),ignore=TRUE)]))
    

    # ensure map and data table contain the same samples in the same order
    ix <- intersect(rownames(map), rownames(otus))
    otus <- otus[ix,]
    map <- map[ix,]

    # ensure there are no "empty" samples
    
    rs <- rowSums(otus)
    num.empties <- sum(rs == 0) 
    if(num.empties > 0){
        stop(sprintf("The following %d samples have zero sequences: %s",
                num.empties,
                paste(rownames(otus)[rs==0],collapse=', '))
            )
    }

    # extract metadata from mapping file
    if(!is.null(predictlist)){
        if(sourceonly){
            # there is a specific list to predict, and we're doing source only
            # therefore ignore sink samples
            sourcesink <- map$SourceSink
            sourcesink[sourcesink=='sink'] <- NA
            names(sourcesink) <- rownames(map)
            sourcesink[predictlist] <- 'sink'
            source.ix <- which(sourcesink=='source')
            sink.ix <- which(sourcesink=='sink')
            sourceonly <- FALSE
        } else {
            source.ix <- which(map$SourceSink=='source')
            sink.ix <- which(map$SourceSink=='sink')
            names(sink.ix) <- rownames(map)[sink.ix]
            sink.ix <- sink.ix[predictlist]
        }
    } else {
        source.ix <- which(map$SourceSink=='source')
        sink.ix <- which(map$SourceSink=='sink')
    }
    envs <- map$Env

	if(length(source.ix) < 1) stop("No samples are identified as sources")
	if(length(sink.ix) < 1 && !sourceonly) stop("No samples are identified as sinks")

    # train SourceTracker object on training data
    st <- sourcetracker(otus[source.ix,,drop=F], envs[source.ix], rarefaction_depth=train.rarefaction)

    # if tuning is requested, obtain alpha values by cross-validation
    if(tune.alphas.ntrials > 0){
        verbosity <- 0
        if(arglist[['-v']]) verbosity <- 2
        tune.res <- tune.st(otus[source.ix,,drop=F], envs[source.ix], ntrials=tune.alphas.ntrials,
                            rarefaction_depth=rarefaction, verbosity=verbosity, beta=beta)
        alpha1 <- tune.res$best.alpha1
        alpha2 <- tune.res$best.alpha2
        cat(sprintf('After tuning: alpha1 = %f, alpha2 = %f, with RMSE= %.3f\n', alpha1, alpha2, tune.res$best.rmse))
        # save alphas
        sink(sprintf('%s/tuned.alphas.txt', outdir))
        cat(sprintf('# Final alpha1=%f, alpha2=%f, RMSE=%.5f, pseudo-R2=%.5f +/- %.5f\n', 
            alpha1, alpha2, tune.res$best.rmse,
            tune.res$best.pseudo.r2.sem, tune.res$best.pseudo.r2.sem))
        cat('alpha1\talpha2\tRMSE\tRMSE s.e.m.\n')
        cat(sprintf('%f\t%f\t%.5f\t%.5f\n',
            tune.res$alphas[,1], tune.res$alphas[,2], 
            tune.res$rmse, tune.res$rmse.sem))
        sink(NULL)
    }
    # plot tuning results
    if(tune.alphas.ntrials > 0 || eval.fit.ntrials > 0) {
        verbosity <- 0
        if(arglist[['-v']]) verbosity <- 2
        if(eval.fit.ntrials > tune.alphas.ntrials) {
            if(arglist[['-v']]) cat(sprintf('Evaluating fit at alpha1=%f, alpha2=%f with %d trials\n',
                                alpha1, alpha2, eval.fit.ntrials))
            results.to.plot <- eval.fit(otus[source.ix,,drop=F], envs[source.ix],
                                ntrials=eval.fit.ntrials, rarefaction_depth=rarefaction,
                                alpha1=alpha1, alpha2=alpha2, beta=beta, verbosity=verbosity-1)
            sink(sprintf('%s/eval.fit.txt',outdir))
            cat(sprintf('alpha1\t%.5f\nalpha2\t%.5f\nrmse\t%.5f\nrmse.sem\t%.5f\n',alpha1, alpha2, results.to.plot$rmse,results.to.plot$rmse.sem))
            sink(NULL)
        } else {
            results.to.plot <- tune.res$results[[which.min(tune.res$rmse)]]
        }
        plot.eval(results.to.plot,plot.type=1,
            filename=sprintf('%s/confusion_scatterplot_combined.pdf', outdir))
        plot.eval(results.to.plot,plot.type=2,
            filename=sprintf('%s/confusion_scatterplot_pairwise.pdf', outdir))
    }

    if(sourceonly){
        # Estimate leave-one-out source proportions in training data 
        results <- predict(st, rarefaction_depth=rarefaction, nrestarts=nrestarts, burnin=burnin, alpha1=alpha1, alpha2=alpha2, beta=beta, full=!arglist[['--suppress_full_results']])
        filebase <- 'source_predictions'
    } else {
        # Estimate source proportions in test data
        testdata <- otus[sink.ix,,drop=F]
        if(length(sink.ix)==1){
            testdata <- matrix(testdata,nrow=1)
            rownames(testdata) <- rownames(otus)[sink.ix]
        }
        results <- predict(st,testdata, rarefaction_depth=rarefaction, nrestarts=nrestarts, burnin=burnin, alpha1=alpha1, alpha2=alpha2, beta=beta, full=!arglist[['--suppress_full_results']])
        filebase <- 'sink_predictions'
    }
    # save full results object
    save(results,file=sprintf('%s/results.RData',outdir))
}

# save results file
sink(sprintf('%s/%s.txt', outdir, filebase))
cat('SampleID\t')
write.table(results$proportions,quote=F,sep='\t')
sink(NULL)

sink(sprintf('%s/%s_stdev.txt', outdir, filebase))
cat('SampleID\t')
write.table(results$proportions_sd,quote=F,sep='\t')
sink(NULL)

if(!arglist[['--suppress_full_results']]){

    # get average of full results across restarts
    res.mean <- apply(results$full.results,c(2,3,4),mean)

    # Get depth of each sample for relative abundance calculation
    sample.depths <- apply(results$full.results[1,,,,drop=F],4,sum)
    
    # create dir
    subdir <- paste(outdir,'full_results',sep='/')
    dir.create(subdir,showWarnings=FALSE, recursive=TRUE)
    # write each env separate file
    for(i in 1:length(results$train.envs)){
        env.name <- results$train.envs[i]
        filename.fractions <- sprintf('%s/%s_%s_contributions.txt', subdir, filebase, env.name)
        res.mean.i <- res.mean[i,,]
		# handle the case where there is only one sink sample
        if(is.null(dim(res.mean.i))) res.mean.i <- matrix(res.mean.i,ncol=1)

        # make rows be samples, columns be features
        res.mean.i <- t(res.mean.i)

        # ensure proper names
        colnames(res.mean.i) <- colnames(otus)
        rownames(res.mean.i) <- results$samplenames
        
        # calculate and save relative abundance
        res.mean.i.ra <- sweep(res.mean.i,1,sample.depths,'/')
        sink(filename.fractions)
        cat('SampleID\t')
        write.table(res.mean.i.ra,quote=F,sep='\t')
        sink(NULL)
    }
}

save.mapping.file(results, map,
        filename=sprintf('%s/map.txt',outdir),
        include.contamination.predictions=sourceonly)

# make plots
if(dim(results$draws)[2] > 1) {
    plot.types <- c('pie', 'bar', 'dist')
} else plot.types <- c('pie', 'bar')
envs <- as.factor(map[rownames(results$proportions),'Env'])
labels = sprintf('%s_%s',envs, rownames(results$proportions))
plotixs <- sort(as.numeric(envs),index=TRUE)$ix
for(plot.type in plot.types){
    # plot each env separately
    for(env in unique(envs)){
        plotixs <- which(envs == env)
        pdf(sprintf('%s/%s_%s_%s.pdf',outdir,filebase,plot.type,env),width=5,height=5)
        plot(results, type=plot.type, labels=labels, include.legend=TRUE, indices=plotixs, env.colors=env.colors)
        dev.off()    
    }
}
