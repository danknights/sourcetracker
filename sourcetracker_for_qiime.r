# USAGE
# get help:
# time R --slave --vanilla --args -h < sourcetracker_for_qiime.r 
#
# run sink predictions using QIIME taxon abundance file:
# time R --slave --vanilla --args -t taxa.txt -m map.txt < sourcetracker_for_qiime.r 
#
# run leave-one-out source-sample predictions using QIIME taxon abundance file:
# time R --slave --vanilla --args -t taxa.txt -m map.txt -s < sourcetracker_for_qiime.r 
#
# run sink predictions using QIIME OTU table:
# time R --slave --vanilla --args -i otutable.txt -m map.txt < sourcetracker_for_qiime.r 
#
# run sink predictions using QIIME OTU table with 1000 burnins, 25 random restarts, and rarefaction depth of 100:
# time R --slave --vanilla --args -i otutable.txt -m map.txt -b 1000 -n 25 -r 100< sourcetracker_for_qiime.r 
#
# run sink predictions using QIIME taxon abundance file and and input file listing the sampleids to predict:
# time R --slave --vanilla --args -t taxa.txt -m map.txt -f sampleid_file.txt < sourcetracker_for_qiime.r 
#
# Note: you must add the path to your SourceTracker.r file to your path, e.g.:
# echo "" >> ~/.bash_profile; echo "export SOURCETRACKER_PATH=$HOME/path/to/your/SourceTracker.r" >> ~/.bash_profile; source ~/.bash_profile

# load SourceTracker package
envvars <- as.list(Sys.getenv())
if(is.element('SOURCETRACKER_PATH', names(envvars))){
    source(envvars[['SOURCETRACKER_PATH']])
} else {
    stop("Please add SOURCETRACKER_PATH environment variable pointing to the source file 'SourceTracker.r'")
}


helpstr <- c(
"-i otu table: QIIME-formatted OTU table (first line is a comment line starting with '#', second line starts with '#OTU ID' followed by sample names). You must supply either this or the taxon table via '-t'.",
"-t taxon table: output from QIIME script summarize_taxa.py. You must supply either this or the otu table via '-i'.",
"-m mapfile: mapping file with an 'Env' column giving the source environments, and a 'SourceSink' column giving 'source' for source samples and 'sink' for sink samples.",
"-n number of restarts of Gibbs sampling (default 10)",
"-b number of burn-in iterations for Gibbs sampling (default 100)",
"-r rarefaction depth, 0 for none (default 1000)",
"-f sampleid_file, file containing list of samples to predict. Useful for parallel processing (default None).",
"-o outdir: output directory; default is '.'",
"-s predict source samples using leave-one-out predictions (default: FALSE)",
"--alpha1 alpha1: Dirichlet hyperparameter for known environments (default: 1e-3)",
"--alpha2 alpha2: Dirichlet hyperparameter for unknown environments (default: 1e-3)",
"-R results file from previous run. If given, no predictions are made, only plotting and output files)",
"--tune_alphas: Tune alpha values using cross-validation on the training set.",
"-v: verbose output (default FALSE)")

allowed.args <- list('-i'=NULL,'-t'=NULL,'-m'=NULL,'-n'=10,'-b'=100,'-r'=1000,
                     '-o'=NULL, '-v'=FALSE, '-s'=FALSE, '-f'=NULL,'-R'=NULL,
                     '--alpha1'=1e-3, '--alpha2'=1e-3, '--tune_alphas'=FALSE)

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
sourceonly <- arglist[['-s']]
outdir <- arglist[['-o']]
predictfile <- arglist[['-f']]
resultsfile <- arglist[['-R']]
alpha1 <- as.numeric(arglist[['--alpha1']])
alpha2 <- as.numeric(arglist[['--alpha2']])
tune.alphas <- arglist[['--tune_alphas']]
if(rarefaction==0) rarefaction <- NULL

# create output directory
if(!is.null(outdir)) {
    dir.create(outdir,showWarnings=FALSE, recursive=TRUE)
} else outdir <- '.'

# load list of samples to predict
predictlist <- NULL
if(!is.null(predictfile)){
    predictlist <- as.character(read.table(predictfile)[,1])
}

# load mapping file
map <- read.table(arglist[['-m']],sep='\t',comment='',head=T,row.names=1,check=FALSE)
if(sum(colnames(map)=='Env')==0) stop("The mapping file must contain an 'Env' column naming the source environment for each sample.")
if(sum(colnames(map)=='SourceSink')==0) stop("The mapping file must contain a 'SourceSink' column indicating 'source' or 'sink' for each sample.")
map$SourceSink <- factor(map$SourceSink,levels=c('source','sink','ignore'))

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
    otus <- as.matrix(t(otus[,!grepl('Consensus|Metadata',colnames(otus))]))

    # ensure map and data table contain the same samples in the same order
    ix <- intersect(rownames(map), rownames(otus))
    otus <- otus[ix,]
    map <- map[ix,]

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

    # train SourceTracker object on training data
    st <- sourcetracker(otus[source.ix,], envs[source.ix], rarefaction_depth=rarefaction)

    # if tuning is requested, obtain alpha values by cross-validation
    if(tune.alphas){
        tune.res <- tune.st(otus[source.ix,], envs[source.ix], rarefaction_depth=1000, verbosity=2)
        alpha1 <- tune.res$best.alpha1
        alpha2 <- tune.res$best.alpha2
        cat(sprintf('After tuning: alpha1 = %f, alpha2 = %f\n', alpha1, alpha2))
    }

    if(sourceonly){
        # Estimate leave-one-out source proportions in training data 
        results <- predict(st, rarefaction_depth=rarefaction, nrestarts=nrestarts, burnin=burnin, alpha1=alpha1, alpha2=alpha2)
        filebase <- 'source_predictions'
    } else {
        # Estimate source proportions in test data
        testdata <- otus[sink.ix,]
        if(length(sink.ix)==1){
            testdata <- matrix(testdata,nrow=1)
            rownames(testdata) <- rownames(otus)[sink.ix]
        }
        results <- predict(st,testdata, rarefaction_depth=rarefaction, nrestarts=nrestarts, burnin=burnin, alpha1=alpha1, alpha2=alpha2)
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

# also append results to a new version of the mapping file in output directory
# add columns filled with NA to the mapping file, add values where available
for(j in 1:ncol(results$proportions)){
    colname <- sprintf('Proportion_%s',colnames(results$proportions)[j])
    colname_sd <- sprintf('Proportion_SD_%s',colnames(results$proportions)[j])
    map[,colname] <- rep(NA,nrow(map))
    map[rownames(results$proportions), colname] <- results$proportions[,j]
    map[rownames(results$proportions), colname_sd] <- results$proportions_sd[,j]
}

# add this environment and its s.d. to mapping file
thisenv <- rep(NA,nrow(results$proportions))
thisenv_sd <- rep(NA,nrow(results$proportions))
for(i in 1:nrow(results$proportions)){
    sampleid <- rownames(results$proportions)[i]
    whichenv <- which(results$train.envs==map[sampleid,'Env'])
    if(length(whichenv) > 0){
        thisenv[i] <-  results$proportions[i,whichenv]
        thisenv_sd[i] <- results$proportions_sd[i,whichenv]
    }
}
map[,'Proportion_This_Env'] <- rep(NA, nrow(map))
map[,'Proportion_SD_This_Env'] <- rep(NA, nrow(map))
map[rownames(results$proportions), 'Proportion_This_Env'] <- thisenv
map[rownames(results$proportions), 'Proportion_SD_This_Env'] <- thisenv_sd

# write new mapping file
sink(sprintf('%s/map_all_samples.txt',outdir))
cat('#SampleID\t')
write.table(map,sep='\t',quote=F)
sink(NULL)

# write new mapping file
sink(sprintf('%s/map.txt',outdir))
cat('#SampleID\t')
write.table(map[rownames(results$proportions),],sep='\t',quote=F)
sink(NULL)


# make plots
if(dim(results$draws)[2] > 1) {
    plot.types <- c('pie', 'bar', 'dist')
} else plot.types <- c('pie', 'bar')
envs <- map[rownames(results$proportions),'Env']
labels = sprintf('%s_%s',envs, rownames(results$proportions))
plotixs <- sort(as.numeric(envs),index=TRUE)$ix
for(plot.type in plot.types){
    # plot each env separately
    for(env in unique(envs)){
        plotixs <- which(envs == env)
        pdf(sprintf('%s/%s_%s_%s.pdf',outdir,filebase,plot.type,env),width=8,height=8)
        plot(results, type=plot.type, labels=labels, include.legend=TRUE, indices=plotixs)
        dev.off()    
    }
}
