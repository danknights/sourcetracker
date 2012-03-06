# This runs SourceTracker on the original "contamination" data set
# (data included in 'data' folder)

# load sample metadata
metadata <- read.table('data/metadata.txt',sep='\t',h=T,row.names=1,check=F,comment='')

# load OTU table
# This 'read.table' command is designed for a 
# QIIME-formatted OTU table.
# namely, the first line begins with a '#' sign
# and actually _is_ a comment; the second line
# begins with a '#' sign but is actually the header
otus <- read.table('data/otus.txt',sep='\t', header=T,row.names=1,check=F,skip=1,comment='')
otus <- t(as.matrix(otus))

# extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
    message <- paste(sprintf('Error: there are %d sample ids in common '),
                    'between the metadata file and data table')
    stop(message)
}

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env
if(is.element('Description',colnames(metadata))) desc <- metadata$Description


# load SourceTracker package
source('src/SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix])

# Estimate source proportions in test data
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

# plot results
labels <- sprintf('%s %s', envs,desc)
plot(results, labels[test.ix], type='pie')

# other plotting functions
# plot(results, labels[test.ix], type='bar')
# plot(results, labels[test.ix], type='dist')
# plot(results.train, labels[train.ix], type='pie')
# plot(results.train, labels[train.ix], type='bar')
# plot(results.train, labels[train.ix], type='dist')

# plot results with legend
# plot(results, labels[test.ix], type='pie', include.legend=TRUE, env.colors=c('#47697E','#5B7444','#CC6666','#79BEDB','#885588'))