# This example usage code will reproduce the results from the SourceTracker paper.
# Note that there is always some stochasticity in the results.

# load sample metadata
metadata <- read.table('data/metadata.txt',sep='\t',h=T,row.names=1,check=F)
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env
desc <- metadata$Description

# load OTU table
otus <- t(read.table('data/otus.txt',sep='\t',header=T,row.names=1,check=F))
otus <- otus[rownames(metadata),]

# load SourceTracker package
source('src/SourceTracker.r')

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix])

# Estimate source proportions in test data
results <- predict(st,otus[test.ix,])

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st)

# plot results
labels <- sprintf('%s %s', envs,desc)
plot(results, labels[test.ix], type='pie')
plot(results, labels[test.ix], type='bar')
plot(results, labels[test.ix], type='dist')
plot(results.train, labels[train.ix], type='pie')
plot(results.train, labels[train.ix], type='bar')
plot(results.train, labels[train.ix], type='dist')
