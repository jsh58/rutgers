# JMG 12/1/16

# Finding DMRs using the R package 'DSS'.

usage <- function() {
  cat('Usage: Rscript findDMRs.r  [options]  -i <input>  -o <output>  \\    \n',
    '         <groupList1>  <groupList2>  [...]                             \n',
    '    <groupList>  Comma-separated list of sample names (as found in     \n',
    '                   the header of <input>, with "-N" and "-X")          \n',
    '    <input>      File listing genomic regions and methylation results  \n',
    '                   (output from combineRegions5.py)                    \n',
    '  Options:                                                             \n',
    '    -n <str>     Comma-separated list of group names (in the same      \n',
    '                   order as the <groupLists>)                          \n',
    '    -c <int>     Minimum number of CpGs in a region (def. 1)           \n',
    '    -d <float>   Minimum methylation difference between sample groups  \n',
    '                   ([0-1]; def. 0 [all results reported])              \n',
    '    -p <float>   Maximum p-value ([0-1]; def. 1 [all results reported])\n',
#    -up          Report only regions hypermethylated in later group
#    -down        Report only regions hypomethylated in later group
#    -dna         Report regions whose diff is 'NA' (occurs when one
#                   group has no methylation data)
#    -pna         Report regions whose p-value is 'NA' (occurs when one
#                   group does not have multiple data points)

    '\n')

  q()
}


# get CL args
args <- commandArgs(trailingOnly=T)
infile <- outfile <- groups <- NULL
names <- list()
i <- 1
while (i < length(args) + 1) {
  if (substr(args[i], 1, 1) == '-' && i < length(args)) {
    if (args[i] == '-i') {
      infile <- args[i + 1]
    } else if (args[i] == '-o') {
      outfile <- args[i + 1]
    } else if (args[i] == '-n') {
      groups <- strsplit(args[i + 1], '[ ,]')[[1]]
    } else if (args[i] == '-h') {
      usage()
    }
    i <- i + 1
  } else {
    names <- c(names, strsplit(args[i], '[ ,]'))
  }
  i <- i + 1
}
if (is.null(infile) || is.null(outfile)) {
  cat('Error! Must specify input and output files\n')
  usage()
}
if (length(names) < 2) {
  cat('Error! Must specify at least two groups of samples\n')
  usage()
}

# group samples into a named list
samples <- list()
for (i in 1:length(names)) {
  if (! is.null(groups) && i <= length(groups)) {
    group <- groups[i]
  } else {
    group <- paste(names[i][[1]], collapse='_')
  }
  samples[[ group ]] <- names[i][[1]]
}

# load data, and reformat to meet the DSS expectations
library(DSS)
cat('Loading data from', infile, '\n')
data <- read.csv(infile, sep='\t', header=T, check.names=F)

# determine columns for samples
idx <- list()
idx[[ 'N' ]] <- list()
idx[[ 'X' ]] <- list()
for (i in names(samples)) {
  idx[[ 'N' ]][[ i ]] <- rep(NA, length(samples[[ i ]]))
  idx[[ 'X' ]][[ i ]] <- rep(NA, length(samples[[ i ]]))
}
for (i in names(samples)) {
  for (j in 1:length(samples[[ i ]])) {
    for (k in 5:ncol(data)) {
      spl <- strsplit(colnames(data)[k], '-')[[1]]
      if (spl[-length(spl)] == samples[[ i ]][ j ]) {
        if (spl[length(spl)] == 'N') {
          idx[[ 'N' ]][[ i ]][ j ] <- k
        } else if (spl[length(spl)] == 'X') {
          idx[[ 'X' ]][[ i ]][ j ] <- k
        }
      }
    }
    if ( is.na(idx[[ 'N' ]][[ i ]][ j ])
        || is.na(idx[[ 'X' ]][[ i ]][ j ]) ) {
      stop('Missing information from input file ', infile, ':\n',
        '  For sample "', samples[[ i ]][ j ], '", need both "',
        samples[[ i ]][ j ], '-N" and "',
        samples[[ i ]][ j ], '-X" columns')
    }
  }
}

# create data frames for each sample
frames <- list()
for (i in names(samples)) {
  for (j in 1:length(samples[[ i ]])) {
    tab <- data.frame('chr'=data$chr, 'pos'=data$start,
      'N'=data[, idx[[ 'N' ]][[ i ]][ j ] ],
      'X'=data[, idx[[ 'X' ]][[ i ]][ j ] ])
    frames[[ samples[[ i ]][ j ] ]] <- tab
  }
}
#for (l in 1:length(frames)) {
#  print(names(frames[l]))
#  print(head(frames[[names(frames)[l]]]))
#}
#stop()

bsdata <- makeBSseqData(frames, names(frames))

# perform DML pairwise test
for (i in 1:(length(samples)-1)) {
  for (j in (i+1):length(samples)) {
    cat('Comparing group "', names(samples)[i],
      '" to group "', names(samples)[j], '"\n', sep='')
    dml <- DMLtest(bsdata, group1=samples[[ i ]], group2=samples[[ j ]])
    write.table(dml, outfile, sep='\t', quote=F, row.names=F, append=T)
  }
}
