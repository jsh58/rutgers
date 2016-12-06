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
    '                   order as the <groupLists>; def. group names are     \n',
    '                   constructed from sample names joined by "_")        \n',
    '    -x <str>     Comma-separated list of column names in <input> to    \n',
    '                   copy to <output>, in addition to the default        \n',
    '                   ("chr", "start", "end", "CpG")                      \n',
    '    -s <str>     Comma-separated list of column names in DSS output    \n',
    '                   to copy to <output>, in addition to the default     \n',
    '                   ("mu", "diff", "pval")                              \n',
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

# default args/parameters
infile <- outfile <- groups <- NULL
names <- list()
minCpG <- 1        # min. number of CpGs
minDiff <- 0       # min. methylation difference
maxPval <- 1       # max. p-value
up <- down <- 0    # report only hyper- or hypo-methylated results
keep <- c('chr', 'start', 'end', 'CpG')  # columns of input to keep
dss <- c('chr', 'pos', 'mu', 'diff', 'pval') # columns of DSS output to keep

# get CL args
args <- commandArgs(trailingOnly=T)
if (length(args) < 4) {
  usage()
}
i <- 1
while (i < length(args) + 1) {
  if (substr(args[i], 1, 1) == '-' && i < length(args)) {
    if (args[i] == '-i') {
      infile <- args[i + 1]
    } else if (args[i] == '-o') {
      outfile <- args[i + 1]
    } else if (args[i] == '-n') {
      groups <- strsplit(args[i + 1], '[ ,]')[[1]]
    } else if (args[i] == '-x') {
      keep <- c(keep, strsplit(args[i + 1], '[ ,]')[[1]])
    } else if (args[i] == '-s') {
      dss <- c(dss, strsplit(args[i + 1], '[ ,]')[[1]])
    } else if (args[i] == '-c') {
      minCpG <- args[i + 1]
    } else if (args[i] == '-d') {
      minDiff <- args[i + 1]
    } else if (args[i] == '-p') {
      maxPval <- args[i + 1]
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
cat('Loading DSS package\n')
suppressMessages(library(DSS))
cat('Loading methylation data from', infile, '\n')
data <- read.csv(infile, sep='\t', header=T, check.names=F,
  stringsAsFactors=F)

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
    for (k in 1:ncol(data)) {
      spl <- strsplit(colnames(data)[k], '-')[[1]]
      if (length(spl) < 2) { next }
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

# perform DML pairwise tests using DSS
Sys.time()
bsdata <- makeBSseqData(frames, names(frames))
res <- data[, keep]  # results table
diffCols <- pvalCols <- c()
for (i in 1:(length(samples)-1)) {
  for (j in (i+1):length(samples)) {

    # perform DML test
    cat('Comparing group "', names(samples)[i],
      '" to group "', names(samples)[j], '"\n', sep='')
    dml <- DMLtest(bsdata, group1=samples[[ i ]], group2=samples[[ j ]])

    # remove extraneous columns
    for (col in colnames(dml)) {
      if (! col %in% dss && ! substr(col, 1, nchar(col)-1) %in% dss) {
        dml[, col] <- NULL
      }
    }

    # add results to res
    start <- ncol(res) + 1
    res <- suppressWarnings( merge(res, dml,
      by.x=c('chr', 'start'), by.y=c('chr', 'pos'),
      all.x=T) )

    # add groups to column names
    comp <- paste(names(samples)[j], names(samples)[i], sep='->')
    for (k in start:ncol(res)) {
      col <- colnames(res)[k]
      if (substr(col, nchar(col), nchar(col)) == '1') {
        colnames(res)[k] <- paste(names(samples)[i],
          substr(col, 1, nchar(col)-1), sep=':')
      } else if (substr(col, nchar(col), nchar(col)) == '2') {
        colnames(res)[k] <- paste(names(samples)[j],
          substr(col, 1, nchar(col)-1), sep=':')
      } else {
        colnames(res)[k] <- paste(comp, col, sep=':')
      }
    }

  }
}

# for repeated columns, average the values
repCols <- c()
sampleCols <- c()
groupCols <- c()
for (i in 1:ncol(res)) {
  if (i %in% repCols) { next }

  # find duplicated columns
  repNow <- c(i)
  if (i < ncol(res)) {
    for (j in (i + 1):ncol(res)) {
      if (colnames(res)[i] == colnames(res)[j]) {
        repNow <- c(repNow, j)
      }
    }
  }

  # average duplicates
  if (length(repNow) > 1) {
    res[, i] <- rowMeans(res[, repNow], na.rm=T)
    repCols <- c(repCols, repNow)
    sampleCols <- c(sampleCols, colnames(res)[i])
  } else if (! colnames(res)[i] %in% keep) {
    groupCols <- c(groupCols, colnames(res)[i])
  }
}
res <- res[, c(keep, sampleCols, groupCols)]

# determine which rows are valid
mat <- matrix(T, nrow=nrow(res), 3)
if (minCpG > 1) {
  mat[, 1] <- res[, 'CpG'] >= minCpG
}
if (minDiff > 0) {
  cols <- grep(':diff$', colnames(res), value=T)
  if (length(cols) > 1) {
    mat[, 2] <- rowSums(abs(res[, cols]) >= minDiff, na.rm=T) > 0
  } else {
    mat[, 2] <- ! is.na(res[, cols]) & abs(res[, cols]) >= minDiff
  }
}
if (maxPval < 1) {
  cols <- grep(':pval$', colnames(res), value=T)
  mat[, 3] <- rowSums(abs(res[, cols]) <= maxPval, na.rm=T) > 0
}
print(head(res))
print(head(mat))
stop()

# write output results
write.table(res, outfile, sep='\t', quote=F, row.names=F)
#write.table(format(res, digits=7, scientific=F), outfile, sep='\t', quote=F, row.names=F)
Sys.time()
stop()


###################################################################


# write output header
Sys.time()
cat('Producing output file', outfile, '\n')
for (i in names(dmls)) {
  for (j in names(dmls[[ i ]])) {

    # fix column names
    cols <- ncol(dmls[[ i ]][[ j ]])
#print(head(dmls[[i]][[j]]))
#print(cols)
#stop()
    for (k in 1:(ncol(res) - cols)) {
      idx <- ncol(res) - cols - 1 + k
      col <- colnames(res)[idx]
      if (substr(col, nchar(col), nchar(col)) == '1') {
      #if (k % 5 == 1) {
        colnames(res)[idx] <- paste(i, substr(col, 1, nchar(col)-1), sep=':')
      #} else if (k % 5 == 2) {
      } else if (substr(col, nchar(col), nchar(col)) == '2') {
        colnames(res)[idx] <- paste(j, substr(col, 1, nchar(col)-1), sep=':')
      } else {
        colnames(res)[idx] <- paste(comp, col, sep=':')
      }
    }
#print(colnames(data))
print(head(res))
stop()
  }
}
print(head(data))
Sys.time()
stop()



header <- colnames(data)[1:4]
for (i in names(samples)) {
  header <- c(header, paste(i, 'mu', sep=':'))
}
for (i in names(dmls)) {
  for (j in names(dmls[[ i ]])) {
    comp <- paste(i, j, sep='->')
    for (val in c('diff', 'pval')) {
      header <- c(header, paste(comp, val, sep=':'))
    }
  }
}
resMat <- matrix()
colnames(resMat, header)
#cat(header, file=outfile, sep='\t')
#cat('\n', file=outfile, append=T)

# write output results
for (key in head(row.names(data), n=3)) {
#for (key in row.names(data)) {
  if (data[key, 'CpG'] < minCpG) {
    next
  }
  valid <- F
  res <- c()
  mus <- list()
  for (i in names(dmls)) {
    for (j in names(dmls[[ i ]])) {
      mus[[ i ]] <- c(mus[[ i ]], dmls[[ i ]][[ j ]][ key, 'mu1' ])
      mus[[ j ]] <- c(mus[[ j ]], dmls[[ i ]][[ j ]][ key, 'mu2' ])
      diff <- dmls[[ i ]][[ j ]][ key, 'diff' ]
      pval <- dmls[[ i ]][[ j ]][ key, 'pval' ]
      if (! is.na(diff) && abs(diff) >= minDiff
          && ! is.na(pval) && pval <= maxPval) {
        valid <- T
      }
      res <- c(res, -diff, pval)
      #cat(i, j, key, diff, pval)
      #print(dmls[[ i ]][[ j ]][key,])
    }
  }
  if (valid) {
    # average mu values for each group
    mu <- c()
    for (i in names(samples)) {
      mu <- c(mu, mean(mus[[ i ]], na.rm=T))
    }
    rbind(resMat, c(unlist(data[key, 1:4], use.names=F), mu, res))
    #cat(unlist(data[key, 1:4], use.names=F), mu, res,
    #  file=outfile, sep='\t', append=T)
    #cat('\n', file=outfile, append=T)
  }
}
write.table(resMat, outfile, sep='\t', quote=F, row.names=F)
Sys.time()
