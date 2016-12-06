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
    '    -q <float>   Maximum q-value ([0-1]; def. 1 [all results reported])\n',
    '                   ("fdr" automatically added to -s list)              \n',
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
maxQval <- 1       # max. q-value (fdr)
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
      minCpG <- as.integer(args[i + 1])
    } else if (args[i] == '-d') {
      minDiff <- as.double(args[i + 1])
    } else if (args[i] == '-p') {
      maxPval <- as.double(args[i + 1])
    } else if (args[i] == '-q') {
      maxQval <- as.double(args[i + 1])
      dss <- c(dss, 'fdr')
    } else if (args[i] == '-h') {
      usage()
    } else {
      cat('Error! Unknown parameter:', args[i], '\n')
      usage()
    }
    i <- i + 1
  } else if (args[i] == '-h') {
    usage()
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
keep <- unique(keep)
dss <- unique(dss)

# group samples into a named list
samples <- list()
for (i in 1:length(names)) {
  if (! is.null(groups) && i <= length(groups)) {
    group <- groups[i]
  } else {
    group <- paste(names[i][[1]], collapse='_')
  }
  if (group %in% names(samples)) {
    stop('Duplicated group name: ', group)
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
comps <- c()  # group comparison strings
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
    comp <- paste(names(samples)[i], names(samples)[j], sep='->')
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
    comps <- c(comps, comp)

  }
}

# remove regions without min. number of CpG sites
Sys.time()
cat('Filtering results\n')
if (minCpG > 1) {
  res <- res[which(res$CpG >= minCpG), ]
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

# determine which rows are valid -- need only one
#   comparison to meet threshold(s)
rows <- rep(T, nrow(res))
for (n in 1:nrow(res)) {
  valid <- T
  for (comp in comps) {
    diff <- res[n, paste(comp, 'diff', sep=':')]
    if (is.na(diff) || abs(diff) < minDiff) {
      valid <- F
      next
    }
    pval <- res[n, paste(comp, 'pval', sep=':')]
    if (is.na(pval) || pval > maxPval) {
      valid <- F
      next
    }
    if (maxQval < 1) {
      qval <- res[n, paste(comp, 'fdr', sep=':')]
      if (is.na(qval) || qval > maxQval) {
        valid <- F
        next
      }
    }
    valid <- T
    break
  }
  if (! valid) {
    rows[n] <- F
  }
}

# write output results
Sys.time()
cat('Producing output file', outfile, '\n')
options(scipen=999)
for (col in c(sampleCols, groupCols)) {
  # limit results to 7 digits; reverse sign on diffs
  spl <- strsplit(col, ':')[[1]]
  if (spl[length(spl)] == 'diff') {
    res[, col] <- -round(res[, col], digits=7)
  } else {
    res[, col] <- round(res[, col], digits=7)
  }
}
write.table(res[rows, ], outfile, sep='\t', quote=F, row.names=F)
Sys.time()
