# JMG 12/1/16

# Finding DMRs using the R package 'DSS'.

usage <- function() {
  cat('Usage: Rscript findDMRs.r  [options]  -i <input>  -o <output>  \\
             <groupList1>  <groupList2>  [...]
    <groupList>  Comma-separated list of sample names (as found in
                   the header of <input>, with "-N" and "-X") --
                   at least two such lists must be provided
    <input>      File listing genomic regions and methylation counts
                   (output from combineRegions.py)
    <output>     File listing genomic regions and methylation results:
                   - for each group: average methylation fraction
                   - for each group->group comparison: methylation
                     difference and p-value
                   Other DSS results can be added using -s (see below)
  Options:
    -m <str>     Comma-separated list of group names (in the same
                   order as the <groupLists>; def. group names are
                   constructed from sample names joined by "_")
    -k <str>     Comma-separated list of column names in <input> to
                   copy to <output>, in addition to the default
                   ("chr", "start", "end", "CpG")
    -s <str>     Comma-separated list of column names in DSS output
                   to copy to <output>, in addition to the default
                   ("mu", "diff", "pval")
    -c <int>     Minimum number of CpGs in a region (def. 1)
    -d <float>   Minimum methylation difference between sample groups
                   ([0-1]; def. 0 [all results reported])
    -p <float>   Maximum p-value ([0-1]; def. 1 [all results reported])
    -q <float>   Maximum q-value ([0-1]; def. 1 [all results reported])
                   ("fdr" automatically added to -s list)
    -up          Report only regions hypermethylated in later group
    -down        Report only regions hypomethylated in later group
    -na          Report regions even if DSS results (diff / p-value /
                   q-value) are "NA" (ignoring -d, -p, -q, -up, -down)
')
  q()
}

# default args/parameters
infile <- outfile <- groups <- NULL
names <- list()        # list of sample names
minCpG <- 1            # min. number of CpGs
minDiff <- 0           # min. methylation difference
maxPval <- 1           # max. p-value
maxQval <- 1           # max. q-value (fdr)
na <- up <- down <- F  # report NA/hyper-/hypo- methylated results
keep <- c('chr', 'start', 'end', 'CpG')  # columns of input to keep
dss <- c('chr', 'pos', 'mu', 'diff', 'pval') # columns of DSS output to keep

# get CL args
args <- commandArgs(trailingOnly=T)
if (length(args) < 6) {
  usage()
}
i <- 1
while (i < length(args) + 1) {
  if (substr(args[i], 1, 1) == '-') {
    if (args[i] == '-na') {
      na <- T
    } else if (args[i] == '-up') {
      up <- T
    } else if (args[i] == '-down') {
      down <- T
    } else if (args[i] == '-h') {
      usage()
    } else if (i < length(args)) {
      if (args[i] == '-i') {
        infile <- args[i + 1]
      } else if (args[i] == '-o') {
        outfile <- args[i + 1]
      } else if (args[i] == '-m') {
        groups <- strsplit(args[i + 1], '[ ,]')[[1]]
      } else if (args[i] == '-k') {
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
      } else {
        cat('Error! Unknown parameter:', args[i], '\n')
        usage()
      }
      i <- i + 1
    } else {
      cat('Error! Unknown parameter with no arg:', args[i], '\n')
      usage()
    }
  } else {
    names <- c(names, strsplit(args[i], '[ ,]'))
  }
  i <- i + 1
}

# check for errors, repeated columns
if (is.null(infile) || is.null(outfile)) {
  cat('Error! Must specify input and output files\n')
  usage()
}
input <- tryCatch( file(infile, 'r'), warning=NULL,
  error=function(e) stop('Cannot read input file ',
  infile, '\n', call.=F) )
output <- tryCatch( file(outfile, 'w'), warning=NULL,
  error=function(e) stop('Cannot write to output file ',
  outfile, '\n', call.=F) )
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

# load DSS
cat('Loading DSS package\n')
suppressMessages(library(DSS))

# load data
cat('Loading methylation data from', infile, '\n')
data <- read.csv(input, sep='\t', header=T, check.names=F)
for (k in keep) {
  if (! k %in% colnames(data)) {
    stop('Missing column in input file ', infile, ': ', k, '\n')
  }
}
if (colnames(data)[1] != 'chr') {
  stop('Improperly formatted input file ', infile, ':\n',
    '  Must have "chr" as first column\n')
}

# sort chromosome names by number/letter
level <- levels(data$chr)
intChr <- strChr <- intLev <- strLev <- c()
for (i in 1:length(level)) {
  if (substr(level[i], 1, 3) == 'chr') {
    sub <- substr(level[i], 4, nchar(level[i]))
    if (!is.na(suppressWarnings(as.integer(sub)))) {
      intChr <- c(intChr, as.numeric(sub))
    } else {
      strChr <- c(strChr, sub)
    }
  } else {
    sub <- level[i]
    if (!is.na(suppressWarnings(as.integer(sub)))) {
      intLev <- c(intLev, as.numeric(sub))
    } else {
      strLev <- c(strLev, sub)
    }
  }
}
# put numeric chroms first, then strings
chrOrder <- c(paste('chr', levels(factor(intChr)), sep=''),
  levels(factor(intLev)),
  paste('chr', levels(factor(strChr)), sep=''),
  levels(factor(strLev)))
data$chr <- factor(data$chr, levels=chrOrder)

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

# for each sample, create data frames to meet DSS requirements
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
bsdata <- makeBSseqData(frames, names(frames))
res <- data[, keep]  # results table
comps <- c()  # group comparison strings
for (i in 1:(length(samples)-1)) {
  for (j in (i+1):length(samples)) {

    # perform DML test
    cat('Comparing group "', names(samples)[i],
      '" to group "', names(samples)[j], '"\n', sep='')
    dml <- DMLtest(bsdata, group1=samples[[ i ]], group2=samples[[ j ]])

    # make sure necessary columns are present, remove extraneous
    col <- colnames(dml)
    if (any( ! dss %in%col & ! paste(dss, '1', sep='') %in%col)) {
      stop('Missing column(s) from DSS result: ',
        paste(dss[! dss %in%col & ! paste(dss, '1', sep='') %in%col],
        collapse=', '), '\n')
    }
    dml[, ! col %in% dss & ! substr(col, 1, nchar(col)-1) %in% dss] <- NULL

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
rows <- logical(length=nrow(res))
if (na) {
  rows <- rep(T, length(rows))
} else {
  # for each comparison, check diff, pval, fdr (optional)
  for (n in 1:nrow(res)) {
    valid <- T
    for (comp in comps) {
      diff <- res[n, paste(comp, 'diff', sep=':')]
      if (is.na(diff) || abs(diff) < minDiff ||
          (up && diff > 0) || (down && diff < 0)) {
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
    rows[n] <- valid
  }
}

# write output results
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
write.table(res[rows, ], output, sep='\t', quote=F, row.names=F)
