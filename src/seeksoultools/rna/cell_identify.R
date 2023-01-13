#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-i", "--inputdir"), type="character",
                help="path of input(full path of raw_feature_bc_matrix)"),
    make_option(c("-o", "--outputdir"), type="character",
                help="path of output, default=filtered_feature_bc_matrix under the sample parent dir of inputdir"),
    make_option(c("-e", "--expect_cells"), type="integer",
                help="expect recovered cells"),
    make_option(c("-p", "--pvalue"), type="double", default='0.01',
                help="pvalue cutoff")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$inputdir)){
  print_help(opt_parser)
  stop("No inputfile/expect_cells set", call.=FALSE)
}

if (!dir.exists(opt$inputdir)) {
    message("input dir does not exist")
    stop
}
inputdir = opt$inputdir

outputdir = ''
if(!is.null(opt$outputdir)){
    outputdir = opt$outputdir
}else{
    outputdir = paste0(dirname(opt$inputdir), '/filtered_feature_bc_matrix')
}

recovered_cells = 3000
if (opt$expect_cells){
    recovered_cells =  opt$expect_cells
}

pvalue = opt$pvalue

print(paste0("inputdir:", opt$inputdir, "    outputdir:", opt$outputdir,
             "    expect_cells: ", recovered_cells, "    pvalue:", pvalue))

suppressMessages({
    library(DropletUtils)
    library(Matrix)
    library(edgeR)
})

rmultinomial <- function(n, size, prob)
#ISALIAS dmultinomial
#--------------------------------------------
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rmultinomial")
  
  if(is.vector(prob) || (dim(prob)[1]) == 1) {
    if(length(size)==1) return(t(rmultinom(n,size,prob)))      # classical
    prob <- matrix(prob,nrow=1)
    }

  nrp <- nrow(prob)
  mnr <- min( max(nrp ,length(size)), n)
  ss  <- rep(size,length.out=mnr)              # recycling size
  
  if(nrp != mnr) prob <- matrix(t(prob),ncol=ncol(prob),nrow=mnr,byrow=TRUE)    # recycling prob

  n1 <- n%/%mnr
  n2 <- n%%mnr

  res <- sapply(1:mnr,function(x) rmultinom(n1,ss[x],prob[x,]))
  res <- matrix(res,ncol=ncol(prob),byrow=TRUE)
  index <- as.vector(matrix(1:(mnr*n1),ncol=mnr,byrow=TRUE))
  res <- res[index,]

  if (n2 != 0){
    res <- rbind(res,t(sapply(1:n2,function(x) rmultinom(1,ss[x],prob[x,]))))
    }
  return(res)
}
                              
dmultinomial <-function(x, size = NULL, prob, log = FALSE){
    if(length(x) == 0) return(numeric(0))
    if(is.vector(prob)) prob <- t(prob)
    if(is.vector(x)) x <- t(x)
    K <- ncol(prob)
    maxn <- max(nx <- nrow(x), np <- nrow(prob))
    if (ncol(x) != K) stop("x[] and prob[] must be equal length vectors or equal col matrix.")
    if(nx != maxn) x <- matrix(t(x), ncol=K, nrow=maxn, byrow=TRUE)
    if(np != maxn) prob <- matrix(t(prob), ncol=K, nrow=maxn, byrow=TRUE)

    N <- apply(x,1,sum)
    if (is.null(size)) size <- N
    size <- rep(size, length = maxn)
    if (any(size != N)) stop("at least one size != sum(x), i.e. one is wrong")

    if (any(prob < 0) || any((s <- apply(prob,1,sum)) == 0))
        stop("probabilities cannot be negative nor all 0.")
    prob <- prob/s
    if (any(as.integer(x + 0.5) < 0)) stop("'x' must be non-negative")

    r <- sapply(1:maxn, function(y) {
                                        xx <- as.integer(x[y,]+.5)
                                        pp <- prob[y,]
                                        i0 <- pp == 0
                                        if (any(i0)) {
                                            if (any(xx[i0] != 0)) return(0)
                                            xx <- xx[!i0]
                                            pp <- pp[!i0]
                                                    }
                                        return(lgamma(size[y] + 1) + sum(xx * log(pp) - lgamma(xx + 1)))
                                        }
                ) 
    if (log) return(r) else return(exp(r))
}

find_top_umi_barcode <- function(barcode_counts, recovered_cells,  
                                 ORDMAG_NUM_BOOTSTRAP_SAMPLES, quantile)
{

    nonzero_index <- which(barcode_counts  > 0)
    bc_counts_nonzero <- barcode_counts [nonzero_index]
    bc_counts_nonzero_length <- length(bc_counts_nonzero)
    baseline_idx = min(round(recovered_cells * (1 - quantile)), bc_counts_nonzero_length -1 )
    find_within_ordmag <- function(bc_counts_nonzero, baseline_idx)
    {
        x = sample(bc_counts_nonzero, bc_counts_nonzero_length, replace=TRUE)
        x_sort = sort(x)
        baseline = x_sort[length(x_sort)-baseline_idx+1]
        cutoff = max(1, round(0.1 * baseline))
        # Return the index corresponding to the cutoff in descending order
        indices_cutoff <- (min(which(x_sort >= cutoff)) - 1)
        return (length(x) - indices_cutoff)
    }

    top_n_boot <- c()
    for (i in 1:ORDMAG_NUM_BOOTSTRAP_SAMPLES) {
        test <- find_within_ordmag(bc_counts_nonzero, baseline_idx)
        top_n_boot <- c(top_n_boot, find_within_ordmag(bc_counts_nonzero, baseline_idx))
    }
    top_n = round(mean(top_n_boot))
    top_n_indices <- nonzero_index[order(bc_counts_nonzero)[(bc_counts_nonzero_length - top_n + 1):bc_counts_nonzero_length]]
    
    
    return(top_n_indices)
}


find_ambient_background <- function(matrix_counts, barcode_counts,  N_PARTITIONS){    
    N_PARTITIONS_HALF = as.integer(N_PARTITIONS/2)
    # select ambient cells
    #barcode_partitions_indices <- which(barcode_sort_indices >= p_reverse & barcode_sort_indices <= p_half_reverse & barcode_sort_indices <= ordmag_reverse )
    false_all <- rep(FALSE, length(barcode_counts))
    barcode_sort_indices <- order(barcode_counts) # get indice of sorted barcode_counts  
    barcode_sort_indices_rev <- rev(barcode_sort_indices)
    false_all[barcode_sort_indices_rev[(N_PARTITIONS_HALF+1):N_PARTITIONS]] = TRUE
    ambient = false_all
    barcode_partitions_indices <- which(ambient)
    return(ambient)
}

compute_ambient_prop <- function(nonzero_transcript_matrix){
    ncells <- ncol(nonzero_transcript_matrix)
    # Computing the average profile from the ambient cells.
    umi.sum <- as.integer(round(colSums(nonzero_transcript_matrix)))
    ambient.cells <- nonzero_transcript_matrix[,ambient]
    ambient.prof <- rowSums(ambient.cells)
    ambient_prop <- edgeR::goodTuringProportions(ambient.prof)
    return(ambient_prop)
} 


find_ambient_candidate <- function(barcode_counts, top_n_indices, min_umis_nonambient, N_CANDIDATE_BARCODES){
    # final[umi > 500 & partitions barcode] 
    median_initial_umis <- median(barcode_counts[top_n_indices])
    min_umis = as.integer(max(min_umis_nonambient, round(ceiling(median_initial_umis * min_umi_frac_of_median))))

    ambient_candidate <- barcode_counts  >= min_umis
    ambient_candidate[top_n_indices] = FALSE
    # reduce the candidate num to N_CANDIDATE_BARCODES if candidate_num > N_CANDIDATE_BARCODES
    candidate_num = length(which(ambient_candidate))
    if (candidate_num > N_CANDIDATE_BARCODES){
        ambient_candidate_indices = which(ambient_candidate)
        ambient_candidate_indices_order = order(barcode_counts[ambient_candidate_indices])
        candidate_start <- length(ambient_candidate_indices) - N_CANDIDATE_BARCODES + 1
        candidate_end <- length(ambient_candidate_indices)
        ambient_candidate_indice <- ambient_candidate_indices[ambient_candidate_indices_order[candidate_start:candidate_end]]
        ambient_candidate <- rep(FALSE, length(barcode_counts))
        ambient_candidate[ambient_candidate_indice] = TRUE
    }
    print(paste0('Number of candidate bcs: ', length(ambient_candidate)))
    print(paste0('Range candidate bc umis: ', min(barcode_counts[ambient_candidate]), '-', max(barcode_counts[ambient_candidate])))
    return(ambient_candidate)
}

compute_simulate_pvalues <- function(barcode_counts, ambient_prop,
                                     barcode_umi_candidate_sorted, num_sims, jump, n_sample_feature_block)
{
    ambient_prop <- ambient_prop[,1]
    
    umi_per_bc_candidate <- barcode_counts[barcode_umi_candidate_sorted]

    unique_n <- unique(umi_per_bc_candidate)

    #loglk <- c(rep(c(rep(0, num_sims)), length(unique_n)))
    loglk <- matrix(0, ncol=num_sims, nrow=length(unique_n))

    num_all_n = max(unique_n) - min(unique_n)

    print(paste0('Number of distinct N supplied: ', length(unique_n))) 
    print(paste0('Range of N: ', num_all_n))
    print(paste0('Number of features: ',  length(ambient_prop)))

    sampled_features <- sample(length(ambient_prop), n_sample_feature_block, prob=ambient_prop, replace=TRUE)

    k = 1

    log_ambient_prop <- log(ambient_prop)

    for (sim_idx in (1:num_sims)){
        curr_counts <- rmultinomial (1, unique_n[1],  ambient_prop)
        curr_loglk = dmultinomial (curr_counts, unique_n[1], ambient_prop, TRUE)
        loglk[1, sim_idx] = curr_loglk
        for (i in 2:length(unique_n)){
            step <- unique_n[i] - unique_n[i-1]
            #if (step > jump){
            if (1 > 2){
                # Instead of iterating for each n, sample the intermediate ns all at once
                curr_counts = curr_counts +  rmultinom(length(ambient_prop), unique_n[1],  ambient_prop)
                curr_loglk = dmultinom(curr_counts, unique_n[1], ambient_prop)
            }else{
                # Iteratively sample between the two distinct values of n
                for (n in ((unique_n[i-1]+1) : (unique_n[i]))){
                    j <- sampled_features[k]
                    k = k + 1
                    if (k >= n_sample_feature_block){
                        # Amortize this operation
                        sampled_features <- sample(length(ambient_prop), n_sample_feature_block, prob=ambient_prop, replace=TRUE)
                        k = 1
                    }
                    curr_counts[j] <- curr_counts[j] + 1
                    curr_loglk =  curr_loglk + (log_ambient_prop[j] + log(n/curr_counts[j]))
              }
            }
            loglk[i, sim_idx] = curr_loglk
        }
    }
    return(loglk)
}


compute_ambient_candidate_pvalues <- function(candidate_loglk_sorted,  sim_loglk, unique_n,  umi_per_bc_candidate){
    sim_n_idx = match(umi_per_bc_candidate, unique_n)
    num_sims = dim(sim_loglk)[2]
    num_barcodes = length(umi_per_bc_candidate)
    pvalues = rep(0, num_barcodes)
    for(i in 1:num_barcodes){
        num_lower_loglk <- sum(sim_loglk[sim_n_idx[i],] < candidate_loglk_sorted[i])
        pvalues[i] <- (1 + num_lower_loglk) / (1 + num_sims)    
    }  
    return(pvalues)
}
                             
output_10x_format_results <-function(outputdir,  matrix_read, matrix_counts, cell_barcode){
    # get cell only data
    gi <- rownames(matrix_read)
    gs <- rowData(matrix_read)[,2:2]
    bc <- colData(matrix_read)[,2][cell_barcode]
    umi.counts <- matrix_counts[,cell_barcode]     
    
    # output 10X format filtered_feature_bc_matrix
    write10xCounts(outputdir, umi.counts, gene.id=gi, gene.symbol=gs, barcodes=bc, version='3', overwrite=TRUE)

    ## add line 2 of matrix.mtx.gz
    ed.matrix <- paste(outputdir, 'matrix.mtx.gz', sep="/")
    middle.matrix <- paste(outputdir, 'matrix.mtx.gz.format', sep="/")
    tenx.matrix <- paste(outputdir, 'matrix.mtx.gz', sep="/")
}                             

#read_10x_and_compare_plot 
read_10x_and_compare_plot <- function(matrix_counts, ambient_prop, barcode_umi_candidate_sorted, sim_logl, pvalues, p_adjust){
    obs_loglk_10x <- read.csv('obs_loglk.txt')[,1]
    barcode_10x_ambient <- read.csv('10x_eval_bcs.txt', header = FALSE)
    pvalue_10x_adj  <- read.csv('10x_pvalues_adj.txt', header = FALSE)[,1]
    pvalue_10x <- read.csv('pvalues.txt', header = FALSE)[,1]
    barcode_10x_ambient <- barcode_10x_ambient[,1] + 1
    indices_10x_ambient <- read.csv('/PROJ/development/maxingyong/EmptyDrops/NP-6/10x_use_bcs.txt')[,1] + 1
    indices_10x_ordmag <- read.csv('10x_barcode_index_ordmag.txt')[,1] + 1
    indices_10x_candidate <- read.csv('10x_eval_bcs.txt')[,1] + 1
    filter_10x_ambient <- read.csv('10x_nonambient.indices.txt')
    indices_10x_nonambient<- indices_10x_candidate[which(filter_10x_ambient == 'True')]
    indices_10x_cell <- c(indices_10x_ordmag, indices_10x_nonambient)
   # print(length(indices_10x_cell))
    indices_ambient_50k <- read.csv('10x_use_bcs.txt', header=FALSE)
    features <- read.table('raw_feature_bc_matrix/features.tsv.gz', sep='\t')
    bc_order <- read.csv('bc_order.txt', header=FALSE)
    #print(bc_order)
    print('bc_order output')

    #compare with 10x and plot
    #plot ambient_profile_p overlap with 10X
    all_gene <- rownames(matrix_counts)
    gene_indices_10x <- (read.csv('eval_features.txt'))[,1] + 1
    ed_gene <- all_gene[gene_indices_10x]
    ambient_profile_p_10x <- read.csv('ambient_profile_p.txt')[,1]
    ambient_prop <- as.numeric(ambient_prop[,1])
    plot(ambient_prop[2:length(ambient_prop)], ambient_profile_p_10x )

    #plot obs pvalue
    candidate_intersect_indices <-  intersect(barcode_10x_ambient, barcode_umi_candidate_sorted)
    pvalue_indices_10x <- match(candidate_intersect_indices, barcode_10x_ambient)
    pvalue_indices_ed <- match(candidate_intersect_indices, barcode_umi_candidate_sorted)
    plot(candidate_loglk_sorted[pvalue_indices_ed], obs_loglk_10x[pvalue_indices_10x])


    #plot simulate prop
    ed_lk <- rowMeans(sim_loglk)
    lk_10x <- read.table('10x_loglk_mean.txt', header=FALSE)[,1]
    plot(ed_lk, lk_10x[1:length(ed_lk)], xlim=range(-1000:-12200), ylim=range(-1000:-12200))

    #plot  pvalues
    barcode_p_intersect <- intersect(barcode_umi_candidate_sorted, barcode_10x_ambient)
    pvalues_ed_indices <- match(barcode_p_intersect, barcode_umi_candidate_sorted)
    pvalues_10x_indices <- match(barcode_p_intersect, barcode_10x_ambient)
    plot(pvalues[pvalues_ed_indices], pvalue_10x[pvalues_10x_indices])

    # plot adj pvalues
    plot(p_adjust[pvalues_ed_indices], pvalue_10x_adj[pvalues_10x_indices])    
}


set.seed(0)

# read and count raw matrix
print("read and count raw matrix")
matrix_read <- read10xCounts(inputdir)
matrix_counts <- counts(matrix_read)
barcode_counts <- colSums(matrix_counts)

#print(barcode_counts)
print("raw_barcode")
#  find cells by top umi num barcode
# Number of additional barcodes to consider after the initial cell calling
print("find cells by top umi num barcode")
ORDMAG_NUM_BOOTSTRAP_SAMPLES = 100
quantile = 0.99
max_filtered_bcs = recovered_cells * 6 # not used
top_n_indices <- find_top_umi_barcode(barcode_counts, recovered_cells,
                                      ORDMAG_NUM_BOOTSTRAP_SAMPLES = 100, quantile = 0.99)
print(paste0("top_umi_barcode:", length(top_n_indices)))


# find ambient candidate barcode
print("find ambient candidate barcode")
min_umis_nonambient = 500 
min_umi_frac_of_median=0.01
N_CANDIDATE_BARCODES = 20000
ambient_candidate <- find_ambient_candidate(barcode_counts,  top_n_indices, min_umis_nonambient, N_CANDIDATE_BARCODES)
barcode_umi_candidate <- which(ambient_candidate)
barcode_umi_candidate_sorted <- barcode_umi_candidate[order(barcode_counts[barcode_umi_candidate])]
umi_sort_indices <- match(barcode_umi_candidate_sorted, barcode_umi_candidate)
print(paste0("ambient_candidate_barcode:", length(umi_sort_indices)))


cell_barcode <- top_n_indices
#cell_barcode
#print(length(cell_barcode))
output_10x_format_results(paste0(outputdir,'_'), matrix_read, matrix_counts, cell_barcode)

if (length(barcode_umi_candidate) == 0){
    # merge ambient and top_umi cell barcode
    cell_barcode <- top_n_indices
    print(paste0('number of all identifed cells: ', length(cell_barcode)))
    #output_10x_format_results
    output_10x_format_results(outputdir,  matrix_read, matrix_counts, cell_barcode)
    stop("no ambient candidate awailable and only top umi barcode was used")
}


# find_ambient_background barcodes
# Number of partitions (max number of barcodes to consider for ambient estimation) 
print("find_ambient_background barcodes")
N_PARTITIONS=90000
ambient <- find_ambient_background(matrix_counts, barcode_counts, N_PARTITIONS)

# compute prop for ambient
# select nonzero transcript

print("compute prop for ambient")
discard <- rowSums(matrix_counts) == 0
nonzero_transcript_matrix <- matrix_counts[!discard,,drop=FALSE]
ambient.cells <- nonzero_transcript_matrix[,ambient]
ambient.prof <- rowSums(ambient.cells)
print('start judge')
if (sum(ambient.prof)==0) {
    # merge ambient and top_umi cell barcode
    cell_barcode <- top_n_indices
    output_10x_format_results(outputdir,  matrix_read, matrix_counts, cell_barcode)    
    stop("no counts available to estimate the ambient profile and only top umi barcode was used")
}
ambient_prop <- compute_ambient_prop(nonzero_transcript_matrix)


# compute candidate prop(for ambient barcode)
print("compute candidate prop(for ambient barcode)")
nonzero_transcript_candidate_matrix <- nonzero_transcript_matrix[,ambient_candidate,drop=FALSE]
nonzero_transcript_candidate_t_matrix=t(as.matrix(nonzero_transcript_candidate_matrix))
ambient_prop_t = t(ambient_prop)
candidate_loglk = dmultinomial(nonzero_transcript_candidate_t_matrix, NULL, ambient_prop_t, TRUE)
candidate_loglk_sorted <- candidate_loglk[umi_sort_indices]


# compute simlate prop
print("compute simlate prop")
num_sims=10000
jump=1000
n_sample_feature_block=1000001
sim_loglk <- compute_simulate_pvalues(barcode_counts, ambient_prop,
                                  barcode_umi_candidate_sorted, num_sims, jump, n_sample_feature_block)

# compute pvalues
print("compute pvalues")
umi_per_bc_candidate <- barcode_counts[barcode_umi_candidate_sorted]
unique_n <- unique(umi_per_bc_candidate)
pvalues <- compute_ambient_candidate_pvalues(candidate_loglk_sorted,  sim_loglk, unique_n,  umi_per_bc_candidate)

# compute adj pvalues
print("compute adj pvalues")
p_adjust <- p.adjust(pvalues, method="BH")

# get cell barcode below p value cutoff
ambient_barcode <- barcode_umi_candidate_sorted[which(p_adjust < pvalue)]
print(paste0('number of cells identify from ambient: ', length(ambient_barcode)))

# merge ambient and top_umi cell barcode
cell_barcode <- c(ambient_barcode, top_n_indices)
head(cell_barcode)
print(length(cell_barcode))
print(paste0('number of all identifed cells: ', length(cell_barcode)))

#output_10x_format_results
output_10x_format_results(outputdir,  matrix_read, matrix_counts, cell_barcode)
