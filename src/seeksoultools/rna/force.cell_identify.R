#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-i", "--inputdir"), type="character",
                help="path of input(full path of raw_feature_bc_matrix)"),
    make_option(c("-o", "--outputdir"), type="character",
                help="path of output, default=filtered_feature_bc_matrix under the sample parent dir of inputdir"),
    make_option(c("-f", "--force_cells"), type="integer",
                help="force cells")
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


print(paste0("inputdir:", opt$inputdir, "    outputdir:", opt$outputdir,
             "    force_cells: ", opt$forceCell))

suppressMessages({
    library(DropletUtils)
    library(Matrix)
    library(edgeR)
})

force_cells =  opt$force_cells
# read and count raw matrix
print("read and count raw matrix")
matrix_read <- read10xCounts(inputdir)
matrix_counts <- counts(matrix_read)
barcode_counts <- colSums(matrix_counts)


nonzero_index <- which(barcode_counts > 0)
bc_counts_nonzero <- barcode_counts[nonzero_index]
bc_counts_length <- length(bc_counts_nonzero)
head(nonzero_index)
need_barcode_counts <- nonzero_index[order(bc_counts_nonzero)[(bc_counts_length - force_cells + 1):bc_counts_length]]
#need_barcode_counts <- nonzero_index[order(bc_counts_nonzero)[1:force_cells]]
#o_barcode_counts=order(barcode_counts,decreasing=T)
#need_barcode_counts <- o_barcode_counts[1:force_cells]
#need_barcode_counts
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


output_10x_format_results(outputdir,  matrix_read, matrix_counts, need_barcode_counts)
