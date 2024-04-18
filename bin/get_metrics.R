#! /usr/bin/env Rscript

suppressMessages(library(optparse))

option_list <- list(
    make_option(c("-p", "--peak"), type = "character", default = NULL,
                help="macs2 output peak results"),
    make_option(c("-f", "--fragment"), type = "character", default = NULL,
                help="fragment input files"),
    make_option(c("-c", "--cell"), type = "character", default = NULL,
                help="tsv file of the cells")
)

opt <- parse_args(OptionParser(option_list=option_list, description = "Summarize medianFragInPeaksPerCell, medianFragPerCell and total_reads from a fragment file"))

suppressMessages(library(tidyverse))
suppressMessages(library(Signac))
suppressMessages(library(rtracklayer))
suppressMessages(library(Seurat))

cells <- read_tsv(opt$cell, col_names = FALSE) %>% pull() %>% na.omit()
peaks <- rtracklayer::import(opt$peak)

frag <- CreateFragmentObject(
    path = opt$frag,
    cells = cells
)

counts <- FeatureMatrix(
    fragments = frag,
    features = peaks,
    cells = cells
)
atac.assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    min.features = 0,
    min.cells = 0,
    fragments = opt$frag
)

obj <- CreateSeuratObject(counts = atac.assay, assay = "ATAC")

## Calculate fragments in peaks for each cell
total_fragments <- CountFragments(opt$frag)
rownames(total_fragments) <- total_fragments$CB
obj$fragments <- total_fragments[colnames(obj), "frequency_count"]

medianFragsInCells <- median(obj$fragments)
medianFragsInPeaksInCells <- median(obj$nCount_ATAC)
totalFragsInCells <- sum(obj$fragments)
total_frags <- sum(total_fragments$frequency_count)
total_reads <- sum(total_fragments$reads_count)

cat(medianFragsInPeaksInCells, "\t", medianFragsInCells, "\t", totalFragsInCells, "\t", total_frags, "\t", total_reads, "\n", sep="")
