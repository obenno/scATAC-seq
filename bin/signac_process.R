#! /usr/bin/env Rscript

suppressMessages(library(optparse))


option_list <- list(
    make_option(c("-f", "--fragment"), type = "character", default = NULL,
                help="fragment input file"),
    make_option(c("-g", "--refGenome"), type = "character", default = "hg38",
                help="reference genome version, support hg38, hg19, mm9 and mm10"),
    make_option(c("-s", "--genomeSize"), type = "double", default = 2.7e9,
                help="define genome size, will overwrite the default ones"),
    make_option(c("-b", "--blacklistBED"), type = "character", default = NULL,
                help="blacklist bed file, will overwrite the default ones"),
    make_option(c("--gtf"), type = "character", default  = NULL,
                help="genome GTF file"),
    make_option(c("--minCell"), type = "integer", default = 10,
                help = "minimum cells"),
    make_option(c("--minFeature"), type = "integer", default = 200,
                help = "minimum featuers"),
    ##make_option(c("--TSS_enrichment_cutoff"), type = "integer", default = 3,
    ##            help = "TSS enrichment cutoff"),
    ##make_option(c("--nucleosome_signal_cutoff"),type = "integer", default = 4,
    ##            help = "Nucleosome signal cutoff"),
    make_option(c("--nCount_min"), type = "integer", default = 0,
                help = "Allowed min fragment counts in cells when selecting cells"),
    ##make_option(c("--nCount_max"), type = "integer", default = 30000,
    ##            help = "Allowed max fragment counts in cells when selecting cells"),
    ##make_option(c("--FRiP_cutoff"), type = "double", default = 0.15,
    ##            help = "Naive FRiP cutoff used when selecting cells"),
    ##make_option(c("--blacklist_fraction_cutoff"), type = "double", default = 0.05,
    ##            help = "Naive blacklist cutoff used when selecting cells"),
    make_option(c("--emptyDrops_fdr"), type = "double", default = 0.001,
                help = "FDR threshold for selecting cells different with ambientProfile of EmptyDrops algorithm"),
    make_option(c("--raw_cells_out"), type = "character", default = "cells.tsv",
                help = "Output tsv file containing barcodes of the cells detected"),
    make_option(c("--raw_meta_metrics"), type = "character", default = "raw_meta.tsv",
                help = "Output tsv file containing meta.data (FRiP, TSS.enrichment, fragmentCounts) of all barcodes"),
    make_option(c("--obj_out"), type = "character", default = "scATAC_obj.rds",
                help = "Output rds file of signac processed object"),
    make_option(c("-t", "--thread"), type = "integer", default = 4,
                help = "Threads to be used"),
    make_option(c("-m", "--memory"), type = "double", default = 10737418240,
                help = "Memory to be used, default 10GB")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

suppressMessages(library(tidyverse))
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(rtracklayer))
suppressMessages(library(future))
suppressMessages(library(reticulate))
suppressMessages(library(DropletUtils))

cmd <- commandArgs(trailingOnly = FALSE)
script_path <- cmd[str_detect(cmd, "^--file")] %>%
    str_remove("--file=") %>% dirname()
source_python(file.path(script_path, "otsu.py"))

message("Used threads: ", opt$thread)
message("Used memory: ", opt$memory)
set.seed(1234)
plan("multisession", workers = opt$thread)
options(future.globals.maxSize = opt$memory)
RhpcBLASctl::blas_set_num_threads(1) # https://github.com/satijalab/seurat/issues/3991

if(opt$refGenome == "hg38"){
    genomeSize <- 2.7e9
    blackListObj <- blacklist_hg38_unified
}else if(opt$refGenome == "hg19"){
    genomeSize <- 2.7e9
    blackListObj <- blacklist_hg19
}else if(opt$refGenome == "mm9"){
    genomeSize <- 1.87e9
    blackListObj <- blacklist_mm9
}else if(opt$refGenome == "mm10"){
    genomeSize <- 1.87e9
    blackListObj <- blacklist_mm10
}else if(opt$refGenome == "hg38-mm10"){
    genomeSize <- 2.7e9+1.87e9
    seqlevels(blacklist_hg38_unified) <- paste0("GRCh38_", seqlevels(blacklist_hg38_unified))
    seqlevels(blacklist_mm10) <- paste0("mm10___", seqlevels(blacklist_mm10))
    blackListObj <- c(blacklist_hg38_unified, blacklist_mm10)
}else{
    genomeSize <- 2.7e9
    blackListObj <- blacklist_hg38_unified
}

## overwrite genomeSize with input parameter
if(!is.null(opt$genomeSize)){
    genomeSize <- opt$genomeSize
}

if(!is.null(opt$blacklistBED)){
    blackListObj <- rtracklayer::import(opt$blacklistBED)
}

## Read gene annotation from GTF
annotations <- rtracklayer::import(opt$gtf)
annotations$gene_biotype <- annotations$gene_type
## Filter annotations by gene_biotype
annotations <- annotations[annotations$gene_biotype %in% c("protein_coding", "lincRNA", "rRNA", "processed_transcript")]


## Perform analysis
frag <- CreateFragmentObject(opt$fragment)
macs_peaks <- CallPeaks(frag, macs2.path = NULL,
                        outdir = "./",
                        effective.genome.size = genomeSize,
                        cleanup = FALSE)
counts <- FeatureMatrix(fragments = frag, features = macs_peaks,
                        sep= c(":", "-"))

## cell calling based fragments overlap peaks threshold selected by ostu algorithm
fragmentCutoff <- threshold_otsu(array(colSums(counts)))
message("otsu method threshold: ", fragmentCutoff)

if(opt$nCount_min > fragmentCutoff){
    fragmentCutoff <- opt$nCount_min
}
if(fragmentCutoff < 200){
    message("otsu method threshold too small (< 200), set nCount cutoff to 200")
    fragmentCutoff <- 200
}
otsu_cells <- colSums(counts)[colSums(counts) >= fragmentCutoff] %>% names

## Use emptyDrops method to adjust cells
ed_out<- emptyDrops(counts)
emptyDrops_cells <- colnames(counts)[ed_out$FDR <= 0.001]
selectedCells <- union(otsu_cells, emptyDrops_cells) %>% na.omit()
write_tsv(tibble(cells = selectedCells), file = opt$raw_cells_out, col_names = FALSE)

chrom_assay_raw <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = opt$fragment,
  min.cells = 0,
  min.features = 0
)

scATAC_raw <- CreateSeuratObject(
  counts = chrom_assay_raw,
  assay = "peaks"
)
Annotation(scATAC_raw) <- annotations
total_fragments_raw <- CountFragments(opt$fragment) %>%
    dplyr::filter(CB %in% colnames(scATAC_raw)) %>% dplyr::select(CB, frequency_count)

scATAC_raw@meta.data <- scATAC_raw@meta.data %>%
    rownames_to_column("CB") %>%
    left_join(total_fragments_raw, by="CB") %>%
    dplyr::rename("fragments" = frequency_count) %>%
    column_to_rownames("CB")

scATAC_raw <- FRiP(
  object = scATAC_raw,
  assay = 'peaks',
  total.fragments = 'fragments'
)
scATAC_raw <- TSSEnrichment(object = scATAC_raw, fast = FALSE)

d <- scATAC_raw[[c("FRiP", "fragments", "TSS.enrichment")]] %>%
    tibble::rownames_to_column("CB") %>%
    mutate(fragsOverlapPeaks = colSums(GetAssayData(scATAC_raw, assay ="peaks", slot="counts"))) %>%
    mutate(
        group = case_when(
            CB %in% selectedCells ~ "cell",
            TRUE ~ "background"
        )
    )

write_tsv(d, file = opt$raw_meta_metrics, col_names = TRUE)
## reconstruct count matrix using selected cells
counts <- FeatureMatrix(fragments = frag,
                        features = macs_peaks,
                        cells = selectedCells,
                        sep= c(":", "-"))
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = opt$fragment,
  min.cells = opt$minCell,
  min.features = opt$minFeature
)

obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

obj$orig.ident <- "orig"

##annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
##seqlevelsStyle(annotations) <- 'UCSC'
Annotation(obj) <- annotations

# compute nucleosome signal score per cell
obj <- NucleosomeSignal(object = obj)

# compute TSS enrichment score per cell
obj <- TSSEnrichment(object = obj, fast = FALSE)

obj$high.tss <- ifelse(obj$TSS.enrichment > 3, 'High', 'Low')

obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

## Calculate fragments in peaks for each cell
total_fragments <- CountFragments(opt$fragment)
rownames(total_fragments) <- total_fragments$CB
obj$fragments <- total_fragments[colnames(obj), "frequency_count"]

obj <- FRiP(
  object = obj,
  assay = 'peaks',
  total.fragments = 'fragments'
)

obj$blacklist_fraction <- FractionCountsInRegion(
  object = obj,
  assay = 'peaks',
  regions = blackListObj
)

## Non naive filters were applied, only select cells by ostu on fragments overlapping peaks
## This will allow users to further filter cells by their own

## Apply naive cutoffs
## Do not constrain TSS and NS
##message("cell selection cutoffs:")
##message("nCount_peaks > ", fragmentCutoff)
##message("nCount_peaks < ", opt$nCount_max)
##message("FRiP > ", opt$FRiP_cutoff)
##message("blacklist_fraction < ", opt$blacklist_fraction_cutoff)
##message("nucleosome_signal < ", opt$nucleosome_signal_cutoff)
##message("TSS.enrichment > ", opt$TSS_enrichment_cutoff)
##scATAC_obj <- subset(
##  x = obj,
##  subset = nCount_peaks > fragmentCutoff ##&
##    nCount_peaks < opt$nCount_max &
##    FRiP > opt$FRiP_cutoff &
##    blacklist_fraction < opt$blacklist_fraction_cutoff &
##    ##nucleosome_signal < opt$nucleosome_signal_cutoff &
##    ##TSS.enrichment > opt$TSS_enrichment_cutoff
##)

scATAC_obj <- obj
filtered_cells <- Cells(scATAC_obj)
message("Final cells: ", nrow(scATAC_obj[[]]))

scATAC_obj <- RunTFIDF(scATAC_obj)
scATAC_obj <- FindTopFeatures(scATAC_obj, min.cutoff = 'q0')
scATAC_obj <- RunSVD(scATAC_obj)

if(length(filtered_cells)>=150){
    p <- DepthCor(scATAC_obj)
    gb <- ggplot_build(p)
    pc1_cor <- gb$data[[1]] %>%
        filter(x==1) %>%
        pull(y)

    ## Exclude the first component if its correlation with sequencing depth is strong
    message("pc1_cor: ", pc1_cor)
    if(pc1_cor > 0.5 || pc1_cor < -0.5){
        reduction_dims <- 2:30
    }else{
        reduction_dims <- 1:30
    }
    scATAC_obj <- RunUMAP(object = scATAC_obj, reduction = 'lsi', dims = reduction_dims)
    scATAC_obj <- FindNeighbors(object = scATAC_obj, reduction = 'lsi', dims = reduction_dims)
    scATAC_obj <- FindClusters(object = scATAC_obj, verbose = FALSE, algorithm = 3)
}

saveRDS(scATAC_obj, file = opt$obj_out)
