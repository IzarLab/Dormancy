# `nmf` kernel
library(tidyverse)
library(patchwork)
# library(viridis)
library(NMF)
library(Seurat)


indir <- '/home/ubuntu/dormancy_nmf'
outdir <- '/home/ubuntu/dormancy_nmf/nmf'

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive=T)
}

message_ <- function(...) {
    message(paste0(Sys.time(), ': ', paste(...)))
}

malig <- readRDS("dormancy.rds")
malig <- subset(malig, subset = major_celltypes_final=='Malignant')

sample_order <- malig@meta.data %>% 
                group_by(sample) %>% 
                summarize(n=n()) %>% 
                arrange(n) %>%
                mutate(index = seq(nrow(.)))
sample_order


produce_mat <- function(orig) {
    sub <- subset(malig, subset = sample == orig)
    message_(toString(dim(sub@assays$RNA$counts)))


    # Now extract the matrix
    all_genes <- rownames(sub@assays$RNA$counts)
    mito.genes <- grep(pattern = "^MT-", x = all_genes, value = TRUE)
    rbl.genes <- grep(pattern = "^RB-", x = all_genes, value = TRUE)
    rsl.genes <- grep(pattern = "^RS-", x = all_genes, value = TRUE)
    rpl.genes <- grep(pattern = "^RPL-", x = all_genes, value = TRUE)
    rbl.genes <- grep(pattern = "^RBL-", x = all_genes, value = TRUE)
    rps.genes <- grep(pattern = "^RPS-", x = all_genes, value = TRUE)
    rbs.genes <- grep(pattern = "^RBS-", x = all_genes, value = TRUE)
    rbl1.genes <- grep(pattern = "^RB", x = all_genes, value = TRUE)
    rsl1.genes <- grep(pattern = "^RS", x = all_genes, value = TRUE)
    rpl1.genes <- grep(pattern = "^RPL", x = all_genes, value = TRUE)
    rbl1.genes <- grep(pattern = "^RBL", x = all_genes, value = TRUE)
    rps1.genes <- grep(pattern = "^RPS", x = all_genes, value = TRUE)
    rbs1.genes <- grep(pattern = "^RBS", x = all_genes, value = TRUE)

    remove_list <- c(mito.genes, rbl.genes, rsl.genes, rpl.genes, rbl.genes,
                    rps.genes,rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,
                    rbl1.genes,rps1.genes,rbs1.genes)

    # Taking top 7k genes by mean expression
    good_genes <- rowSums(sub@assays$RNA$counts) %>% 
                sort(decreasing=T) %>% 
                as.data.frame() %>%
                rownames_to_column('gene') %>%
                filter(!(gene %in% remove_list)) %>%
                head(7000) %>%
                .$gene
    seu <- subset(sub, features = good_genes) 

    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    # seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    # seu <- ScaleData(seu, features = rownames(seu))

    mat <- as.matrix(seu@assays$RNA$data)
    mat <- mat[rowSums(mat)>0,]

    message_('Generated matrix for', orig)
    message_(toString(dim(seu@assays$RNA$counts)))
    return(mat)
}


#ct <- 0

# First run just to initialize all possible matrices
for (orig in c("SC78_9")) {
    #ct <- ct + 1
    #if (ct %% 4 != 1) next();  # TODO: EDIT BASED ON WHICH MACHINE IT'S RUNNING

    # outfile <- file.path(outdir, paste0(orig, '.pdf'))
    sample_outdir <- file.path(outdir, orig)
    if (!dir.exists(sample_outdir)) {
        dir.create(sample_outdir, recursive=T)
    }
    estim_outdir <- file.path(outdir, 'all_ranks')
    if (!dir.exists(estim_outdir)) {
        dir.create(estim_outdir, recursive=T)
    }
    

    ranks <- c(4, 5, 6, 7, 8, 9)
    outfiles <- lapply(ranks, 
                      function(x) file.path(sample_outdir, 
                                            paste0(orig, '_r', x, '.RDS')))
    names(outfiles) <- ranks
    rank_estim_rds <- file.path(estim_outdir, paste0(orig, '_all.RDS'))
    rank_estim_pdf <- file.path(estim_outdir, paste0(orig, '.pdf'))
    

    message_(orig)

    # For the individual rank calculationss
    donefile <- file.path(sample_outdir, 'done.checkpoint')
    if (!file.exists(donefile)) {
        mat <- produce_mat(orig)

        # Individual NMF ranks
        for (rank in ranks) {
            t_0 <- Sys.time()
            message_('Starting individual NMF for', orig, 'rank', rank)

            if (!file.exists(outfiles[[as.character(rank)]])) {
                # For 2700 cells, took between 1.5-2h at rank 7
                # Optimistic estimate of 12ish hours for full estimate
                # assuming the 10 runs are perfectly run in parallel
                res <- nmf(mat, 
                        rank,
                        method='brunet', 
                        nrun=10,
                        seed=123456,
                        rng=123456, 
                        .opt='vP60')
                saveRDS(res, outfiles[[as.character(rank)]])
            }

            t_1 <- Sys.time()
            message_('Done with individual NMF for', orig, 'rank', rank,
                     '. Elapsed time: ', t_1 - t_0)   
        }

    } else {
        message_('Already done with individual NMFs for', orig)
    }
    file.create(donefile)


    # For the rank estimation
    
    if (!file.exists(rank_estim_rds)) {
        message_('Now running rank estimation for', orig)
        
        mat <- produce_mat(orig)

        res <- nmf(mat, 
                    ranks,
                    method='brunet', 
                    nrun=10,
                    seed=123456,
                    rng=123456, 
                    .opt='vP60')
        saveRDS(res, rank_estim_rds)
        message_('Finished rank estimation for', orig)
    } else {
        message_('Already done with rank estimation for', orig)
    }
    
    message_('Fully done with', orig)
}
