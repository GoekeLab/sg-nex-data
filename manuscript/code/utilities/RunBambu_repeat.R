#!/usr/bin/env Rscript
# script for repeat analysis -run bambu with varied NDR thresholds ===============

rm(list = ls())


require(docopt)
'Usage:
RunBambu_repeat.R [-g <g> -t <t>]

Options:
-g index
-t type
]' -> doc

opts <- docopt(doc)
print(opts)
nnn <- as.integer(opts$g)
ttt <- as.integer(opts$t)

align_type <- c("best","noClose","primary")[ttt]


library(bambu)
rcSaveDir <- paste0("RunBambu12June/rc/",align_type)

rcfiles <- dir(rcSaveDir, pattern = "rds", full.names = TRUE)
bambuAnnotations <- readRDS("bambuAnnotations.rds")
genome.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
se <- bambu(reads = rcfiles,
            annotations = bambuAnnotations,
            genome = genome.file,
            ncore = 4,
            returnDistTable = TRUE,
            NDR = nnn/10,
            opt.em = list(degradationBias = FALSE),
            verbose=TRUE)
saveRDS(se, file = paste0("bambuOutput_24Aug2023_",align_type,"alignments_NDR",nnn/10,".rds"))
