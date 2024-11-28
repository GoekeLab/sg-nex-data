calculateFragmentData <- function(sePathList, geneLengths.max, widthThreshold=0.25, short_read = FALSE){
    geneLengthsTbl <- tibble(id=names(geneLengths.max))
    countList <- list()
    countUnsplicedList = list()
    countFragmentList=list()
    meanCoverageList=list()
    names(sePathList) <- basename(sePathList)
    for(sePathName in names(sePathList)){
        seDist1 <- readRDS(file=sePathList[sePathName])
        if(short_read){
            rcCountByGene <- tapply(assay(seDist1)[,1]/(geneLengths.max[rowData(seDist1)$GENEID]/1000), rowData(seDist1)$GENEID, sum)
        }else{
            rcCountByGene <- tapply(assay(seDist1)[,1], rowData(seDist1)$GENEID, sum)
        }
        
        countList[[sePathName]] <- left_join(geneLengthsTbl,
                                             tibble(id=names(rcCountByGene),counts=rcCountByGene))[,2] ## total read class count for the gene
        setUnspliced = (rowData(seDist1)$confidenceType=='unsplicedWithin')
        rcCountByGeneUnspliced <- tapply(assay(seDist1)[setUnspliced,1], rowData(seDist1)$GENEID[setUnspliced], sum) 
        countUnsplicedList[[sePathName]] <-left_join(geneLengthsTbl,tibble(id=names(rcCountByGeneUnspliced),counts.unspliced=rcCountByGeneUnspliced))[,2]
        relWidth <- (sum(width(rowRanges(seDist1)))/geneLengths.max[rowData(seDist1)$GENEID]) # relative width defined as total read class width/gene width
        setWidth <- (relWidth<widthThreshold) # threshold set as 0.25 
        rcCountByRcWidth <- tapply(assay(seDist1)[setWidth,1], rowData(seDist1)$GENEID[setWidth], sum)
        countFragmentList[[sePathName]] <-left_join(geneLengthsTbl,tibble(id=names(rcCountByRcWidth),counts.shortFragment=rcCountByRcWidth))[,2]
        rcMeanCountByRcWidth <- tapply(relWidth*assay(seDist1)[,1], rowData(seDist1)$GENEID, sum)
        rcMeanCountByRcWidth <- tapply(relWidth*assay(seDist1)[,1], rowData(seDist1)$GENEID, sum)
        meanCoverageList[[sePathName]] <-left_join(geneLengthsTbl,tibble(id=names(rcMeanCountByRcWidth),meanCoverage=rcMeanCountByRcWidth))[,2]
    }
    countMatrix <- do.call(cbind, countList)
    countUnsplicedMatrix <- do.call(cbind, countUnsplicedList)
    countFragmentMatrix <- do.call(cbind, countFragmentList)
    meanCoverageMatrix <- do.call(cbind, meanCoverageList)
    names(countMatrix) <- names(sePathList)
    names(countUnsplicedMatrix) <- names(sePathList)
    names(countFragmentMatrix) <- names(sePathList)
    names(meanCoverageMatrix) <- names(sePathList)
    return(list(countMatrix,countUnsplicedMatrix,countFragmentMatrix,meanCoverageMatrix))
}



process_heatmap_data <- function(dt, methodNames = c("bambu_lr"), allCellLine = TRUE, gene= TRUE,
    genevec,ensemblAnnotations.transcripts,samples, runnamevec){
    cellLines <- c("A549","K562","HepG2","Hct116","MCF7","H9","HEYA8")
    protocols <- c("cDNA","directcDNA","directRNA","Illumina")
    pro_types <- c("protein_coding","antisense_RNA","lincRNA","non_coding","macro_lncRNA")
    rl_data <- dt[grepl("^ENSG",gene_name)&(gene_name %in% genevec)&(method %in% methodNames)]
    rl_data <- unique(ensemblAnnotations.transcripts[,.(gene_name, gene_biotype)])[rl_data, on = c("gene_name")]
    rl_data[, protocol_general := gsub("RandomPrimer", "", protocol_general)]
    if(allCellLine){
        filtered_rl_data <- rl_data[runname %in% runnamevec&(gene_biotype %in% pro_types)] 
    }else{
        filtered_rl_data <- rl_data[runname %in% runnamevec&(cellLine %in% cellLines)&(gene_biotype %in% pro_types)]
    }
    if(!gene){
        plotdata <- dcast(filtered_rl_data, tx_name ~ runname, value.var = 'normEst') 
    }else{
        plotdata <- dcast(filtered_rl_data, gene_name ~ runname, value.var = 'normEst')
    } 
    plotdata[is.na(plotdata)] <- 0
    tmp <- plotdata[,-1,with=FALSE]
    tmp <- log2(tmp+1)
    tmp[is.na(tmp)] <- 0
    runInfo <- samples[match(colnames(tmp),runname)]
    return(list(runInfo, tmp))
}





complexHeatmap_plot <- function(countMatrix,runInfo, number_of_genes = 1000){
    sdvec <- apply(countMatrix,1,sd)
    corMatrix <- cor(countMatrix[rank(-sdvec)<=number_of_genes,],method = 'spearman')
    col_fun = colorRamp2(seq(0.7,1,length.out = 8), brewer.pal(8,"BuGn"))
    colCellLines <- c(brewer.pal(8,"Dark2"),adjustcolor(brewer.pal(8,"Dark2"),alpha = 0.5), rev(brewer.pal(9,"Paired")))[seq_along(unique(runInfo$cellLine))]
    colProtocol <-  adjustcolor(brewer.pal(8,"Dark2")[1:4],0.7)
    colCancer <- brewer.pal(12,'Paired')[seq_along(unique(runInfo$cancer_type))]
    names(colCancer) <- unique(runInfo$cancer_type)
    names(colCellLines) <- unique(runInfo$cellLine)
    names(colProtocol) <- unique(runInfo$protocol_type)
    cellLine_anno = columnAnnotation(
        cellLine = as.factor(runInfo$cellLine),
        protocol = as.factor(runInfo$protocol_type),
        cancer_type = as.factor(runInfo$cancer_type),
        col = list(cellLine = colCellLines,
                   protocol = colProtocol,
                   cancer_type = colCancer),
        annotation_name_side = "left")
    colnames(corMatrix) <- NULL
    rownames(corMatrix) <- NULL
    p <- Heatmap(corMatrix, name = "Cor", col = col_fun, 
                 cluster_rows = TRUE,
                 cluster_columns =  TRUE,
                 top_annotation = cellLine_anno)
    return(p)
}

library(ggfortify)
plot_pca <- function(countMatrix, runInfo){
    sdvec <- apply(countMatrix,1,sd)
    corMatrix <- cor(countMatrix[rank(-sdvec)<=1000,],method = 'spearman')
    
    cor_pca <- prcomp(corMatrix)
    p <- ggplot2::autoplot(cor_pca, data = runInfo, colour = "cellLine", shape = "protocol_type") + theme_classic()
    print(p)
}

plot_heatmap <- function(dt, methodNames,gene, data_type,genevec,ensemblAnnotations.transcripts , samples, number_of_genes, runnamevec, plot_type = "PCA"){
    dataList <- process_heatmap_data(dt, methodNames,allCellLine = TRUE, gene = gene, genevec,ensemblAnnotations.transcripts , samples, runnamevec)
    runInfo <- dataList[[1]]
    data <- dataList[[2]]
    protocol_general <- runInfo$protocol_type
    cellLine <- runInfo$cellLine
    
    
    if(data_type == "all"){
        final_data <- data
        
    }else if(data_type == "cellline_batch_removal"){
        design0 <- model.matrix(~protocol_general)
        final_data <- limma::removeBatchEffect(data,batch = cellLine,design = design0) # ,
    }else if(data_type == "protocol_batch_removal"){
        design0 <- model.matrix(~cellLine)
        final_data <- limma::removeBatchEffect(data,batch = protocol_general, design = design0)
    }
    if(tolower(plot_type) == "pca"){
        return(plot_pca(final_data, runInfo))
    }else{
         return(complexHeatmap_plot(final_data, runInfo, number_of_genes))
    }
}

process_replicate <- function(dt, methodNames, gene = TRUE, 
samples, ensemblAnnotations.transcripts, genevec, runnamevec, 
majorMinor = TRUE, scatterPlot = TRUE, complexity = TRUE, 
expressionLevel = TRUE, metric_type_id = 1,bpParameters){ #reproducibility_check = TRUE, 
   
    rl_data <- dt[grepl("^ENSG",gene_name)&(gene_name %in% genevec)&(method %in% methodNames)] # already filtered 
    rl_data <- unique(ensemblAnnotations.transcripts[,.(gene_name, gene_biotype)])[rl_data, on = c("gene_name")]
    rl_data[, protocol_general := gsub("RandomPrimer", "", protocol_general)]
    filtered_rl_data <- rl_data[runname %in% runnamevec] 
    cellLines <- c('Hct116','HepG2','K562','A549','MCF7','H9','HEYA8')
    
    source(paste0('gene_cluster_code.R'))
    filtered_rl_data[, gene_cluster:=ifelse(gene_biotype %in% tr_gene_list,'TR gene',
                                    ifelse(gene_biotype %in% long_noncoding_rna_list, 'lncRNA',
                                           ifelse(gene_biotype %in% noncoding_rna_list, 'ncRNA',
                                                  ifelse(gene_biotype %in% pseudogene_list,'Pseudogene', ifelse(gene_biotype %in% ig_gene_list, 'IG gene',gene_biotype)))))]
    filtered_rl_data[, gene_cluster := ifelse(gene_cluster %in% c("IG gene","TR gene", "Mt_tRNA","Mt_rRNA","ribozyme"),"others", gene_cluster)]
    
    filtered_rl_data[, agg_gene_cluster := ifelse(gene_cluster == "processed_transcript", "lncRNA",
                                                 ifelse(gene_cluster %in% c("ncRNA","Pseudogene"), "others", gene_cluster))]

    protocolVec <- gsub("RandomPrimer","",unique(filtered_rl_data$protocol_general))
    cellLineList <- c(as.list(cellLines), list(cellLines))
    geneClusterList <- unique(filtered_rl_data$agg_gene_cluster)
    protein_coding_id <- grep("protein",geneClusterList)
    geneClusterList <- c(as.list(geneClusterList), list(geneClusterList))
   
        combMat <- combn(1:4,2)
   
    vv <- CJ(p = seq_len(ncol(combMat)), t = seq_along(cellLineList), g = seq_along(geneClusterList))
    if(majorMinor|complexity){
        vvIds <-  which(vv$g %in% c(protein_coding_id,length(geneClusterList))) # when major minor is considered, use proteining coding genes only 
    }else{
        vvIds <- seq_len(nrow(vv))
    }
    
    # filtering by at least highly expressed in either one, cause if not,
    # there might be a lot transcripts expression noise for lowly expressed transcripts,
    # and this will limit the finding of relationship
    if(scatterPlot){
        source("utility_function.R")
        np <- bplapply(vvIds[which(vv[vvIds]$t != 8)],pairwise_scatterplot_function , vv = vv, samples = samples,
                                            cellLineList = cellLineList, protocolVec = protocolVec, geneClusterList = geneClusterList, 
                                            combMat = combMat, filtered_rl_data = filtered_rl_data, gene = gene,
                                            majorMinor = majorMinor, 
                       complexity = complexity, expressionLevel = expressionLevel,
                       BPPARAM=bpParameters)
    }else{
        source("utility_function.R")
        mat_cor <- do.call("rbind",bplapply(vvIds ,pairwise_function , vv = vv, samples = samples,
                                            cellLineList = cellLineList, protocolVec = protocolVec, geneClusterList = geneClusterList, 
                                            combMat = combMat, filtered_rl_data = filtered_rl_data, gene = gene,
                                            majorMinor = majorMinor, 
                                            complexity = complexity, expressionLevel = expressionLevel,
                                            expression_t = 2,metric_type = c("cor","mae","mard","mard_mod","rmse")[metric_type_id],
                                            BPPARAM=bpParameters))
        return(mat_cor)
    }
       
    
    
    
    
}


pairwise_function <- function(v, vv, cellLineList, protocolVec, geneClusterList, combMat,
    filtered_rl_data, gene,  samples, majorMinor, complexity, expressionLevel, expression_t, metric_type){
    # cause the set of transcripts different from salmon and bambu, only use ENSG and ENST transcripts
    p <- vv[v]$p
    t <- vv[v]$t
    g <- vv[v]$g
    cellLineV <- cellLineList[[t]]
    protocolV <- protocolVec[combMat[,p]]
    geneCluster <- geneClusterList[[g]]
    print(paste(p,t,v))
    print(paste(cellLineV,protocolV,geneCluster))
    tmp <- filtered_rl_data[(cellLine %in% cellLineV)&(protocol_general %in% protocolV)&(agg_gene_cluster %in% geneCluster)]
    if(gene){
        tmp_wide <- dcast(tmp, gene_name ~ runname, value.var = "normEst")
    }else{
        tmp_wide <- dcast(tmp, tx_name + gene_name + ntx ~ runname, value.var = "normEst")
        if(majorMinor){
            tmp_wide <- unique(dominant_typeData[cellLine %in% cellLineV,
                .(tx_name, gene_name,majorBoth, majorEither, majorEitherOnly, majorSecBoth, majorLongRead, majorShortRead
                )])[tmp_wide, on = "tx_name"]
        }
        
    }
    tmp_wide[is.na(tmp_wide)] <- 0

    nameMat <- CJ(v1 = colnames(tmp_wide)[grep(paste0("_",protocolV[1]),colnames(tmp_wide))],
                  v2 = colnames(tmp_wide)[grep(paste0("_",protocolV[2]),colnames(tmp_wide))])
    # remove the one with the same name
    nameMat <- nameMat[v1 != v2]
    # remove the duplciated pair
    
    nameMat[, rep_status := (samples[which(samples$runname == v1)]$bioRep ==  samples[which(samples$runname == v2)]$bioRep), by = list(v1,v2)]
    nameMat[, rep := samples[which(samples$runname == v1)]$bioRep, by = v1]

  
    if(gene|(t>7)){
        get_corValues <- get("get_corValues1")
    }else{
        get_corValues <- get("get_corValues2")
    }
    if(t>7){ # remove exactly same pair
        nameMat[, cellline_status := (samples[runname == v1]$cellLine ==  samples[runname == v2]$cellLine), by = list(v1,v2)]
        nameMat <- nameMat[which(!cellline_status)]
        nameMat[, cellline_status := NULL]
    }
    
     temp_gene_cluster <- geneCluster
    if(length(geneCluster)>1) temp_gene_cluster <- "all"
   
    if(metric_type == "cor"){
        corValues <- NULL
        corValues <- get_corValues(tmp_wide, corValues, nameMat, typeName = "all", cellLineV, temp_gene_cluster, protocolV, expression_t)
        if(majorMinor){
            corValues <- get_corValues(tmp_wide[majorBoth == TRUE], corValues, nameMat, typeName = "majorBoth",cellLineV, temp_gene_cluster, protocolV, expression_t)
            corValues <- get_corValues(tmp_wide[majorEither == TRUE], corValues, nameMat, typeName = "majorEither", cellLineV, temp_gene_cluster, protocolV,expression_t)
            corValues <- get_corValues(tmp_wide[majorEitherOnly == TRUE], corValues, nameMat, typeName = "majorEitherOnly",cellLineV,temp_gene_cluster, protocolV,expression_t)
            corValues <- get_corValues(tmp_wide[majorSecBoth == TRUE], corValues, nameMat, typeName = "majorSecBoth", cellLineV, temp_gene_cluster, protocolV,expression_t)
            corValues <- get_corValues(tmp_wide[majorEitherOnly&majorLongRead], corValues, nameMat, typeName = "majorLongReadOnly",cellLineV,temp_gene_cluster, protocolV,expression_t)
            corValues <- get_corValues(tmp_wide[majorEitherOnly&majorShortRead], corValues, nameMat, typeName = "majorShortReadOnly",cellLineV,temp_gene_cluster, protocolV,expression_t)
            corValues <- get_corValues(tmp_wide[!(majorEitherOnly|majorEither|majorEitherOnly|majorSecBoth)], corValues, nameMat, typeName = "others", cellLineV, temp_gene_cluster, protocolV,expression_t)
            
        }
    
        if(complexity){
            if(majorMinor){
                corValues <- get_corValues(tmp_wide[ntx <=3&(majorBoth == TRUE)], corValues, nameMat, typeName = "<=3",cellLineV, temp_gene_cluster, protocolV, expression_t)
                corValues <- get_corValues(tmp_wide[ntx >3 &(ntx<=9)&(majorBoth == TRUE)], corValues, nameMat, typeName = "(3,9]",cellLineV, temp_gene_cluster, protocolV, expression_t)
                corValues <- get_corValues(tmp_wide[ntx >9 &(ntx<=15)&(majorBoth == TRUE)], corValues, nameMat, typeName = "(9,15]",cellLineV, temp_gene_cluster, protocolV, expression_t)
                corValues <- get_corValues(tmp_wide[ntx >15&(majorBoth == TRUE)], corValues, nameMat, typeName = "(15,193]",cellLineV, temp_gene_cluster, protocolV, expression_t)
                corValues <- get_corValues(tmp_wide[majorBoth == TRUE], corValues, nameMat, typeName = "all", cellLineV, temp_gene_cluster, protocolV, expression_t)
            }else{
                corValues <- get_corValues(tmp_wide[ntx <=3], corValues, nameMat, typeName = "<=3",cellLineV, temp_gene_cluster, protocolV, expression_t)
                corValues <- get_corValues(tmp_wide[ntx >3 &(ntx<=9)], corValues, nameMat, typeName = "(3,9]",cellLineV, temp_gene_cluster, protocolV, expression_t)
                corValues <- get_corValues(tmp_wide[ntx >9 &(ntx<=15)], corValues, nameMat, typeName = "(9,15]",cellLineV, temp_gene_cluster, protocolV, expression_t)
                corValues <- get_corValues(tmp_wide[ntx >15], corValues, nameMat, typeName = "(15,193]",cellLineV, temp_gene_cluster, protocolV, expression_t)
                
            }
        }
        return(corValues)
    }
    if(metric_type == "mae"){
        maeValues <- NULL
        if(majorMinor){
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorBoth == TRUE],"majorBoth", maeValues, expressionLevel, expression_t)
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEither == TRUE],"majorEither", maeValues, expressionLevel, expression_t)
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly == TRUE],"majorEitherOnly", maeValues, expressionLevel, expression_t)
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorSecBoth == TRUE], "majorSecBoth", maeValues, expressionLevel, expression_t)
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorLongRead], "majorLongReadOnly", maeValues, expressionLevel, expression_t)
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorShortRead], "majorShortReadOnly", maeValues, expressionLevel, expression_t)
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[!(majorEitherOnly|majorEither|majorEitherOnly|majorSecBoth)], "others", maeValues, expressionLevel, expression_t)
           
        }
        if(complexity){
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx <=3],"<=3", maeValues, expressionLevel, expression_t)
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >3 &(ntx<=9)],"(3,9]", maeValues, expressionLevel, expression_t)
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >9 &(ntx<=15)],"(9,15]", maeValues, expressionLevel, expression_t)
            maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >15],"(15,193]", maeValues, expressionLevel, expression_t)
        }
        maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide,"all", maeValues, expressionLevel, expression_t)
        return(maeValues)
    }
    if(metric_type == "mard"){
        mardValues <- NULL
        if(majorMinor){
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorBoth == TRUE],"majorBoth", mardValues, expressionLevel, expression_t)
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEither == TRUE],"majorEither", mardValues, expressionLevel, expression_t)
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly == TRUE],"majorEitherOnly",mardValues, expressionLevel, expression_t)
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorSecBoth == TRUE], "majorSecBoth",mardValues, expressionLevel, expression_t)
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorLongRead],"majorLongReadOnly",mardValues, expressionLevel, expression_t)
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorShortRead],"majorShortReadOnly",mardValues, expressionLevel, expression_t)
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[!(majorEitherOnly|majorEither|majorEitherOnly|majorSecBoth)], "others", mardValues, expressionLevel, expression_t)
        }
        if(complexity){
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx <=3],"<=3", mardValues, expressionLevel, expression_t)
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >3 &(ntx<=9)],"(3,9]", mardValues, expressionLevel, expression_t)
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >9 &(ntx<=15)],"(9,15]", mardValues, expressionLevel, expression_t)
            mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >15],"(15,193]", mardValues, expressionLevel, expression_t)
            
        }
       mardValues <-calc_mard(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide,"all", mardValues, expressionLevel, expression_t)
        return(mardValues)
    }
    if(metric_type == "mard_mod"){
        mardModValues <- NULL
        if(majorMinor){
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorBoth == TRUE],"majorBoth",mardModValues, expressionLevel, expression_t)
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEither == TRUE],"majorEither", mardModValues, expressionLevel, expression_t)
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly == TRUE],"majorEitherOnly",mardModValues, expressionLevel, expression_t)
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorSecBoth == TRUE], "majorSecBoth",mardModValues, expressionLevel, expression_t)
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorLongRead],"majorLongReadOnly",mardModValues, expressionLevel, expression_t)
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorShortRead],"majorShortReadOnly",mardModValues, expressionLevel, expression_t)
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[!(majorEitherOnly|majorEither|majorEitherOnly|majorSecBoth)], "others", mardModValues, expressionLevel, expression_t)
        }
        if(complexity){
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx <=3],"<=3", mardModValues, expressionLevel, expression_t)
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >3 &(ntx<=9)],"(3,9]", mardModValues, expressionLevel, expression_t)
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >9 &(ntx<=15)],"(9,15]", mardModValues, expressionLevel, expression_t)
            mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >15],"(15,193]", mardModValues, expressionLevel, expression_t)
           
        }
        mardModValues <-calc_mard_mod(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide,"all", mardModValues, expressionLevel, expression_t)
        return(mardModValues)
    }
    
    if(metric_type == "rmse"){
        rmseValues <- NULL
        if(majorMinor){
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorBoth == TRUE],"majorBoth", rmseValues, expressionLevel, expression_t)
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEither == TRUE],"majorEither", rmseValues, expressionLevel, expression_t)
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly == TRUE],"majorEitherOnly",rmseValues, expressionLevel, expression_t)
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorSecBoth == TRUE], "majorSecBoth",rmseValues, expressionLevel, expression_t)
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorLongRead],"majorLongReadOnly",rmseValues, expressionLevel, expression_t)
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorShortRead],"majorShortReadOnly",rmseValues, expressionLevel, expression_t)
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[!(majorEitherOnly|majorEither|majorEitherOnly|majorSecBoth)], "others", rmseValues, expressionLevel, expression_t)
        }
        if(complexity){
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx <=3],"<=3", rmseValues, expressionLevel, expression_t)
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >3 &(ntx<=9)],"(3,9]", rmseValues, expressionLevel, expression_t)
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >9 &(ntx<=15)],"(9,15]", rmseValues, expressionLevel, expression_t)
            rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >15],"(15,193]", rmseValues, expressionLevel, expression_t)
            
        }
        rmseValues <-calc_rmse(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide,"all", rmseValues, expressionLevel, expression_t)
        return(rmseValues)
    }
    
}


## version 1 : by replicates
get_corValues1 <- function(dd, corValues, nameMat, typeName = "majorBoth", cellLineV, geneCluster, protocolV, expression_t){
    if(nrow(dd)==0) return(corValues)
    for( x in seq_len(nrow(nameMat))){
    corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
    setnames(corData, c(1:2), c("V1","V2"))
    expressed <- which(apply(corData>expression_t,1,any)) # expressed in both
    log2CorData <- log2(corData+1)
    
    r2 <- summary(lm(V2~V1, data = log2CorData))$r.squared
    
    print(protocolV)
    temp_corValues <- data.table(r_expressed = cor(log2CorData[expressed], method = "spearman")[1,2],
                                 r = cor(log2CorData, method = "spearman")[1,2],
                                 r2 = r2,
                                 ne = length(expressed),
                                 match_status = nameMat[x]$rep_status,
                                 common_type = typeName,
                                 cellLine = paste(cellLineV, collapse = "_"),
                                 agg_gene_cluster = geneCluster,
                                 protocol_comparison = paste(protocolV, collapse = " vs "),
                                 n = nrow(dd))
    corValues <- do.call("rbind", list(corValues, temp_corValues))
    }
    return(corValues)
}



## version 2 : take mean across replicates
get_corValues2 <- function(dd, corValues, nameMat, typeName = "majorBoth",cellLineV, geneCluster, protocolV, expression_t){
    corData <- cbind(apply(dd[,nameMat$v1,with = FALSE],1,mean, na.rm = TRUE),
                   apply(dd[,nameMat$v2,with = FALSE],1,mean, na.rm = TRUE))
    expressed <- which(apply(corData>expression_t,1,any))
    temp_corValues <- data.table(r = cor(log2(corData+1))[1,2],
                                 spr = cor(log2(corData+1), method = "spearman")[1,2],
                                 r_expressed = cor(log2(corData+1)[expressed,])[1,2],
                                 spr_expressed = cor(log2(corData+1)[expressed,], method = "spearman" )[1,2],
                                 common_type = typeName,
                                 cellLine = paste(cellLineV, collapse = "_"),
                                 agg_gene_cluster = geneCluster,
                                 protocol_comparison = paste(protocolV, collapse = " vs "))
   
        corValues <- do.call("rbind", list(corValues, temp_corValues))
    return(corValues)
}

pairwise_scatterplot_function <- function(v, vv, cellLineList, protocolVec, geneClusterList, combMat,
                              filtered_rl_data, gene,  samples, majorMinor, complexity, expressionLevel){
    # cause the set of transcripts different from salmon and bambu, only use ENSG and ENST transcripts
    p <- vv[v]$p
    t <- vv[v]$t
    g <- vv[v]$g
    cellLineV <- cellLineList[[t]]
    protocolV <- protocolVec[combMat[,p]]
    geneCluster <- geneClusterList[[g]]
    print(paste(p,v))
    
    tmp <- filtered_rl_data[(cellLine %in% cellLineV)&(protocol_general %in% protocolV)&(gene_cluster %in% geneCluster)]
    if(gene){
        tmp_wide <- dcast(tmp, gene_name ~ runname, value.var = "normEst")
    }else{
        tmp_wide <- dcast(tmp, tx_name + ntx + gene_name ~ runname, value.var = "normEst")
        if(majorMinor){
            
            tmp_wide <- unique(dominant_typeData[,.(tx_name, gene_name,majorBoth, majorEither, majorEitherOnly,
                                                   majorSecBoth,majorFirstSecBoth, majorLongRead, majorShortRead
            )])[tmp_wide, on = "tx_name"]
        }
    }
    tmp_wide[is.na(tmp_wide)] <- 0
    
    nameMat <- CJ(v1 = colnames(tmp_wide)[grep(paste0("_",protocolV[1]),colnames(tmp_wide))],
                  v2 = colnames(tmp_wide)[grep(paste0("_",protocolV[2]),colnames(tmp_wide))])
    nameMat[, rep_status := (samples[runname == v1]$bioRep ==  samples[runname == v2]$bioRep), by = list(v1,v2)]
    nameMat[, rep := samples[runname == v1]$bioRep, by = v1]
    
    if(t>7){ # remove exactly same pair
        nameMat[, cellline_status := (samples[runname == v1]$cellLine ==  samples[runname == v2]$cellLine), by = list(v1,v2)]
        nameMat <- nameMat[which(!cellline_status)]
    }
    
    if(majorMinor){
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorBoth == TRUE],"majorBoth")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEither == TRUE],"majorEither")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly == TRUE],"majorEitherOnly")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorSecBoth == TRUE],"majorSecBoth")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorFirstSecBoth == TRUE],"majorFirstSecBoth")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorLongRead],"majorLongReadOnly")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[majorEitherOnly&majorShortRead],"majorShortReadOnly")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[!(majorEitherOnly|majorEither|majorEitherOnly|majorSecBoth)],"others")
    }
    
    if(complexity){
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx <=3],"<=3")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >3 &(ntx<=9)],"(3,9]")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >9 &(ntx<=15)],"(9,15]")
        plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >15],"(15,193]")
    }
    plot_scatter(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide,"")
}


plot_scatter <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix){
    if(nrow(dd)==0) return(NULL)
    pNames <- paste0(c(cellLineList[[8]],"all")[t], "_", 
                     
                     c(geneClusterList[[length(geneClusterList)]],"all")[g], "_",
                     nameMat$rep_status)
    pNamesOverall <- paste0(c(cellLineList[[8]],"all")[t], "_", 
                            protocolV[1],
                            protocolV[2],
                            c(geneClusterList[[length(geneClusterList)]],"all")[g],post_fix)
    pList <- vector("list",length(pNames))
    for( x in seq_len(nrow(nameMat))){
        corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
        varnames <- colnames(corData)
        setnames(corData, c(1:2), c("V1","V2"))
        expressed <- which(apply(corData,1,sum)>0)
        
            p_scatter <- ggplot(data = corData, aes(x = log2(V1+1), y = log2(V2+1)))+
               
                xlab(varnames[1])+
                ylab(varnames[2])+
                geom_hex(aes(fill = stat(count)),
                         binwidth = 0.25) +
                scale_fill_gradient(name = 'count', low = "grey", high = "black")+
               
                stat_cor(aes(label = ..r.label..),method = "spearman", cor.coef.name = "Sp.R")+
                ggtitle(pNames[x])+
               
                theme_classic()
        pList[[x]] <- p_scatter
    }
    
    names(pList) <- pNames
    saveplot.dir <- "replicate_plot_complexity/" 
    if(!dir.exists(saveplot.dir)) dir.create(saveplot.dir)
    png(paste0(saveplot.dir,pNamesOverall,".png"), width = 10, height = ceiling(length(pList)/4)*2, units = "in", res = 300)
    do.call("grid.arrange", c(pList, ncol=4))
    dev.off()
}

# this is essentially the same as the one used in salmon: median(log2(estimated/truth))
calc_mae <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
    diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
    for( x in seq_len(nrow(nameMat))){
        corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
        varnames <- colnames(corData)
        setnames(corData, c(1:2), c("V1","V2"))
        expressed <- which(apply(corData,1,sum)>0)
        diffMat[expressed,x] = abs(log2(corData[expressed]$V1+1)-log2(corData[expressed]$V2+1))
        
    }
    temp_maeValues <- copy(nameMat)
    expressed <- apply(corData>expression_t,1,sum)>0
    if(expressionLevel){
        temp_maeValues[,mae := apply(diffMat[expressed],2,function(x) mean(x, na.rm = TRUE))]
    }else{
        temp_maeValues[,mae := apply(diffMat,2,function(x) mean(x, na.rm = TRUE))]
    }
    setnames(temp_maeValues, "rep_status","match_status")
    
    temp_maeValues[, `:=`(
        cellLine = c(cellLineList[[8]],"all")[t],
        gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
         common_type = post_fix,
        protocol_comparison = paste(protocolV, collapse = " vs "))
    ]
    
    
    maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
    return(maeValues)
}

calc_rmse <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
    diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
    for( x in seq_len(nrow(nameMat))){
        corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
        varnames <- colnames(corData)
        setnames(corData, c(1:2), c("V1","V2"))
        diffMat[(corData$V1+corData$V2)>0,x] = (log2(corData[(corData$V1+corData$V2)>0]$V1+1)-log2(corData[(corData$V1+corData$V2)>0]$V2+1))^2
        
    }
    temp_maeValues <- copy(nameMat)
    expressed <- apply(corData>expression_t,1,sum)>0
    if(expressionLevel){
        temp_maeValues[,mae := apply(diffMat[expressed],2,function(x) sqrt(mean(x, na.rm = TRUE)))]
    }else{
        temp_maeValues[,mae := apply(diffMat,2,function(x) sqrt(mean(x, na.rm = TRUE)))]
    }
    setnames(temp_maeValues, "rep_status","match_status")
    
    temp_maeValues[, `:=`(
        cellLine = c(cellLineList[[8]],"all")[t],
        gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
        common_type = post_fix,
        protocol_comparison = paste(protocolV, collapse = " vs "))
    ]
    
    maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
    return(maeValues)
}


# mean absolute relative difference used by both salmon and kallisto: 2*abs(estimated-truth)/(estimated+truth)
calc_mard <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
    diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
    for( x in seq_len(nrow(nameMat))){
        corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
        varnames <- colnames(corData)
        setnames(corData, c(1:2), c("V1","V2"))
        corData[, V1 := log2(V1+1)]
        corData[, V2 := log2(V2+1)]
        expressed <- which(apply(corData[,c("V1","V2"),with = TRUE],1,sum)>0)
        diffMat[expressed,x] = 2*abs(corData$V1-corData$V2)[expressed]/(corData$V1+corData$V2)[expressed]
        
    }
    
    temp_maeValues <- copy(nameMat)
    expressed <- apply(corData>expression_t,1,sum)>0
    if(expressionLevel){
        temp_maeValues[,mae := apply(diffMat[expressed],2,function(x) mean(x, na.rm = TRUE))]
    }else{
        temp_maeValues[,mae := apply(diffMat,2,function(x) mean(x, na.rm = TRUE))]
    }
    setnames(temp_maeValues, "rep_status","match_status")
    
    temp_maeValues[, `:=`(
        cellLine = c(cellLineList[[8]],"all")[t],
        gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
        common_type = post_fix,
        protocol_comparison = paste(protocolV, collapse = " vs "))
    ]
    
    
    maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
    return(maeValues)
}


calc_mard_mod <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
    diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
    for( x in seq_len(nrow(nameMat))){
        corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
        varnames <- colnames(corData)
        setnames(corData, c(1:2), c("V1","V2"))
        corData[, V1 := log2(V1+1)]
        corData[, V2 := log2(V2+1)]
        expressed <- which(apply(corData[,c("V1","V2"),with = TRUE],1,sum)>0)
        diffMat[expressed,x] = 2*(corData$V1-corData$V2)[expressed]/(corData$V1+corData$V2)[expressed]
        
    }
    
    temp_maeValues <- copy(nameMat)
    expressed <- apply(corData>expression_t,1,sum)>0
    if(expressionLevel){
        temp_maeValues[,mae := apply(diffMat[expressed],2,function(x) mean(x, na.rm = TRUE))]
    }else{
        temp_maeValues[,mae := apply(diffMat,2,function(x) mean(x, na.rm = TRUE))]
    }
    setnames(temp_maeValues, "rep_status","match_status")
    
    temp_maeValues[, `:=`(
        cellLine = c(cellLineList[[8]],"all")[t],
        gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
        common_type = post_fix,
        protocol_comparison = paste(protocolV, collapse = " vs "))
    ]
    maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
    return(maeValues)
}




calc_mard_ave <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
    diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
    for( x in seq_len(nrow(nameMat))){
        corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
        varnames <- colnames(corData)
        setnames(corData, c(1:2), c("V1","V2"))
        diffMat[(corData$V1+corData$V2)>0,x] = abs(corData[(corData$V1+corData$V2)>0]$V1-corData[(corData$V1+corData$V2)>0]$V2)/(corData[(corData$V1+corData$V2)>0]$V1+corData[(corData$V1+corData$V2)>0]$V2)
    }
    if("tx_name" %in% colnames(dd)){
        temp_maeValues <- dd[,.(tx_name, gene_name, ntx)]
    }else{
        temp_maeValues <- dd[,.(gene_name)]
    }
    if(expressionLevel) temp_maeValues[, expressed := apply(corData>expression_t,1,sum)>0]
    temp_maeValues[,mae := apply(diffMat,1,mean,na.rm = TRUE)]
    temp_maeValues[,`:=`(match_mae = NA, non_match_mae = NA)]
    if(length(which(nameMat$rep_status))>1){
        temp_maeValues[,match_mae := apply(diffMat[, which(nameMat$rep_status)],1,mean,na.rm = TRUE)]
    }else if(length(which(nameMat$rep_status))==1){
        temp_maeValues[,match_mae := mean(diffMat[, which(nameMat$rep_status)],na.rm = TRUE)]
    }
    if(length(which(!nameMat$rep_status))>1){
        temp_maeValues[,non_match_mae := apply(diffMat[, which(!nameMat$rep_status)],1,mean,na.rm = TRUE)]
    }else if(length(which(!nameMat$rep_status))==1){
        temp_maeValues[,match_mae := mean(diffMat[, which(!nameMat$rep_status)],na.rm = TRUE)]
    }
    
    temp_maeValues[, `:=`(
        cellLine = c(cellLineList[[8]],"all")[t],
        gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
        common_type = post_fix,
        protocol_comparison = paste(protocolV, collapse = " vs "))
    ]
    
    maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
    return(maeValues)
}

# add isoform rank for gene and transcript expression
# it should be a combined dataset
# we fix it for all combinations, i.e., we identify major isoforms in 
# both long read and short read and apply it across all scenarios 
identifyMajorMinorIsoforms <- function(rl_data){#com_data
    
    
    rl_data_gene <- unique(rl_data[, list(tpm=sum(normEst)), by = list(cellLine, protocol_general, gene_name, runname, gene_biotype, method)])
    rl_data_gene_ave <- unique(rl_data_gene[, list(tpm = mean(tpm)), by = list(cellLine, protocol_general, gene_name, gene_biotype, method)])
    
    rl_data_gene_ave[, cellLine_protocol:=paste0(cellLine," ",protocol_general)]
    rl_data_gene_ave[, tpm_log:=log2(tpm+1)]
    
    rl_data_tx_ave <- unique(rl_data[, list(tpm = mean(normEst)), by = list(cellLine, short_read, gene_name, tx_name, method)])
    rl_data_tx_ave_new <- dcast(rl_data_tx_ave, tx_name+gene_name+cellLine ~ short_read, value.var = "tpm")
    cols <- colnames(rl_data_tx_ave_new)[4:5]
    rl_data_tx_ave_new[ , (cols) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols =cols]
    
    rl_data_tx_ave <- melt(rl_data_tx_ave_new, id.vars = colnames(rl_data_tx_ave_new)[1:3], measure.vars = cols)
    setnames(rl_data_tx_ave, old = c("variable","value"),new = c("short_read","tpm"))
    rl_data_tx_ave[, geneExpression := sum(tpm), by = list(gene_name, cellLine, short_read)]
    rl_data_tx_ave[, isoform_rank:=rank(-tpm, ties.method = "random"), 
        by = list(cellLine,short_read,gene_name)]
    rl_data_tx_ave[, geneExpressedInBoth := all(geneExpression>0), by = list(gene_name, cellLine)]
   
    rl_data_tx_ave[, majorBoth := all(isoform_rank==1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorEither := any(isoform_rank==1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorEitherOnly := (sum(isoform_rank==1)==1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorSecBoth := all(isoform_rank==2)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorFirstSecBoth := (all(isoform_rank==2)|all(isoform_rank == 1))&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorLongRead := any(short_read==FALSE&isoform_rank == 1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorShortRead := any(short_read==TRUE&isoform_rank == 1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    return(rl_data_tx_ave)
}

majorMinor_generic <- function(rl_data){
    rl_data_tx_ave <- unique(rl_data[, list(tpm = mean(normEst)), by = list(protocol_general,gene_name, tx_name)])
    rl_data_tx_ave_new <- dcast(rl_data_tx_ave, tx_name+gene_name ~ protocol_general, value.var = "tpm")
    cols <- colnames(rl_data_tx_ave_new)[3:4]
    rl_data_tx_ave_new[ , (cols) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols =cols]
    
    rl_data_tx_ave <- melt(rl_data_tx_ave_new,
                           id.vars = colnames(rl_data_tx_ave_new)[1:2], 
                           measure.vars = cols)
    setnames(rl_data_tx_ave, old = c("variable","value"),new = c("protocol_general","tpm"))
    rl_data_tx_ave[, geneExpression := sum(tpm), by = list(gene_name, protocol_general)]
    rl_data_tx_ave[, isoform_rank:=rank(-tpm, ties.method = "random"), 
                   by = list(protocol_general,gene_name)]
    rl_data_tx_ave[, geneExpressedInBoth := all(geneExpression>0), by = list(gene_name)]
   
    rl_data_tx_ave[, majorBoth := all(isoform_rank==1)&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    rl_data_tx_ave[, majorEither := any(isoform_rank==1)&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    rl_data_tx_ave[, majorEitherOnly := (sum(isoform_rank==1)==1)&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    rl_data_tx_ave[, majorSecBoth := all(isoform_rank==2)&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    rl_data_tx_ave[, majorFirstSecBoth := (all(isoform_rank==2)|all(isoform_rank == 1))&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    return(rl_data_tx_ave)
}


repeat_analysis_function <- function(seOutput, ervRanges, retrotransposonOnly = FALSE){

    com_data <- process_seOutput(seOutput)
    
    extendedAnnotationGRangesList <- rowRanges(seOutput)
    
    isoTEratio <- estimate_repeat_ratio(extendedAnnotationGRangesList,anno_exByTx,ervRanges, txLengths.tbldf, retrotransposonOnly = retrotransposonOnly)
    
    return(list(com_data, isoTEratio))
}


process_seOutput <- function(seOutput){
    tmp <- data.table(as.data.frame(rowData(seOutput)))
   
    fullLengthCounts <- as.data.table(assays(seOutput)$fullLengthCounts, keep.rownames = TRUE)
    totalCounts <- as.data.table(assays(seOutput)$counts, keep.rownames = TRUE)
    uniqueCounts <- as.data.table(assays(seOutput)$uniqueCounts, keep.rownames = TRUE)
    incompatibleCountsTable <- data.table(as.data.frame(metadata(seOutput)$incompatibleCounts))
    incompatibleDt <- data.table(runname = colnames(incompatibleCountsTable)[-1],
    incompatibleCounts = as.numeric(colSums(apply(incompatibleCountsTable[,-1, with = FALSE], c(1,2), as.numeric))))
    geneTxTable <- as.data.table(rowData(seOutput))
    setnames(fullLengthCounts, "rn","tx_name")
    setnames(uniqueCounts, "rn","tx_name")
    setnames(totalCounts, "rn","tx_name")
    setnames(geneTxTable, c("GENEID","TXNAME"),c("gene_name","tx_name"))
    
    fullLengthCounts_lr <- melt(fullLengthCounts, id.vars = "tx_name", measure.vars = colnames(fullLengthCounts)[-1])
    setnames(fullLengthCounts_lr, c("variable","value"),c("runname","fullLengthCounts"))
    
    uniqueCounts_lr <- melt(uniqueCounts, id.vars = "tx_name", measure.vars = colnames(uniqueCounts)[-1])
    setnames(uniqueCounts_lr, c("variable","value"),c("runname","uniqueCounts"))
    
    
    totalCounts_lr <- melt(totalCounts, id.vars = "tx_name", measure.vars = colnames(totalCounts)[-1])
    setnames(totalCounts_lr, c("variable","value"),c("runname","totalCounts"))
    
    genomeem_lr <- merge(fullLengthCounts_lr, totalCounts_lr, by = c("tx_name","runname"))
    genomeem_lr <- merge(uniqueCounts_lr, genomeem_lr, by = c("tx_name","runname"))
    
    genomeem_lr[, runname:=gsub("HCT116","Hct116",gsub("(.genome_alignment.sorted)|(_R1.sorted)|(_sorted)","",gsub("-pre|-Pre","",runname))), by = runname]
    genomeem_lr[runname == "GIS_Hct116_cDNA_Rep2_Run4", runname:="GIS_Hct116_cDNA_Rep2_Run5"]
     genomeem_lr[, `:=`(fullLengthCounts =sum(fullLengthCounts),
                       uniqueCounts =sum(uniqueCounts),
                       totalCounts = sum(totalCounts)), by = list(runname, tx_name)]
    genomeem_lr <- unique(genomeem_lr[,.(tx_name, fullLengthCounts,uniqueCounts, totalCounts, runname)])
    genomeem_lr[, ntotal:=sum(totalCounts), by = runname]
  
    genomeem_lr <- geneTxTable[genomeem_lr, on = "tx_name"]
    genomeem_lr[, method := "genomeem_lr"]
    
    
    merge.colnames <- c('tx_name','gene_name','fullLengthCounts','uniqueCounts','totalCounts','runname','method','ntotal')
    com_data <- genomeem_lr[, merge.colnames, with =FALSE]
 
    
    com_data[, runname := as.character(runname)]
    com_data[, protocol:=unlist(strsplit(runname, 
                                         '_'))[3],by= runname]
    
    com_data[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
    com_data[, protocol_method:=paste0(protocol_general,'.',method)]
    com_data[, short_read:=as.numeric(protocol_general == 'Illumina')]
    com_data[, runname_method:=paste0(runname,'.',method)]
    com_data[, cellLine := ifelse(grepl("uniAln$",runname), "uni_aln_filter", 
                                  ifelse(grepl("priAln$",runname), "pri_uni_aln_filter", "no_filter")), by = runname]
    return(com_data)
}


estimate_repeat_ratio <- function(extendedAnnotationGRangesList,anno_exByTx, ervRanges, txLengths.tbldf,retrotransposonOnly=FALSE){
    
    geneTxTable_extended <- as.data.table(mcols(extendedAnnotationGRangesList))
    setnames(geneTxTable_extended, c("GENEID","TXNAME"),c("gene_name","tx_name"))
    geneTxTable_extended[,nisoform:=length(unique(tx_name)), by = gene_name]
    
    geneTxTable <- txLengths.tbldf[,.(tx_name, gene_id, nisoform)]
    setnames(geneTxTable, 'gene_id', 'gene_name')
    geneTxTable_extended <- geneTxTable[geneTxTable_extended, on = c("gene_name","tx_name")]
    if(any(grepl("txClassDescription", colnames(geneTxTable_extended)))){
        geneTxTable_extended[, newTxClassAggregated:=ifelse(grepl("newFirstExon",txClassDescription)&(grepl("newLastExon",txClassDescription)),"newFirstLastExon",
                                                            ifelse(grepl("newFirstExon",txClassDescription), "newFirstExon",
                                                                   ifelse(grepl("newLastExon",txClassDescription), "newLastExon",
                                                                          ifelse(grepl("Junction|allNew",txClassDescription),"newJunction",
                                                                                 ifelse(grepl("unspliced",txClassDescription),"unspliced",
                                                                                        ifelse(grepl("newGene-",txClassDescription),"newGene",txClassDescription))))))]
        
        m <- table(grepl("unspliced",mcols(extendedAnnotationGRangesList)$txClassDescription))
        m  ## FALSE 203102 TRUE 119
        prop.table(m) # 99.94% are not unspliced about 0.06% of transcripts are removed 
        
        # in terms of novel transcript
        m <- table(grepl("unspliced",mcols(extendedAnnotationGRangesList[mcols(extendedAnnotationGRangesList)$txClassDescription != "annotation"])$txClassDescription))
        m ## 2451 119
        prop.table(m) 
        extendedAnnotationGrangesList <- extendedAnnotationGRangesList[!grepl("unspliced",mcols(extendedAnnotationGRangesList)$txClassDescription)]
        
    }
    # 95.4 and 4.6 this will remove about 4.6% of new transcripts
    
    exons_granges <- unlist(extendedAnnotationGRangesList)
    anno_exons <- unlist(anno_exByTx)
    ov <- findOverlaps(exons_granges, anno_exons, type = "any")
    
    exons_granges$anno_status <- "novel"
    exons_granges[queryHits(ov)]$anno_status <- "annotated"
    seqlevelsStyle(ervRanges) <- 'NCBI'
    ## reduce granges by repeat class 
    grl <- split(ervRanges, ervRanges$repClass)
    
    # reduce each element (independently reduce ranges for each gene)
    grl_redux <- reduce(grl) # element-wise, like lapply(grl, reduce) 
    
    # return single GRanges, with rownames derived from hgnc
    ervRangesReduceByRepClass <- unlist(grl_redux)
    ervRangesReduceByRepClass$repClass <- names(ervRangesReduceByRepClass)
    
    if(retrotransposonOnly){
        ervRangesReduceByRepClass <- ervRangesReduceByRepClass[grep("LINE|SINE|LTR",ervRangesReduceByRepClass$repClass)]
    }
    
    ovByRepClass <- compute_overlap(exons_granges, ervRangesReduceByRepClass,ignore.strand, by_rep_type = TRUE)
    
    isoTEratio_all_byrepclass <- unique(ovByRepClass[, list(repRatio_byisoform = sum(repRatio_all)), by = list(txId, rep_class)])
    isoTEratio_byrepclass <- unique(ovByRepClass[, list(repRatio = sum(repRatio)), by = list(txId, rep_class, anno_status)])
    isoTEratio_byrepclass <-isoTEratio_all_byrepclass[isoTEratio_byrepclass, on = c("txId","rep_class")]
   
    ### by all repeat together 
    reducedErvRanges <- reduce(ervRanges)
    ovOverall <- compute_overlap(exons_granges, reducedErvRanges,ignore.strand, by_rep_type = FALSE)
    isoTEratio_all <- unique(ovOverall[, list(repRatio_byisoform_all = sum(repRatio_all)), by = txId])
    isoTEratio <- unique(ovOverall[,list(repRatio_all = sum(repRatio)), by = list(txId,anno_status)])
    isoTEratio <-isoTEratio_all[isoTEratio, on = c("txId")]
    
    isoTEratio <- isoTEratio[isoTEratio_byrepclass, on = c("txId","anno_status")]
    
    setnames(isoTEratio, "txId","tx_name")
    isoTEratio <- geneTxTable_extended[isoTEratio, on = "tx_name"]
   
    return(isoTEratio)
}

gffcompare_function <- function(ref.gtf,query.gtf){
    system(paste0("gffcompare -TNRQS -e 15 -d 15 -o gffCompare ", 
                  query.gtf, " -r ", ref.gtf),ignore.stderr = TRUE)
    queryTx <- read.delim(paste0("./gffCompare.combined.gtf"),header=FALSE,comment.char='#')
    colnames(queryTx) <- c("seqname","source","type","start","end","score","strand","frame","attribute")
    queryTx <- data.table(queryTx)
    queryTx <- queryTx[type=='transcript']
    
    queryTx[, `:=`(tx_id = gsub('.*transcript_id (.*?);.*', '\\1',attribute),
                   TXNAME = gsub('.*oId (.*?);.*', '\\1',attribute),
                   CONTAIN_id = ifelse(grepl("contained_in",attribute),gsub('.*contained_in (.*?);.*', '\\1',attribute),NA),
                   REFTXNAME = ifelse(grepl("cmp_ref",attribute),gsub('.*cmp_ref (.*?);.*', '\\1',attribute),NA),
                   class_code = ifelse(grepl("class_code",attribute), gsub('.*class_code (.*?);.*', '\\1',attribute),NA))]
    
    queryTx[class_code == "=", NEWTXNAME := REFTXNAME]
    queryTx[class_code != "=", NEWTXNAME := paste0(REFTXNAME,"_",TXNAME,"_",seqname)]
    
    return(queryTx)
}

compute_overlap <- function(exons_granges, tmpErvRanges,ignore.strand, by_rep_type = FALSE){
    ov <- findOverlaps(exons_granges,tmpErvRanges, ignore.strand = ignore.strand)
    p <- Pairs(exons_granges, tmpErvRanges, hits=ov)
    hitIntersect <- pintersect(p, ignore.strand = ignore.strand)
    rm(p)
    gc()
    
    overlapWidth <- width(hitIntersect)
    rm(hitIntersect)
    gc()
    
    qHits <- queryHits(ov)
    sHits <- subjectHits(ov)
    rm(ov)
    gc()
    
    nnovel <- data.table(table(names(exons_granges[exons_granges$anno_status=="novel"])))
    nannotated <- data.table(table(names(exons_granges[exons_granges$anno_status=="annotated"])))
    exons_granges$nnovel <- nnovel[match(names(exons_granges),V1)]$N
    exons_granges$nannotated <- nannotated[match(names(exons_granges),V1)]$N  
    
    txDt <- data.table(txId = names(exons_granges), ## overlapping with repeats at least
                       exon_rank = exons_granges$exon_rank,
                       anno_status = exons_granges$anno_status,
                       nNovelExon = exons_granges$nnovel,
                       nAnnotatedExon = exons_granges$nannotated,
                       exon_width = width(exons_granges))
    txDt[, tx_width := sum(exon_width), by = txId]
    txDt[, tx_width_by_status := sum(exon_width), by = list(txId,anno_status)]
    
    if(by_rep_type){
        ov <- data.table(txId = names(exons_granges)[qHits], ## overlapping with repeats at least
                         exon_rank = exons_granges[qHits]$exon_rank,
                         rep_class = tmpErvRanges[sHits]$repClass,
                         rep_width = width(tmpErvRanges)[sHits],
                         ov_width = overlapWidth)
        
        txDt <- txDt[txId %in% unique(ov$txId)]
        repDt <- data.table(rep_class = unique(tmpErvRanges$repClass))
        tx_rep_dt <- optiRum::CJ.dt(txDt, repDt)
        ov <- ov[tx_rep_dt, on = c("txId","exon_rank","rep_class")]
        ov[is.na(ov_width), ov_width := 0]
    }else{
        ov <- data.table(txId = names(exons_granges)[qHits], ## overlapping with repeats at least
                         exon_rank = exons_granges[qHits]$exon_rank,
                         rep_width = width(tmpErvRanges)[sHits],
                         ov_width = overlapWidth)
        txDt <- txDt[txId %in% unique(ov$txId)]
        ov <- ov[txDt, on = c("txId","exon_rank")]
        ov[is.na(ov_width), ov_width := 0]
    }
    ov[, repRatio := ov_width/tx_width_by_status] #exon_width
    ov[, repRatio_all := ov_width/tx_width]
    
    ov[is.na(nNovelExon), nNovelExon := 0]
    ov[is.na(nAnnotatedExon), nAnnotatedExon := 0]
    return(ov)
}

process_salmonOutput <- function(filePaths, txLengths){
    tx2gene <- txLengths[,c(2,3)]
    x <- 1
    
    salmon_sr <- do.call('rbind',lapply(filePaths,function(filePath){
        print(filePath)
        txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
        names(txi)
        
        short_read <- data.table(tx_name = rownames(txi$abundance),
                                 abundance = txi$abundance[,1],
                                 counts = txi$counts[,1],
                                 length = txi$length[,1],
                                 countsFromAbundance = txi$countsFromAbundance)
        short_read[, runname:=basename(gsub("\\/transcripts_quant","",dirname(filePath)))]
        return(short_read)
    }))
    salmon_sr[, ntotal:=sum(counts), by = runname]
    setnames(salmon_sr, 'abundance','estimates')
    salmon_sr[, `:=`(
        length = NULL,
        countsFromAbundance = NULL)]
    salmon_sr[, TPM:=estimates]
    salmon_sr[, estimates:=TPM/1000000*ntotal]
    return(salmon_sr)
}


