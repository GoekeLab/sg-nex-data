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
        
        # if(short_read){
        #     rcMeanCountByRcWidth <- tapply(relWidth*assay(seDist1)[,1], rowData(seDist1)$GENEID, sum) # how to normalize gene read coverage??
        # }else{
            rcMeanCountByRcWidth <- tapply(relWidth*assay(seDist1)[,1], rowData(seDist1)$GENEID, sum)
        # }
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
    
    #&(gene_biotype %in% pro_types)],
    
    
    # txvec <- fread(paste0("/mnt/projects/SGNExManuscript/output/txList_matchingToGTF_wtChrIs.txt"), header = FALSE)
    # txvec <- gsub("\\..*","",txvec$V1)
    # ensemblAnnotations.transcripts <- copy(general_list$ensemblAnnotations.transcripts)
    # setnames(ensemblAnnotations.transcripts, "ensembl_gene_id","gene_name")
    # ensemblAnnotations.transcripts <- data.table(tx_name = txvec, status = TRUE)[ensemblAnnotations.transcripts, on = "tx_name"]
    # ensemblAnnotations.transcripts[is.na(status), status := FALSE]
    # ensemblAnnotations.transcripts[, all_in := all(status), by = gene_name]
    # genevec <- unique(ensemblAnnotations.transcripts[which(all_in)]$gene_name)
    #rl_data <- com_data[grepl("^ENSG",gene_name)&(gene_name %in% genevec)]
    rl_data <- dt[grepl("^ENSG",gene_name)&(gene_name %in% genevec)&(method %in% methodNames)]
    rl_data <- unique(ensemblAnnotations.transcripts[,.(gene_name, gene_biotype)])[rl_data, on = c("gene_name")]
    rl_data[, protocol_general := gsub("RandomPrimer", "", protocol_general)]
    
    
    # rl_data_gene <- unique(rl_data[, list(tpm=sum(normEst)), by = list(cellLine, protocol_general, gene_name, runname, gene_biotype, method)])
    # rl_data_gene_ave <- unique(rl_data_gene[, list(tpm = mean(tpm)), by = list(cellLine, protocol_general, gene_name, gene_biotype, method)])
    
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
    # #
    # samples <- copy(unique(general_list$samples[,.(runname, publicName,cellLine, protocol_type, cancer_type, Platform)]))
    # samples[protocol_type != "Illumina", runname := publicName]
    # samples[, publicName := NULL]
    runInfo <- samples[match(colnames(tmp),runname)]
    #runInfo <- unique(rl_data[,.(runname,cellLine, protocol_general)])[match(colnames(plotdata)[-1],runname)]
    
   
    return(list(runInfo, tmp))
}





complexHeatmap_plot <- function(countMatrix,runInfo, number_of_genes = 1000){
    sdvec <- apply(countMatrix,1,sd)
    #sd0.25 <- quantile(sdvec, prob = 0.75)
    # set.seed(2222)
    corMatrix <- cor(countMatrix[rank(-sdvec)<=number_of_genes,],method = 'spearman')
    #corMatrix <- cor(countMatrix[which(sdvec>sd0.25),],method = 'spearman')
    #corMatrix <- cor(countMatrix[sample(which(sdvec>sd0.25), 1000, replace = FALSE),],method = 'spearman')
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
    #plotMatrix <- as.matrix(tmp)
    colnames(corMatrix) <- NULL
    rownames(corMatrix) <- NULL
    #plotMatrix[sample(seq_len(nrow(plotMatrix)),5000)]
    p <- Heatmap(corMatrix, name = "Cor", col = col_fun, 
                 cluster_rows = TRUE,
                 cluster_columns =  TRUE,
                 # split = clusters$cluster,
                 # clustering_method_rows = "centroid" ,
                 # clustering_method_columns = "centroid",
                 #column_km = 7,
                 #clustering_distance_columns = "maximum",
                 top_annotation = cellLine_anno)#,
    #row_split = as.factor(rowInfo$cellLine),
    #cluster_row_slices = FALSE,
    #right_annotation = hgncTypes)
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
    # tmp2 <- limma::removeBatchEffect(tmp,batch = (runInfo$protocol_general=='Illumina'),batch2 = (runInfo$protocol_general=='cDNA')) # ,
    
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
    
    # print(p_heatmap)
}





process_replicate <- function(dt, methodNames, gene = TRUE, 
samples, ensemblAnnotations.transcripts, genevec, runnamevec, 
majorMinor = TRUE, scatterPlot = TRUE, complexity = TRUE, 
expressionLevel = TRUE, metric_type_id = 1,bpParameters){ #reproducibility_check = TRUE, 
   
    #rl_data <- com_data[grepl("^ENSG",gene_name)&(gene_name %in% genevec)]
    rl_data <- dt[grepl("^ENSG",gene_name)&(gene_name %in% genevec)&(method %in% methodNames)] # already filtered 
    rl_data <- unique(ensemblAnnotations.transcripts[,.(gene_name, gene_biotype)])[rl_data, on = c("gene_name")]
    rl_data[, protocol_general := gsub("RandomPrimer", "", protocol_general)]
    
    # rl_data_gene <- unique(rl_data[, list(tpm=sum(normEst)), by = list(cellLine, protocol_general, gene_name, runname, gene_biotype, method)])
    # rl_data_gene_ave <- unique(rl_data_gene[, list(tpm = mean(tpm)), by = list(cellLine, protocol_general, gene_name, gene_biotype, method)])
    filtered_rl_data <- rl_data[runname %in% runnamevec] 
    
    #vv <- CJ(k = 1:5, p = "directcDNA", s = seq_len(20), t = "MCF7")
    # kvar <- c("CPM", "uniqueCounts", "fullLengthCounts","uniqueCountsCPM", "fullLengthCountsCPM")
   
    cellLines <- c('Hct116','HepG2','K562','A549','MCF7','H9','HEYA8')
    
    source(paste0('/mnt/projects/SGNExManuscript/R/gene_cluster_code.R'))
    filtered_rl_data[, gene_cluster:=ifelse(gene_biotype %in% tr_gene_list,'TR gene',
                                    ifelse(gene_biotype %in% long_noncoding_rna_list, 'lncRNA',
                                           ifelse(gene_biotype %in% noncoding_rna_list, 'ncRNA',
                                                  ifelse(gene_biotype %in% pseudogene_list,'Pseudogene', ifelse(gene_biotype %in% ig_gene_list, 'IG gene',gene_biotype)))))]
    filtered_rl_data[, gene_cluster := ifelse(gene_cluster %in% c("IG gene","TR gene", "Mt_tRNA","Mt_rRNA","ribozyme"),"others", gene_cluster)]
    
    filtered_rl_data[, agg_gene_cluster := ifelse(gene_cluster == "processed_transcript", "lncRNA",
                                                 ifelse(gene_cluster %in% c("ncRNA","Pseudogene"), "others", gene_cluster))]
     # further group: IG gene, TR gene, IG gene, Mt gene, ribozyme together
    # IG gene               lncRNA              Mt_rRNA 
    # 213                14157                    2 
    # Mt_tRNA                ncRNA processed_transcript 
    # 22                 7557                  543 
    # protein_coding           Pseudogene             ribozyme 
    # 19847                14690                    8 
    # TR gene 
    # 197 
    protocolVec <- gsub("RandomPrimer","",unique(filtered_rl_data$protocol_general))
    cellLineList <- c(as.list(cellLines), list(cellLines))
    geneClusterList <- unique(filtered_rl_data$agg_gene_cluster)
    protein_coding_id <- grep("protein",geneClusterList)
    geneClusterList <- c(as.list(geneClusterList), list(geneClusterList))
    # if(reproducibility_check){
    #     combMat <- rbind(combn(1:4,1),
    #                      combn(1:4,1))
    # }else{
        combMat <- combn(1:4,2)
    # }
    
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
        source("/mnt/projects/SGNExManuscript/R/utility_function.R")
        np <- bplapply(vvIds[which(vv[vvIds]$t != 8)],pairwise_scatterplot_function , vv = vv, samples = samples,
                                            cellLineList = cellLineList, protocolVec = protocolVec, geneClusterList = geneClusterList, 
                                            combMat = combMat, filtered_rl_data = filtered_rl_data, gene = gene,
                                            majorMinor = majorMinor, 
                       complexity = complexity, expressionLevel = expressionLevel,
                       BPPARAM=bpParameters)
    }else{
        source("/mnt/projects/SGNExManuscript/R/utility_function.R")
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
    #print(paste(paste(cellLineV, collapse = " "), paste(protocolV, collapse = " "), collapse = ","))
    tmp <- filtered_rl_data[(cellLine %in% cellLineV)&(protocol_general %in% protocolV)&(agg_gene_cluster %in% geneCluster)]
    if(gene){
        tmp_wide <- dcast(tmp, gene_name ~ runname, value.var = "normEst")
    }else{
        tmp_wide <- dcast(tmp, tx_name + gene_name + ntx ~ runname, value.var = "normEst")
        if(majorMinor){
            #dominantTypeData <- majorMinor_generic(tmp)
            tmp_wide <- unique(dominant_typeData[cellLine %in% cellLineV,
                .(tx_name, gene_name,majorBoth, majorEither, majorEitherOnly, majorSecBoth, majorLongRead, majorShortRead
                )])[tmp_wide, on = "tx_name"]
        }
        
    }
    tmp_wide[is.na(tmp_wide)] <- 0
    ## pairwise correlation is calculated for transcripts being expressed in 
    ## either sample
    #combMat <- combn(seq_len(ncol(tmp_wide))[-1], 2)
    nameMat <- CJ(v1 = colnames(tmp_wide)[grep(paste0("_",protocolV[1]),colnames(tmp_wide))],
                  v2 = colnames(tmp_wide)[grep(paste0("_",protocolV[2]),colnames(tmp_wide))])
    # remove the one with the same name
    nameMat <- nameMat[v1 != v2]
    # remove the duplciated pair
    #nameMat[, v1v2 := paste(sort(c(v1,v2)), collapse = ""), by = list(v1,v2)]
    #nameMat[, duplicated_status := duplicated(v1v2)]
    #nameMat <- nameMat[which(!duplicated_status)]
    nameMat[, rep_status := (samples[which(samples$runname == v1)]$bioRep ==  samples[which(samples$runname == v2)]$bioRep), by = list(v1,v2)]
    nameMat[, rep := samples[which(samples$runname == v1)]$bioRep, by = v1]

    #nameMat <- nameMat[which(rep_status)]
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
    # print("debug 1")
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
           
            # gene quantile based: 1, (1,4], (4,193]
            
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
            
            
            # maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx <=1],"<=1", maeValues)
            # maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >1 &(ntx<=4)],"(1,4]", maeValues)
            # maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >4 &(ntx<=193)],"(4,193]", maeValues)
            
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
            
            
            # maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx <=1],"<=1", maeValues)
            # maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >1 &(ntx<=4)],"(1,4]", maeValues)
            # maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >4 &(ntx<=193)],"(4,193]", maeValues)
            
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
            
            
            # maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx <=1],"<=1", maeValues)
            # maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >1 &(ntx<=4)],"(1,4]", maeValues)
            # maeValues <-calc_mae(t,g,nameMat,cellLineList,protocolV, geneClusterList,tmp_wide[ntx >4 &(ntx<=193)],"(4,193]", maeValues)
            
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
    # log2Diff <- apply(log2CorData,1,diff, na.rm = TRUE)
    # dVec <- abs(log2Diff)
    # rd <- dVec/apply(log2CorData,1,mean, na.rm = TRUE)
    # dRankVec <- abs(rank(log2CorData$V1)-rank(log2CorData$V2))
    # diagDist <- dVec/sqrt(2)
    r2 <- summary(lm(V2~V1, data = log2CorData))$r.squared
    # rmse <- sqrt(mean(log2Diff^2))
    print(protocolV)
    temp_corValues <- data.table(r_expressed = cor(log2CorData[expressed], method = "spearman")[1,2],#)[1,2]
                                 r = cor(log2CorData, method = "spearman")[1,2], #)[1,2]
                                 # d_mean = mean(dVec, na.rm = TRUE),
                                 # d_sd = sd(dVec, na.rm = TRUE),
                                 # d_zscore_mean = mean(dVec, na.rm = TRUE)/sd(dVec, na.rm = TRUE),
                                 # d_coefvar = sd(dVec, na.rm = TRUE)/mean(dVec, na.rm = TRUE),
                                 r2 = r2,
                                 # rmse = rmse,
                                 ne = length(expressed),
                                 #bioRep = nameMat[x]$rep,
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
    # for( x in seq_len(nrow(nameMat))){
    #     corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
    #     setnames(corData, c(1:2), c("V1","V2"))
    #     expressed <- which(apply(corData>0,1,any)) # expressed in both
    #     log2CorData <- log2(corData+1)
    #     log2Diff <- apply(log2CorData,1,diff, na.rm = TRUE)
    #     dVec <- abs(log2Diff)
    #     # rd <- dVec/apply(log2CorData,1,mean, na.rm = TRUE)
    #     # dRankVec <- abs(rank(log2CorData$V1)-rank(log2CorData$V2))
    #     # diagDist <- dVec/sqrt(2)
    #     r2 <- summary(lm(V2~V1, data = log2CorData))$r.squared
    #     rmse <- sqrt(mean(log2Diff^2))
    #     temp_corValues <- data.table(r_expressed = cor(log2CorData[expressed])[1,2],# method = "spearman")[1,2],
    #                                  r = cor(log2CorData)[1,2],#, method = "spearman")[1,2],
    #                                  d_mean = mean(dVec, na.rm = TRUE),
    #                                  d_sd = sd(dVec, na.rm = TRUE),
    #                                  d_zscore_mean = mean(dVec, na.rm = TRUE)/sd(dVec, na.rm = TRUE),
    #                                  d_coefvar = sd(dVec, na.rm = TRUE)/mean(dVec, na.rm = TRUE),
    #                                  r2 = r2,
    #                                  rmse = rmse,
    #                                  ne = length(expressed),
    #                                  #bioRep = nameMat[x]$rep,
    #                                  match_status = nameMat[x]$rep_status,
    #                                  common_type = typeName,
    #                                  n = nrow(dd))
        corValues <- do.call("rbind", list(corValues, temp_corValues))
    # }
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
    #print(paste(paste(cellLineV, collapse = " "), paste(protocolV, collapse = " "), collapse = ","))
    tmp <- filtered_rl_data[(cellLine %in% cellLineV)&(protocol_general %in% protocolV)&(gene_cluster %in% geneCluster)]
    if(gene){
        tmp_wide <- dcast(tmp, gene_name ~ runname, value.var = "normEst")
    }else{
        tmp_wide <- dcast(tmp, tx_name + ntx + gene_name ~ runname, value.var = "normEst")
        if(majorMinor){
            #dominantTypeData <- majorMinor_generic(tmp)
            tmp_wide <- unique(dominant_typeData[,.(tx_name, gene_name,majorBoth, majorEither, majorEitherOnly,
                                                   majorSecBoth,majorFirstSecBoth, majorLongRead, majorShortRead
            )])[tmp_wide, on = "tx_name"]
        }
    }
    tmp_wide[is.na(tmp_wide)] <- 0
    ## pairwise correlation is calculated for transcripts being expressed in 
    ## either sample
    #combMat <- combn(seq_len(ncol(tmp_wide))[-1], 2)
    nameMat <- CJ(v1 = colnames(tmp_wide)[grep(paste0("_",protocolV[1]),colnames(tmp_wide))],
                  v2 = colnames(tmp_wide)[grep(paste0("_",protocolV[2]),colnames(tmp_wide))])
    nameMat[, rep_status := (samples[runname == v1]$bioRep ==  samples[runname == v2]$bioRep), by = list(v1,v2)]
    nameMat[, rep := samples[runname == v1]$bioRep, by = v1]
    
    #nameMat <- nameMat[which(rep_status)]
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
    # return(pList)
    
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
        # if(majorMinor){
        #     p_scatter <- ggplot(data = corData, aes(x = log2(V1+1), y = log2(V2+1)))+
        #         # xlim(0,1)+
        #         # ylim(0,1)+
        #         # expand_limits(x = c(0,1), y = c(0,1))+
        #         xlab(varnames[1])+
        #         ylab(varnames[2])+
        #         geom_hex(aes(fill = stat(count)),
        #                  binwidth = 0.01) +
        #         scale_fill_gradient(name = 'count', low = "grey", high = "black")+
        #         #labels = c('0', '1', '2', '3','4+')
        #         #)+
        #         stat_cor(aes(label = ..r.label..),method = "spearman", cor.coef.name = "Sp.R")+
        #         #facet_wrap(~common_type)+
        #         ggtitle(pNames[x])+
        #         theme_classic()
        # }else{
            p_scatter <- ggplot(data = corData, aes(x = log2(V1+1), y = log2(V2+1)))+
                # xlim(0,1)+
                # ylim(0,1)+
                # expand_limits(x = c(0,1), y = c(0,1))+
                xlab(varnames[1])+
                ylab(varnames[2])+
                geom_hex(aes(fill = stat(count)),
                         binwidth = 0.25) +
                scale_fill_gradient(name = 'count', low = "grey", high = "black")+
                #labels = c('0', '1', '2', '3','4+')
                #)+
                stat_cor(aes(label = ..r.label..),method = "spearman", cor.coef.name = "Sp.R")+
                ggtitle(pNames[x])+
                # facet_wrap(~common_type)+
                theme_classic()
        # }
        pList[[x]] <- p_scatter
    }
    
    names(pList) <- pNames
    saveplot.dir <- "/mnt/projects/SGNExManuscript/output/replicate_plot_complexity/" 
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
        # match_status = nameMat$rep_status,
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
        # expressed <- which(apply(corData,1,sum)>0)
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
        # match_status = nameMat$rep_status,
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
        # match_status = nameMat$rep_status,
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
        # match_status = nameMat$rep_status,
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
        # expressed <- which(apply(corData,1,sum)>0)
        diffMat[(corData$V1+corData$V2)>0,x] = abs(corData[(corData$V1+corData$V2)>0]$V1-corData[(corData$V1+corData$V2)>0]$V2)/(corData[(corData$V1+corData$V2)>0]$V1+corData[(corData$V1+corData$V2)>0]$V2)
        
    }
    
    #
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
        # match_status = nameMat$rep_status,
        common_type = post_fix,
        protocol_comparison = paste(protocolV, collapse = " vs "))
    ]
    
    maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
    return(maeValues)
}

# rd_mean = mean(rd, na.rm = TRUE),
# rd_sd = sd(rd, na.rm = TRUE),
# rd_zscore_mean = mean(rd, na.rm = TRUE)/sd(rd, na.rm = TRUE),
# rd_coefvar = sd(rd, na.rm = TRUE)/mean(rd, na.rm = TRUE),
# diagDist_mean = mean(diagDist, na.rm = TRUE),
# diagDist_sd = sd(diagDist, na.rm = TRUE),
# diagDist_zscore_mean = mean(diagDist, na.rm = TRUE)/sd(diagDist, na.rm = TRUE),
# diagDist_coefvar = sd(diagDist, na.rm = TRUE)/mean(diagDist, na.rm = TRUE),
# dRank_mean = mean(dRankVec, na.rm = TRUE),
# dRank_median = median(dRankVec, na.rm = TRUE),


# rd_mean = corValues$rd_mean,
# rd_sd = corValues$rd_sd,
# rd_zscore_mean = corValues$rd_zscore_mean,
# rd_coefvar = corValues$rd_coefvar,
# diagDist_mean = corValues$diagDist_mean,
# diagDist_sd = corValues$diagDist_sd,
# diagDist_zscore_mean = corValues$diagDist_zscore_mean,
# diagDist_coefvar = corValues$diagDist_coefvar,
# dRank_mean = corValues$dRank_mean,
# dRank_median = corValues$dRank_median,


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
    # dominant_typeData <- dcast(unique(rl_data_tx_ave[,
    #     .(tx_name, gene_name, cellLine, isoform_rank, short_read, geneExpressedInBoth)]), 
    #     tx_name + gene_name + cellLine + geneExpressedInBoth ~ short_read, value.var = "isoform_rank")
    
    rl_data_tx_ave[, majorBoth := all(isoform_rank==1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorEither := any(isoform_rank==1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorEitherOnly := (sum(isoform_rank==1)==1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorSecBoth := all(isoform_rank==2)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorFirstSecBoth := (all(isoform_rank==2)|all(isoform_rank == 1))&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorLongRead := any(short_read==FALSE&isoform_rank == 1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    rl_data_tx_ave[, majorShortRead := any(short_read==TRUE&isoform_rank == 1)&(geneExpressedInBoth), by = list(tx_name, gene_name,cellLine)]
    return(rl_data_tx_ave)
}

#dominant_typeData <- identifyMajorMinorIsoforms(com_data, ensemblAnnotations.transcripts)

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
    # dominant_typeData <- dcast(unique(rl_data_tx_ave[,
    #                                                  .(tx_name, gene_name, isoform_rank, protocol_general, geneExpressedInBoth)]), 
    #                            tx_name + gene_name + geneExpressedInBoth ~ protocol_general, value.var = "isoform_rank")
    # 
    rl_data_tx_ave[, majorBoth := all(isoform_rank==1)&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    rl_data_tx_ave[, majorEither := any(isoform_rank==1)&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    rl_data_tx_ave[, majorEitherOnly := (sum(isoform_rank==1)==1)&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    rl_data_tx_ave[, majorSecBoth := all(isoform_rank==2)&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    rl_data_tx_ave[, majorFirstSecBoth := (all(isoform_rank==2)|all(isoform_rank == 1))&(geneExpressedInBoth), by = list(tx_name, gene_name)]
    return(rl_data_tx_ave)
}


repeat_analysis_function <- function(seOutput, ervRanges, retrotransposonOnly = FALSE){
    #seOutput <- readRDS("/mnt/data/spikeinSe/bambuOutput_LR_allSgNexSample.rds")
   # take the alignment type out 
    #com_data[, alignment_type := ifelse(grep("_uniAln",))]
    com_data <- process_seOutput(seOutput)
    
    # com_data_filter_fullLength <- com_data#[ntotal>400000]
    
    extendedAnnotationGRangesList <- rowRanges(seOutput)
    
    isoTEratio <- estimate_repeat_ratio(extendedAnnotationGRangesList,anno_exByTx,ervRanges, txLengths.tbldf, retrotransposonOnly = retrotransposonOnly)
    
    # isoEst <- data.table(assays(seOutput)$counts, keep.rownames = TRUE)
    # setnames(isoEst,"rn","tx_name")
    # isoEst <- melt(isoEst, id.var = "tx_name", measure.vars = colnames(isoEst)[-1])
    # isoEst[, runname:=gsub("HCT116","Hct116",gsub("(.genome_alignment.sorted)|(_R1.sorted)|(_sorted)","",
    #                                               gsub("-pre|-Pre","",variable))), by = variable]
    # isoEst[runname == "GIS_Hct116_cDNA_Rep2_Run4", runname:="GIS_Hct116_cDNA_Rep2_Run5"]
    # 
    # isoEst[,`:=`(estimates = sum(value)), by = list(tx_name,runname)]
    # isoEst <- unique(isoEst[,.(tx_name, runname, estimates)])
    # isoEst[, `:=`(ntotal = sum(estimates)), by = runname]
    # 
    # 
    # cellLineVec <- c("A549","K562","MCF7","Hct116","HepG2","H9","HEYA8")
    # all_samples <- gsub("_sorted|_R1.sorted","",colnames(seOutput))
    # dt <- data.table(runname = all_samples)
    # dt[, cellLine := ifelse(grepl("uniAln$",runname), "uni_aln_filter", 
    #                         ifelse(grepl("priAln$",runname), "pri_uni_aln_filter", "no_filter")), by = runname]
    # #dt[, cellLine:=gsub('k562','K562',strsplit(runname, '\\_')[[1]][2]),by = runname]
    # dt[, protocol:=strsplit(runname, '\\_')[[1]][3], by = runname]
    # dt[, cDNAstranded:=ifelse(protocol %in% c('cDNA','cDNAStranded'), protocol=='cDNAStranded',NA)]
    # dt[, randomPrimer:=grepl('RandomPrimer',protocol)]
    # dt[, protocol_type:=gsub('Stranded|RandomPrimer','',gsub('PromethionD','d', protocol))]
    # dt[, repInfo:=strsplit(runname, '\\_')[[1]][4], by = runname]
    # dt[grep("GIS_Hct116_directRNA_[1-9]|(GIS_Hct116_directcDNA_[1-9])|(GIS_Hct116_cDNA_[1-9])",runname), repInfo:=paste0("Rep1-Run",repInfo)]
    # dt[, bioRep:=strsplit(repInfo, '\\-')[[1]][1], by = repInfo]
    # 
    # dt[, bioRep:=gsub('Rep','',bioRep)]
    # dt[, techRep:=strsplit(repInfo, '\\-')[[1]][2], by = repInfo]
    # dt[is.na(techRep), techRep:=1]
    # dt[, techRep:=gsub('Run','',techRep)]
    # dt[, patient_derived:=(!(cellLine %in% cellLineVec))]
    # 
    # cellLines <- unique(dt$cellLine)
    # 
    # isoEst_filter <- isoEst#[ntotal>400000]
    # isoEst_filter[, normEstimates:=(estimates/ntotal*10^6), by = runname]
    # 
    # # tmp <- isoEst_filter[,.I[which(any(normEstimates>20))], by = tx_name]
    # # plot_tmp <- unique(isoTEratio[!(grepl("tx.",tx_name)&(anno_status=="annotated"))][tx_name %in% tmp$tx_name][,.(tx_name, repRatio_corrected,anno_status, newTxClassAggregated)])
    # # 
    # isoEst_filter <- dt[isoEst_filter, on = "runname"]
    # isoEst_filter_wide <-  dcast(isoEst_filter[cellLine %in% cellLines], tx_name ~ runname, value.var = "normEstimates")
    # 
    # isoTEratio[, repRatioByRepClass:=sum(averageRepRatio), by = list(tx_name, rep_class,strand)]
    # 
    # isoTEratio_agg <- unique(isoTEratio[!(grepl("tx.",tx_name)&(anno_status=="annotated"))][tx_name %in% tmp$tx_name,.(repRatioByRepClass, tx_name, rep_class, strand, anno_status, repRatioIso_all, gene_name)])
    # 
    #return(list(isoTEratio, isoTEratio_agg, isoEst_filter_wide))
    return(list(com_data, isoTEratio))
}


process_seOutput <- function(seOutput){
    tmp <- data.table(as.data.frame(rowData(seOutput)))
    #tmp[!grepl("unspliced",txClassDescription)&(grepl("new",txClassDescription))]
    
    
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
    
    # incompatibleDt[, runname:=gsub("HCT116","Hct116",gsub("(.genome_alignment.sorted)|(_R1.sorted)|(_sorted)","",gsub("-pre|-Pre","",runname))), by = runname]
    # incompatibleDt[runname == "GIS_Hct116_cDNA_Rep2_Run4", runname:="GIS_Hct116_cDNA_Rep2_Run5"]
    # incompatibleDt[, incompatibleCounts:=sum(incompatibleCounts), by = runname]
    # incompatibleDt <- unique(incompatibleDt, by = NULL)
    # totalDt <- unique(genomeem_lr[,.(runname, ntotal)], by = NULL)
    # totalDt <- incompatibleDt[totalDt, on = "runname"]
    # totalDt[, ntotal := ntotal+incompatibleCounts]
    # genomeem_lr[, ntotal := NULL]
    # genomeem_lr <- totalDt[genomeem_lr, on = "runname"]
    
    genomeem_lr <- geneTxTable[genomeem_lr, on = "tx_name"]
    genomeem_lr[, method := "genomeem_lr"]
    
    #genomeem_lr <- genomeem_lr[runname %in% sampleNames]
    
    merge.colnames <- c('tx_name','gene_name','fullLengthCounts','uniqueCounts','totalCounts','runname','method','ntotal')
    com_data <- genomeem_lr[, merge.colnames, with =FALSE]
    #salmon_sr[, merge.colnames, with =FALSE])
    
    com_data[, runname := as.character(runname)]
    com_data[, protocol:=unlist(strsplit(runname, 
                                         '_'))[3],by= runname]
    
    com_data[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
    com_data[, protocol_method:=paste0(protocol_general,'.',method)]
    com_data[, short_read:=as.numeric(protocol_general == 'Illumina')]
    com_data[, runname_method:=paste0(runname,'.',method)]
    #com_data[, cellLine := gsub('k562','K562',unlist(strsplit(runname,'_'))[2]), by = runname]
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
    #if(grepl("newTxClass", colnames(geneTxTable_extended))){
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
    # all(names(grl) == sort(unique(gr$hgnc))
    
    # reduce each element (independently reduce ranges for each gene)
    grl_redux <- reduce(grl) # element-wise, like lapply(grl, reduce) 
    # all(names(grl_redux) == names(grl)) & all(lengths(grl_redux) <= lengths(grl))
    
    # return single GRanges, with rownames derived from hgnc
    ervRangesReduceByRepClass <- unlist(grl_redux)
    ervRangesReduceByRepClass$repClass <- names(ervRangesReduceByRepClass)
    
    #ov <- compute_overlap(exons_granges, ervRanges,ignore.strand, by_rep_type = TRUE)
    if(retrotransposonOnly){
        ervRangesReduceByRepClass <- ervRangesReduceByRepClass[grep("LINE|SINE|LTR",ervRangesReduceByRepClass$repClass)]
    }
    
    ovByRepClass <- compute_overlap(exons_granges, ervRangesReduceByRepClass,ignore.strand, by_rep_type = TRUE)
    
    # isoTEratio_all <- unique(ov[, list(averageRepRatio_byisoform = sum(repRatio_all),
    #                                                              #averageRepRatio = sum(repRatio)/nAnnotatedExon,
    #                                                              # strand = strand,
    #                                                              rep_class = rep_class), by = list(txId, rep_name)])
    ## by rep class 
    isoTEratio_all_byrepclass <- unique(ovByRepClass[, list(repRatio_byisoform = sum(repRatio_all)), by = list(txId, rep_class)])
    isoTEratio_byrepclass <- unique(ovByRepClass[, list(repRatio = sum(repRatio)), by = list(txId, rep_class, anno_status)])
    isoTEratio_byrepclass <-isoTEratio_all_byrepclass[isoTEratio_byrepclass, on = c("txId","rep_class")]
    # isoTEratio_anno[, repRatioIso_all := sum(averageRepRatio), by = txId] # sum by repeat name 
    # isoTEratio_novel[, repRatioIso_all := sum(averageRepRatio), by = txId]# sum by repeat name
    # 
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
    system(paste0("/mnt/projects/gffcompare/gffcompare -TNRQS -e 15 -d 15 -o gffCompare ", 
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
                         # anno_status = exons_granges[qHits]$anno_status,
                         # nNovelExon = exons_granges[qHits]$nnovel,
                         # nAnnotatedExon = exons_granges[qHits]$nannotated,
                         #rep_id = names(ervRanges)[sHits],
                         # strand = as.character(strand(ervRanges[sHits])),
                         #rep_name = ervRanges[sHits]$repName,
                         rep_class = tmpErvRanges[sHits]$repClass,
                         #rep_family = ervRanges[sHits]$repFamily,
                         # exon_width = width(exons_granges)[qHits],
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
                         # anno_status = exons_granges[qHits]$anno_status,
                         # nNovelExon = exons_granges[qHits]$nnovel,
                         # nAnnotatedExon = exons_granges[qHits]$nannotated,
                         # exon_width = width(exons_granges)[qHits],
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
        # print(k)
        # if(grepl("H9|HEYA8",k)){ #v what's the difference between old and new ones?
        #     filePath <- sort(dir(paste0('/mnt/projects/SGNExManuscript/output/sr/02_Mapping/',k,'/transcripts_quant'),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]
        # }else{
        #     filePath <- sort(dir(paste0('/mnt/projects/SGNExManuscript/output/sr_new/02_Mapping_matchedToGTF/',k,'/transcripts_quant_biasCorrected'),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]
        # }
        # 
        # if(length(filePath)==0){
        #     return(NULL)
        # }
        print(filePath)
        txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
        names(txi)
        
        short_read <- data.table(tx_name = rownames(txi$abundance),
                                 abundance = txi$abundance[,1],
                                 counts = txi$counts[,1],
                                 length = txi$length[,1],
                                 countsFromAbundance = txi$countsFromAbundance)
        # short_read <- fread(filePath, header = TRUE)
        short_read[, runname:=basename(gsub("\\/transcripts_quant","",dirname(filePath)))]
        return(short_read)
    }))
    salmon_sr[, ntotal:=sum(counts), by = runname]
    setnames(salmon_sr, 'abundance','estimates')
    salmon_sr[, `:=`(#counts = NULL,
        length = NULL,
        countsFromAbundance = NULL)]
    salmon_sr[, TPM:=estimates]
    salmon_sr[, estimates:=TPM/1000000*ntotal]
    return(salmon_sr)
}


# trim_lr_150bp <- 
#     corData <- lapply(cellLines, function(s){
#         corData <- do.call("rbind",lapply(protocolVec[1:3], function(k){
#             plotdata <- dcast(rl_data_tx_ave[cellLine == s &(protocol_general %in% c(k, "Illumina"))], tx_name + gene_biotype ~ protocol_general, value.var = "tpm_log")
#             mat <- as.matrix(plotdata[gene_biotype %in% pro_types][,c(3,4), with = FALSE])
#             mat[is.na(mat)] <- 0 
#             return(data.table(cellLine = s,
#                               protocol = k, 
#                               r = cor(mat[,1], mat[,2],method = "spearman")))
#         }))
#         return(corData)
#     })
# corData <- do.call("rbind",corData)
# 
# # instead of by protocol, should be long read vs short read 
# isoform_types <- c("major","either","major_complimentary","all")
# txCor <- do.call("rbind",lapply(seq_len(ncol(protocol_combinations)), function(k){
#     pv <- protocolVec[protocol_combinations[,k]]
#     dominant_typeData <- dcast(unique(rl_data_tx_ave[protocol_general %in% pv,.(tx_name, gene_name, cellLine, isoform_rank, protocol_general)]), tx_name + gene_name + cellLine ~ protocol_general, value.var = "isoform_rank")
#     
#     rl_data_tx_ave_pv <- dominant_typeData[rl_data_tx_ave, on = c("tx_name","gene_name","cellLine")]
#     
#     txCor <- do.call("rbind",lapply(cellLines, function(s){
#         print(k)
#         print(s)
#         
#         tt <- copy(rl_data_tx_ave_pv[cellLine == s &(gene_biotype %in% pro_types)&(protocol_general %in% pv)&(geneExpression>0)])
#         setnames(tt, pv, c("protocol1","protocol2"))
#         
#         txCor <- do.call("rbind",lapply(1:4, function(l){
#             if(l==1){
#                 plotdata <- dcast(tt[((protocol1 == 1)&(protocol2 == 1))], gene_name + tx_name + gene_biotype ~ protocol_general, value.var = "tpm")
#                 
#             }
#             if(l == 3){
#                 plotdata <- dcast(tt[!((protocol1 <= 1)&(protocol2 <= 1))], gene_name + tx_name + gene_biotype ~ protocol_general, value.var = "tpm")
#             }
#             if(l == 2){
#                 plotdata <- dcast(tt[((protocol1 == 1)|(protocol2 == 1))], gene_name + tx_name + gene_biotype ~ protocol_general, value.var = "tpm")
#             }
#             if(l == 4){
#                 plotdata <- dcast(tt, gene_name + tx_name + gene_biotype ~ protocol_general, value.var = "tpm")
#             }
#             setnames(plotdata, pv, c("protocol1","protocol2"))
#             mat <- as.matrix(plotdata[,c("protocol1","protocol2"), with = FALSE])
#             mat[is.na(mat)] <- 0 
#             return(data.table(cellLine = s,
#                               protocol_comparison = paste(gsub("Illumina","Illu",gsub("direct","d",pv)), collapse = " vs "),
#                               n_isoform = nrow(mat),
#                               isoform_type = isoform_types[l],
#                               r = cor(mat[,1], mat[,2],method = "spearman")))
#             
#         }))
#         return(txCor)
#     }))
#     return(txCor)
# }))
# 
# noprint <- lapply(1:4, function(l){
#     breaks_set <- list(c(0.85,0.9,0.95),
#                        c(0.65,0.75,0.85,0.95),
#                        c(0.4,0.6,0.8),
#                        c(0.4,0.6,0.8))[[l]]
#     breaks_limit <- list(c(0.85,0.95),
#                          c(0.65,0.95),
#                          c(0.4,0.8),
#                          c(0.4,0.81))[[l]]
#     iso_type <- isoform_types[l]
#     p <- ggplot(txCor[isoform_type == iso_type], aes(x = protocol_comparison, y = r))+
#         geom_boxplot(outlier.shape = NA)+
#         geom_jitter(aes(col = cellLine), size = 2, pch = 1)+
#         scale_y_continuous(breaks = breaks_set, limits = breaks_limit)+
#         scale_x_discrete(limits =  unique(corDataGene$protocol_comparison)[c(1,2,4,3,5,6)])+
#         xlab("")+
#         coord_flip()+
#         scale_color_brewer(type = "qual", palette = 3)+
#         ylab("Spearman correlation of transcript expression estimates")+
#         ggtitle(iso_type)+
#         theme_classic()
#     pdf(paste0("figures/tx_correlation_boxplot_celllineprotocol_comparison_",iso_type,".pdf"), width = 8, height = 6)
#     print(p)
#     dev.off()
# })
# 
# 
# 
# 
# saveRDS(corData, file = paste0("output/LRvsSR_correlationData_tx.rds"))
# 
# corData_gene <- readRDS(paste0("output/LRvsSR_correlationData_gene.rds"))
# corData_tx <- readRDS(paste0("output/LRvsSR_correlationData_tx.rds"))
# 
# corData_gene[, feature := "gene"]
# corData_tx[, feature := "tx"]
# 
# corData <- do.call("rbind", list(corData_gene, corData_tx))
# 
# protocolCol <- adjustcolor(brewer.pal(8,"Dark2")[1:4],0.7)
# protocolVec <-  c("directRNA","directcDNA","cDNA","Illumina")
# protocolLabel <- c("RNA","PCR-free cDNA","cDNA","Illumina")
# 
# p <- ggplot(corData, aes(x = feature, y = r))+
#     geom_boxplot(aes(col = feature))+
#     geom_jitter(aes(col = feature),position=position_jitter(0.2), size = 2)+
#     ylab("Spearman correlation")+
#     xlab("Feature type")+
#     scale_color_brewer(type = "qual", guide = "none")+
#     theme_classic()
# 
# pdf(paste0("figures/cor_plot.pdf"), width = 4, height = 3)
# print(p)
# dev.off()
# 
# pro_types <- c("protein_coding","antisense_RNA","lincRNA","non_coding","macro_lncRNA")
# 
# plotdata <- dcast(rl_data_tx_ave[cellLine == "A549"&(gene_biotype %in% pro_types) &(protocol_general %in% c("directcDNA", "Illumina"))], gene_name + tx_name + gene_biotype ~ protocol_general, value.var = "tpm")
# 
# p1 <- ggplot(plotdata, aes(x=log2(Illumina+1), y=log2(directcDNA+1)))+
#     geom_hex(aes(fill = stat(cut(log(count), breaks = log(c(1, 10, 100, 1000,10000,Inf)), labels = F, right = T, include.lowest = T))), binwidth = 0.1) +
#     scale_fill_gradient(name = 'count(log10)', low = "light blue", high = "steelblue", labels = c('0', '1', '2', '3','4+'))+
#     xlab('Short read estimates (salmon)')+
#     ylab('Long read estimates (Bambu)')+
#     ggpubr::stat_cor(aes(label = ..r.label..),method = "spearman", cor.coef.name = "Sp.R")+
#     theme_classic()
# 
# p1
# pdf(paste0("figures/tx_level_1run_LRvsSR.pdf"), width = 8, height = 6)
# print(p1)
# dev.off()
# 
# 
# protocol_types <- c("directcDNA", "Illumina")
# dominant_typeData <- dcast(unique(rl_data_tx_ave[protocol_general %in% protocol_types,.(tx_name, gene_name, cellLine, isoform_rank, protocol_general)]), tx_name + gene_name + cellLine ~ protocol_general, value.var = "isoform_rank")
# 
# rl_data_tx_ave <- dominant_typeData[rl_data_tx_ave, on = c("tx_name","gene_name","cellLine")]
# saveRDS(rl_data_tx_ave, file = "output/rl_data_tx_ave.rds")
# pro_types <- c("protein_coding","antisense_RNA","lincRNA","non_coding","macro_lncRNA")
# #[gene_biotype %in% c("protein_coding")]
# #cellLineRepVec <- unique(rl_data$cellLineRep)
# plotdata_major <- dcast(rl_data_tx_ave[cellLine == "A549" &(gene_biotype %in% pro_types)&(protocol_general %in% protocol_types)&((directcDNA == 1)&(Illumina == 1))&(geneExpression>0)], gene_name + tx_name + gene_biotype ~ protocol_general, value.var = "tpm")
# plotdata_major[is.na(Illumina), Illumina := 0]
# plotdata_major[is.na(directcDNA), directcDNA := 0]
# 
# plotdata_others <- dcast(rl_data_tx_ave[cellLine == "A549" &(gene_biotype %in% pro_types)&(protocol_general %in% protocol_types)&!((directcDNA <= 1)&(Illumina <= 1))&(geneExpression>0)], gene_name + tx_name + gene_biotype ~ protocol_general, value.var = "tpm")
# plotdata_others[is.na(Illumina), Illumina := 0]
# plotdata_others[is.na(directcDNA), directcDNA := 0]
# 
# plotdata_either <- dcast(rl_data_tx_ave[cellLine == "A549" &(gene_biotype %in% pro_types)&(protocol_general %in% protocol_types)&((directcDNA == 1)|(Illumina == 1))&(geneExpression>0)], gene_name + tx_name + gene_biotype ~ protocol_general, value.var = "tpm")
# plotdata_either[is.na(Illumina), Illumina := 0]
# plotdata_either[is.na(directcDNA), directcDNA := 0]
# 
# p1 <- ggplot(plotdata_major, aes(y=log2(directcDNA+1), x=log2(Illumina+1)))+
#     geom_hex(
#         aes(fill =
#                 #stat(count)),
#                 stat(cut(log(count), breaks = log(c(1, 10, 100, 1000,10000,Inf)), labels = F, right = T, include.lowest = T))),
#         binwidth = 0.1) +
#     scale_fill_gradient(name = 'count(log10)', low = "light blue", high = "steelblue",#)+
#                         labels = c('0', '1', '2', '3','4+'))+#))+#,
#     #labels = c('0', '1', '2', '3'))+
#     xlab('Short read estimates (Salmon)')+
#     ylab('Long read estimates (Bambu)')+
#     ggpubr::stat_cor(aes(label = ..r.label..),method = "spearman", cor.coef.name = "Sp.R")+
#     theme_classic()
# p1
# pdf(paste0("figures/tx_level_1run_LRvsSR_major.pdf"), width = 8, height = 6)
# print(p1)
# dev.off()
# 
# 
# p1 <- ggplot(plotdata_others, aes(y=log2(directcDNA+1), x=log2(Illumina+1)))+
#     geom_hex(
#         aes(fill =
#                 #stat(count)),
#                 stat(cut(log(count), breaks = log(c(1, 10, 100, 1000,10000,Inf)), labels = F, right = T, include.lowest = T))),
#         binwidth = 0.1) +
#     scale_fill_gradient(name = 'count(log10)', low = "light blue", high = "steelblue",#)+
#                         labels = c('0', '1', '2', '3','4+'))+#))+#,
#     #labels = c('0', '1', '2', '3'))+
#     xlab('Short read estimates (Salmon)')+
#     ylab('Long read estimates (Bambu)')+
#     ggpubr::stat_cor(aes(label = ..r.label..),method = "spearman", cor.coef.name = "Sp.R")+
#     theme_classic()
# p1
# pdf(paste0("figures/tx_level_1run_LRvsSR_others.pdf"), width = 8, height = 6)
# print(p1)
# dev.off()


## wrong calculattion:
#metrics like mae, mard and rmse, should be per replicate pair instead of per transcript
# this is essentially the same as the one used in salmon: median(log2(estimated/truth))
# calc_mae <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
#     diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
#     for( x in seq_len(nrow(nameMat))){
#         corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
#         varnames <- colnames(corData)
#         setnames(corData, c(1:2), c("V1","V2"))
#         # expressed <- which(apply(corData,1,sum)>0)
#         diffMat[(corData$V1+corData$V2)>0,x] = abs(log2(corData[(corData$V1+corData$V2)>0]$V1+1)-log2(corData[(corData$V1+corData$V2)>0]$V2+1))
#         
#     }
#     
#     #
#     if("tx_name" %in% colnames(dd)){
#         temp_maeValues <- dd[,.(tx_name, gene_name, ntx)]
#     }else{
#         temp_maeValues <- dd[,.(gene_name)]
#     }
#     if(expressionLevel) temp_maeValues[, expressed := apply(corData>expression_t,1,sum)>0]
#     temp_maeValues[,mae := apply(diffMat,1,mean,na.rm = TRUE)]
#     temp_maeValues[,`:=`(match_mae = NA, non_match_mae = NA)]
#     if(length(which(nameMat$rep_status))>1){
#         temp_maeValues[,match_mae := apply(diffMat[, which(nameMat$rep_status)],1,mean,na.rm = TRUE)]
#     }else if(length(which(nameMat$rep_status))==1){
#         temp_maeValues[,match_mae := mean(diffMat[, which(nameMat$rep_status)],na.rm = TRUE)]
#     }
#     if(length(which(!nameMat$rep_status))>1){
#         temp_maeValues[,non_match_mae := apply(diffMat[, which(!nameMat$rep_status)],1,mean,na.rm = TRUE)]
#     }else if(length(which(!nameMat$rep_status))==1){
#         temp_maeValues[,match_mae := mean(diffMat[, which(!nameMat$rep_status)],na.rm = TRUE)]
#     }
#     
#     temp_maeValues[, `:=`(
#         cellLine = c(cellLineList[[8]],"all")[t],
#         gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
#         # match_status = nameMat$rep_status,
#         common_type = post_fix,
#         protocol_comparison = paste(protocolV, collapse = " vs "))
#     ]
#     
#     maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
#     return(maeValues)
# }
# 
# calc_rmse <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
#     diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
#     for( x in seq_len(nrow(nameMat))){
#         corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
#         varnames <- colnames(corData)
#         setnames(corData, c(1:2), c("V1","V2"))
#         # expressed <- which(apply(corData,1,sum)>0)
#         diffMat[(corData$V1+corData$V2)>0,x] = (log2(corData[(corData$V1+corData$V2)>0]$V1+1)-log2(corData[(corData$V1+corData$V2)>0]$V2+1))^2
#         
#     }
#     
#     #
#     if("tx_name" %in% colnames(dd)){
#         temp_maeValues <- dd[,.(tx_name, gene_name, ntx)]
#     }else{
#         temp_maeValues <- dd[,.(gene_name)]
#     }
#     if(expressionLevel) temp_maeValues[, expressed := apply(corData>expression_t,1,sum)>0]
#     temp_maeValues[,mae := apply(diffMat,1,mean,na.rm = TRUE)]
#     temp_maeValues[,`:=`(match_mae = NA, non_match_mae = NA)]
#     if(length(which(nameMat$rep_status))>1){
#         temp_maeValues[,match_mae := apply(diffMat[, which(nameMat$rep_status)],1,function(x) sqrt(mean(x, na.rm = TRUE)))]
#     }else if(length(which(nameMat$rep_status))==1){
#         temp_maeValues[,match_mae := sqrt(mean(diffMat[, which(nameMat$rep_status)],na.rm = TRUE))]
#     }
#     if(length(which(!nameMat$rep_status))>1){
#         temp_maeValues[,non_match_mae := apply(diffMat[, which(!nameMat$rep_status)],1,function(x) sqrt(mean(x, na.rm = TRUE)))]
#     }else if(length(which(!nameMat$rep_status))==1){
#         temp_maeValues[,match_mae := sqrt(mean(diffMat[, which(!nameMat$rep_status)],na.rm = TRUE))]
#     }
#     
#     temp_maeValues[, `:=`(
#         cellLine = c(cellLineList[[8]],"all")[t],
#         gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
#         # match_status = nameMat$rep_status,
#         common_type = post_fix,
#         protocol_comparison = paste(protocolV, collapse = " vs "))
#     ]
#     
#     maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
#     return(maeValues)
# }
# # calc_mae2 <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
# #     diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
# #     for( x in seq_len(nrow(nameMat))){
# #         corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
# #         varnames <- colnames(corData)
# #         setnames(corData, c(1:2), c("V1","V2"))
# #         # expressed <- which(apply(corData,1,sum)>0)
# #         diffMat[(corData$V1+corData$V2)>0,x] = abs(log2(corData[(corData$V1+corData$V2)>0]$V1+1)-log2(corData[(corData$V1+corData$V2)>0]$V2+1))
# #         
# #     }
# #     
# #     #
# #     if("tx_name" %in% colnames(dd)){
# #         temp_maeValues <- dd[,.(tx_name, gene_name, ntx)]
# #     }else{
# #         temp_maeValues <- dd[,.(gene_name)]
# #     }
# #     if(expressionLevel) temp_maeValues[, expressed := apply(corData>expression_t,1,sum)>0]
# #     temp_maeValues[,mae := apply(diffMat,1,mean,na.rm = TRUE)]
# #     temp_maeValues[,`:=`(match_mae = NA, non_match_mae = NA)]
# #     if(length(which(nameMat$rep_status))>1){
# #         temp_maeValues[,match_mae := apply(diffMat[, which(nameMat$rep_status)],1,mean,na.rm = TRUE)]
# #     }else if(length(which(nameMat$rep_status))==1){
# #         temp_maeValues[,match_mae := mean(diffMat[, which(nameMat$rep_status)],na.rm = TRUE)]
# #     }
# #     if(length(which(!nameMat$rep_status))>1){
# #         temp_maeValues[,non_match_mae := apply(diffMat[, which(!nameMat$rep_status)],1,mean,na.rm = TRUE)]
# #     }else if(length(which(!nameMat$rep_status))==1){
# #         temp_maeValues[,match_mae := mean(diffMat[, which(!nameMat$rep_status)],na.rm = TRUE)]
# #     }
# #     
# #     temp_maeValues[, `:=`(
# #         cellLine = c(cellLineList[[8]],"all")[t],
# #         gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
# #         # match_status = nameMat$rep_status,
# #         common_type = post_fix,
# #         protocol_comparison = paste(protocolV, collapse = " vs "))
# #     ]
# #     
# #     maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
# #     return(maeValues)
# # }
# 
# # mean absolute relative difference used by both salmon and kallisto: 2*abs(estimated-truth)/(estimated+truth)
# calc_mard <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
#     diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
#     for( x in seq_len(nrow(nameMat))){
#         corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
#         varnames <- colnames(corData)
#         setnames(corData, c(1:2), c("V1","V2"))
#         corData[, V1 := log2(V1+1)]
#         corData[, V2 := log2(V2+1)]
#         expressed <- which(apply(corData[,c("V1","V2"),with = TRUE],1,sum)>0)
#         diffMat[expressed,x] = 2*abs(corData$V1-corData$V2)[expressed]/(corData$V1+corData$V2)[expressed]
#         
#     }
#     
#     #
#     if("tx_name" %in% colnames(dd)){
#         temp_maeValues <- dd[,.(tx_name, gene_name, ntx)]
#     }else{
#         temp_maeValues <- dd[,.(gene_name)]
#     }
#     if(!is.null(dim(diffMat))){
#         if(expressionLevel) temp_maeValues[, expressed := apply(corData>expression_t,1,sum)>0]
#         temp_maeValues[,mae := apply(diffMat,1,mean,na.rm = TRUE)]
#     }else{
#         if(expressionLevel) temp_maeValues[, expressed := sum(corData>expression_t)>0]
#         temp_maeValues[,mae := mean(diffMat, na.rm = TRUE)]
#     }
#     
#     temp_maeValues[,`:=`(match_mae = NA, non_match_mae = NA)]
#     match_diffMat <- diffMat[, which(nameMat$rep_status)]
#     if(!is.null(dim(match_diffMat))){
#         temp_maeValues[,match_mae := apply(match_diffMat,1,mean,na.rm = TRUE)]
#     }else if(is.null(dim(match_diffMat))&(!isEmpty(match_diffMat))){
#         temp_maeValues[,match_mae := mean(match_diffMat,na.rm = TRUE)]
#     }
#     non_match_diffMat <- diffMat[, which(!nameMat$rep_status)]
#     if(!is.null(dim(non_match_diffMat))){
#         temp_maeValues[,non_match_mae := apply(non_match_diffMat,1,mean,na.rm = TRUE)]
#     }else if(is.null(dim(non_match_diffMat))&(!isEmpty(non_match_diffMat))){
#         temp_maeValues[,non_match_mae := mean(non_match_diffMat,na.rm = TRUE)]
#     }
#     
#     temp_maeValues[, `:=`(
#         cellLine = c(cellLineList[[8]],"all")[t],
#         gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
#         # match_status = nameMat$rep_status,
#         common_type = post_fix,
#         protocol_comparison = paste(protocolV, collapse = " vs "))
#     ]
#     
#     maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
#     return(maeValues)
# }
# 
# 
# calc_mard_mod <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
#     diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
#     for( x in seq_len(nrow(nameMat))){
#         corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
#         varnames <- colnames(corData)
#         setnames(corData, c(1:2), c("V1","V2"))
#         corData[, V1 := log2(V1+1)]
#         corData[, V2 := log2(V2+1)]
#         expressed <- which(apply(corData[,c("V1","V2"),with = TRUE],1,sum)>0)
#         diffMat[expressed,x] = 2*(corData$V1-corData$V2)[expressed]/(corData$V1+corData$V2)[expressed]
#         
#     }
#     
#     #
#     if("tx_name" %in% colnames(dd)){
#         temp_maeValues <- dd[,.(tx_name, gene_name, ntx)]
#     }else{
#         temp_maeValues <- dd[,.(gene_name)]
#     }
#     if(!is.null(dim(diffMat))){
#         if(expressionLevel) temp_maeValues[, expressed := apply(corData>expression_t,1,sum)>0]
#         temp_maeValues[,mae := apply(diffMat,1,mean,na.rm = TRUE)]
#     }else{
#         if(expressionLevel) temp_maeValues[, expressed := sum(corData>expression_t)>0]
#         temp_maeValues[,mae := mean(diffMat, na.rm = TRUE)]
#     }
#     
#     temp_maeValues[,`:=`(match_mae = NA, non_match_mae = NA)]
#     match_diffMat <- diffMat[, which(nameMat$rep_status)]
#     if(!is.null(dim(match_diffMat))){
#         temp_maeValues[,match_mae := apply(match_diffMat,1,mean,na.rm = TRUE)]
#     }else if(is.null(dim(match_diffMat))&(!isEmpty(match_diffMat))){
#         temp_maeValues[,match_mae := mean(match_diffMat,na.rm = TRUE)]
#     }
#     non_match_diffMat <- diffMat[, which(!nameMat$rep_status)]
#     if(!is.null(dim(non_match_diffMat))){
#         temp_maeValues[,non_match_mae := apply(non_match_diffMat,1,mean,na.rm = TRUE)]
#     }else if(is.null(dim(non_match_diffMat))&(!isEmpty(non_match_diffMat))){
#         temp_maeValues[,non_match_mae := mean(non_match_diffMat,na.rm = TRUE)]
#     }
#     
#     temp_maeValues[, `:=`(
#         cellLine = c(cellLineList[[8]],"all")[t],
#         gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
#         # match_status = nameMat$rep_status,
#         common_type = post_fix,
#         protocol_comparison = paste(protocolV, collapse = " vs "))
#     ]
#     
#     maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
#     return(maeValues)
# }
# 
# 
# 
# 
# calc_mard_ave <- function(t,g,nameMat,cellLineList, protocolV,geneClusterList,dd, post_fix, maeValues, expressionLevel, expression_t){
#     diffMat <- matrix(NA, nrow = nrow(dd), ncol = nrow(nameMat))
#     for( x in seq_len(nrow(nameMat))){
#         corData <- dd[,c(nameMat[x]$v1, nameMat[x]$v2),with = FALSE]
#         varnames <- colnames(corData)
#         setnames(corData, c(1:2), c("V1","V2"))
#         # expressed <- which(apply(corData,1,sum)>0)
#         diffMat[(corData$V1+corData$V2)>0,x] = abs(corData[(corData$V1+corData$V2)>0]$V1-corData[(corData$V1+corData$V2)>0]$V2)/(corData[(corData$V1+corData$V2)>0]$V1+corData[(corData$V1+corData$V2)>0]$V2)
#         
#     }
#     
#     #
#     if("tx_name" %in% colnames(dd)){
#         temp_maeValues <- dd[,.(tx_name, gene_name, ntx)]
#     }else{
#         temp_maeValues <- dd[,.(gene_name)]
#     }
#     if(expressionLevel) temp_maeValues[, expressed := apply(corData>expression_t,1,sum)>0]
#     temp_maeValues[,mae := apply(diffMat,1,mean,na.rm = TRUE)]
#     temp_maeValues[,`:=`(match_mae = NA, non_match_mae = NA)]
#     if(length(which(nameMat$rep_status))>1){
#         temp_maeValues[,match_mae := apply(diffMat[, which(nameMat$rep_status)],1,mean,na.rm = TRUE)]
#     }else if(length(which(nameMat$rep_status))==1){
#         temp_maeValues[,match_mae := mean(diffMat[, which(nameMat$rep_status)],na.rm = TRUE)]
#     }
#     if(length(which(!nameMat$rep_status))>1){
#         temp_maeValues[,non_match_mae := apply(diffMat[, which(!nameMat$rep_status)],1,mean,na.rm = TRUE)]
#     }else if(length(which(!nameMat$rep_status))==1){
#         temp_maeValues[,match_mae := mean(diffMat[, which(!nameMat$rep_status)],na.rm = TRUE)]
#     }
#     
#     temp_maeValues[, `:=`(
#         cellLine = c(cellLineList[[8]],"all")[t],
#         gene_cluster = c(geneClusterList[[length(geneClusterList)]],"all")[g],
#         # match_status = nameMat$rep_status,
#         common_type = post_fix,
#         protocol_comparison = paste(protocolV, collapse = " vs "))
#     ]
#     
#     maeValues <- do.call("rbind", list(maeValues, temp_maeValues))
#     return(maeValues)
# }
