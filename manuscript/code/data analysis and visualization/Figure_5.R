library(data.table)
library(ggplot2)
aa <- c(2702.5,
        1055.4,
        28701.8,
        2133.3,
        685.7,
        4085.4)

dt <- data.table(candidate_id = c(1:13),
                 lr_conc = c( 3405.2,
                              6619.2,
                              6506.8,
                              9762,
                              3183.6,
                              4928.4,
                              6544,
                              9949.6,
                              5726.4,
                              3812,
                              4531.2,
                              2565.6,
                              3560.4
                 ),
                 sr_conc = c( 1.26,
                              6.272,
                              4195.6,
                              6362,
                              3130.4,
                              5100,
                              0.228,
                              4.664,
                              0,
                              0,
                              6.608,
                              0.628,
                              0
                 ))
dt_long <- melt(dt, id.var = "candidate_id", measure.vars = c("sr_conc","lr_conc"))
dt_long[, variable := gsub("_conc","", variable)]

library(readxl)
dtEst <- data.table(as.data.frame(read_xlsx("digital_PCR/candidateTable_dPCR_updated_27Jun2024.xlsx", 
                                            col_names = TRUE)))
dtEst <- dtEst[!is.na(`Candidate Id`)]
dtEst[, candidate_id := as.integer(gsub("Candidate ","",`Candidate Id`))]

dt_long <- dtEst[,.(variable, candidate_id, lrEst, srEst)][dt_long, on = c("candidate_id","variable")]


setnames(dt_long, c("variable","value"), c("major_isoform_type","conc"))
dt_longlong <- melt(dt_long, id.vars = c("candidate_id","major_isoform_type","conc"),measure.vars = c("lrEst","srEst"))
setnames(dt_longlong, c("variable","value"), c("data_type","estimates"))

dt_longlong[, common_status := (candidate_id %in% c(3,4,5,6))]

dt_longlong[common_status == TRUE, 
            revised_estimates := sum(estimates), by = list(candidate_id, data_type)]
dt_longlong[is.na(revised_estimates), revised_estimates := estimates]
dt_longlong[common_status == TRUE&(major_isoform_type =="lr"), revised_estimates := estimates]
dt_longlong[, revised_isoform_type := ifelse(common_status == TRUE&(major_isoform_type=="sr"),"sr+lr",major_isoform_type)]
p2 <- ggplot(dt_longlong, aes(x = log10(conc+1), y = log10(revised_estimates+1)))+
    geom_abline(intercept = 0, slope = 1)+
    geom_text(aes(label = candidate_id, col = candidate_id %in% c(3,4,5,6)))+
    facet_wrap(common_status~data_type, scales = "free")+
    stat_cor(method = "spearman",
             label.x = 0,
             label.y = 4)+
    xlim(0,4)+
    ylim(0,4)+
    theme_classic()+
    theme(legend.position = "top")

## for correlation only combined a and b
ggplot(dt_longlong, aes(x = log10(conc+1), y = log10(revised_estimates+1)))+
    geom_abline(intercept = 0, slope = 1)+
    geom_text(aes(label = candidate_id, col = candidate_id %in% c(3,4,5,6)))+
    facet_wrap(~data_type, scales = "free")+
    stat_cor(method = "pearson",
             label.x = 0,
             label.y = 4)+
    xlim(0,4)+
    ylim(0,4)+
    theme_classic()+
    theme(legend.position = "top")

p1 <- ggplot(dt, aes(x = log10(lr_conc+1), y = log10(sr_conc+1)))+
    geom_abline(intercept = 0, slope = 1)+
    geom_text(aes(label = candidate_id, col = candidate_id %in% c(3,4,5,6)))+
    xlim(0,4)+
    ylim(0,4)+
    facet_wrap(~common_status,scales = "free", ncol = 1, nrow = 2)+
    stat_cor(method = "spearman",
             label.x = 0,
             label.y = 4)+
    theme_classic()+
    theme(legend.position = "top")

library(ggpubr)
pdf("dPCR_results_2Aug2024_text_spcor.pdf")
ggarrange(p1,p2, nrow=1, ncol =2, widths = c(1,2), align = "hv", common.legend = TRUE)
dev.off()


## plot the annotations pdf for two examples: annotations 
library(bambu)
library(dplyr)
customised_gtf <- "dpcr_candidateTx_ucsc.gtf"
customised_annotations <- prepareAnnotations(customised_gtf)
seq_ranges_sr_dt  <- bind_rows(list(
                                 data.frame(chr = "chr5",
                                 end = rev(c(82276255,82275295)),
                                 start = rev(c(82276178,82275244)),
                                 strand = "-",
                                 tx_name = dtEst[candidate_id==1&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr3",
                                            end = rev(c(12841545,12840353)),
                                            start = rev(c(12841494,12840142)),
                                            strand = "-",
                                            tx_name = dtEst[candidate_id==2&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr5",
                                            end = rev(c(40835186,40834606)),
                                            start = rev(c(40835019,40834471)),
                                            strand = "-",
                                            tx_name = dtEst[candidate_id==13&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr8",
                                            end = rev(c(73295789,73292797)),
                                            start = rev(c(73295765,73292689)),
                                            strand = "-",
                                            tx_name = dtEst[candidate_id==12&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr11",
                                            end = rev(c(61967631,61967453)),
                                            start = rev(c(61967548,61967312)),
                                            strand = "-",
                                            tx_name = dtEst[candidate_id==11&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr9",
                                            end = rev(c(19380198,19379618)),
                                            start = rev(c(19380088,19379487)),
                                            strand = "-",
                                            tx_name = dtEst[candidate_id==10&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr5",
                                            end = rev(c(150448913,150447735)),
                                            start = rev(c(150448833,150447585)),
                                            strand = "-",
                                            tx_name = dtEst[candidate_id==9&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr16",
                                            start = c(89560681,89561227),
                                            end = c(89560712,89561368),
                                            strand = "+",
                                            tx_name = dtEst[candidate_id==8&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr2",
                                            start = c(216499086,216499270),
                                            end = c(216499151,216499398),
                                            strand = "+",
                                            tx_name = dtEst[candidate_id==7&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr12",
                                            start = c(6534810,6536494,6536684),
                                            end = c(6534861,6536593,6536790),
                                            strand = "+",
                                            tx_name = dtEst[candidate_id==6&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr7",
                                            start = c(44799247,44799392,44799702),
                                            end = c(44799277,44799480,44799874),
                                            strand = "+",
                                            tx_name = dtEst[candidate_id==5&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr11",
                                            start = c(75402352,75404020,75404672),
                                            end = c(75402446,75404207,75404868),
                                            strand = "+",
                                            tx_name = dtEst[candidate_id==4&(variable == "sr")]$tx_name),
                                 data.frame(chr = "chr2",
                                            start = c(101002702,101004158,101005959),
                                            end = c(101002808,101004283,101006071),
                                            strand = "+",
                                            tx_name = dtEst[candidate_id==3&(variable == "sr")]$tx_name)))
 
seq_ranges_lr_dt  <- bind_rows(list(
                           data.frame(chr = "chr5",
                              end = c(82276255),
                                        start = c(82276061),
                                    strand = "-",
                                      tx_name = dtEst[candidate_id==1&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr3",
                                      end = rev(c(12841545,12840242)),
                                      start = rev(c(12841494,12840142)),
                                      strand = "-",
                                      tx_name = dtEst[candidate_id==2&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr5",
                                      end = rev(c(40835222,40834606)),
                                      start = rev(c(40835183,40834471)),
                                      strand = "-",
                                      tx_name = dtEst[candidate_id==13&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr8",
                                      end = rev(c(73293633,73292797)),
                                      start = rev(c(73293599,73292689)),
                                      strand = "-",
                                      tx_name = dtEst[candidate_id==12&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr11",
                                      end = rev(c(61967634)),
                                      start = rev(c(61967312)),
                                      strand = "-",
                                      tx_name = dtEst[candidate_id==11&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr9",
                                      end = rev(c(19380236,19379618)),
                                      start = rev(c(19380190,19379487)),
                                      strand = "-",
                                      tx_name = dtEst[candidate_id==10&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr5",
                                      end = rev(c(150449747,150447735)),
                                      start = rev(c(150449703,150447585)),
                                      strand = "-",
                                      tx_name = dtEst[candidate_id==9&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr16",
                                      start = c(89560657,89560940,89561227),
                                      end = c(89560712,89561063,89561368),
                                      strand = "+",
                                      tx_name = dtEst[candidate_id==8&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr2",
                                      start = c(216498844,216499270),
                                      end = c(216498877,216499398),
                                      strand = "+",
                                      tx_name = dtEst[candidate_id==7&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr12",
                                      start = c(6534517,6534810),
                                      end = c(6534569,6534861),
                                      strand = "+",
                                      tx_name = dtEst[candidate_id==6&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr7",
                                      start = c(44799702,44801287),
                                      end = c(44799874,44801630),
                                      strand = "+",
                                      tx_name = dtEst[candidate_id==5&(variable == "lr")]$tx_name),
                           data.frame(chr = "chr11",
                                      start = c(75404672,75405614),
                                      end = c(75404868,75405692),
                                      strand = "+",
                                      tx_name = dtEst[candidate_id==4&(variable == "lr")]$tx_name),
                        data.frame(chr = "chr2",
                                 start = c(101005959,101006350),
                                 end = c(101006071,101006421),
                            strand = "+",
                       tx_name = dtEst[candidate_id==3&(variable == "lr")]$tx_name)))
                                                              
sr_seq <- makeGRangesListFromDataFrame(seq_ranges_sr_dt, split.field = "tx_name",
                                       names.field = "tx_name")
lr_seq <- makeGRangesListFromDataFrame(seq_ranges_lr_dt, split.field = "tx_name",
                                       names.field = "tx_name")
library(ggbio)
pdf("dpcr_examples_annotations_allcandidates.pdf")
for(i in 1:13){
    sr_tx <- dtEst[candidate_id==i&(variable == "sr")]$tx_name
    lr_tx <- dtEst[candidate_id==i&(variable == "lr")]$tx_name
    sr_anno <- autoplot(customised_annotations[sr_tx], group.selfish = TRUE)
    lr_anno <- autoplot(customised_annotations[lr_tx], group.selfish = TRUE)
    srseq_anno <- autoplot(sr_seq[sr_tx], group.selfish = TRUE)
    lrseq_anno <- autoplot(lr_seq[lr_tx], group.selfish = TRUE)
    print(tracks(sr_anno,srseq_anno, lr_anno,lrseq_anno))
}
dev.off()



qPCR_c1to6 <- fread("20240719 Run 1 Cand1-6 Testing w RT-Tube 1 10xd -  Quantification Amplification Results_FAM.csv")
qPCR_c1to6_label <- fread("20240719 Run 1 Cand1-6 Testing w RT-Tube 1 10xd -  End Point Results_FAM.csv")

qPCR_c1to6_long <- melt(qPCR_c1to6, id.vars = c("Cycle"), measure.vars = colnames(qPCR_c1to6)[-c(1,2)])
setnames(qPCR_c1to6_long, "variable","Well")


qPCR_c1to6_label[, `:=`(candidate_id = unlist(strsplit(Sample, " "))[1],
                        short_read = grepl("S",unlist(strsplit(Sample, " "))[2])), by = Sample] 

qPCR_c1to6_label[substr(Well,2,2)==0, Well_new := gsub("0","",Well)]
qPCR_c1to6_label[!is.na(Well_new), Well := Well_new]

qPCR_c1to6_long <- qPCR_c1to6_label[,.(candidate_id, short_read, Well,`Sample Type`)][qPCR_c1to6_long, on = "Well"]



qPCR_c7to13 <- fread("20240719 Run 2 Cand7-13 Testing w RT-Tube 1 10xd -  Quantification Amplification Results_FAM.csv")
qPCR_c7to13_label <- fread("20240719 Run 2 Cand7-13 Testing w RT-Tube 1 10xd -  End Point Results_FAM.csv")
qPCR_c7to13_long <- melt(qPCR_c7to13, id.vars = c("Cycle"), measure.vars = colnames(qPCR_c7to13)[-c(1,2)])
setnames(qPCR_c7to13_long, "variable","Well")


qPCR_c7to13_label[, `:=`(candidate_id = unlist(strsplit(Sample, " "))[1],
                        short_read = grepl("S",unlist(strsplit(Sample, " "))[2])), by = Sample] 

qPCR_c7to13_label[substr(Well,2,2)==0, Well_new := gsub("0","",Well)]
qPCR_c7to13_label[!is.na(Well_new), Well := Well_new]

qPCR_c7to13_long <- qPCR_c7to13_label[,.(candidate_id, short_read, Well,`Sample Type`)][qPCR_c7to13_long, on = "Well"]


qPCR_results <- do.call("rbind",list(qPCR_c1to6_long,qPCR_c7to13_long))
qPCR_results[candidate_id == "C9S",`:=`(candidate_id = "C9", short_read = TRUE)]

pdf("qPCR_results_candidate_allcandidates.pdf", width = 10, height = 8)
ggplot(qPCR_results[`Sample Type` != "NTC"], aes(x = Cycle, y = value, group=Well))+
    geom_abline(intercept = 50, slope = 0)+
    geom_line(aes(col = short_read))+
    facet_wrap(~candidate_id, scales = "free")+
    xlab("Number of cycles")+
    ylab("Relative fluorescence units")+
    theme_classic()+
    theme(legend.position = "top")
dev.off()


qPCR_ct_part1 <- fread("20240719 Run 1 Cand1-6 Testing w RT-Tube 1 10xd -  Quantification Cq Results_0.csv")
qPCR_ct_part2 <- fread("20240719 Run 2 Cand7-13 Testing w RT-Tube 1 10xd -  Quantification Cq Results_0.csv")

qPCR_ct <- do.call("rbind",list(qPCR_ct_part1, qPCR_ct_part2))

qPCR_ct <- qPCR_ct[,.(Sample, Cq)]
qPCR_ct <- qPCR_ct[!is.na(Cq)&(grepl("tube",tolower(Sample)))]
qPCR_ct[, `:=`(candidate_id = unlist(strsplit(Sample, " "))[1],
               short_read = grepl("S",unlist(strsplit(Sample, " "))[2])), by = Sample]
qPCR_ct[, cq_mean := mean(Cq), by = Sample]
qPCR_ct_processed <- unique(qPCR_ct[,.(Sample, cq_mean, candidate_id, short_read)])
qPCR_ct_diff <- qPCR_ct_processed[order(candidate_id, -short_read), list(ct_diff = diff(cq_mean)), by = candidate_id]

qPCR_ct_diff[, lfc := 2^(-ct_diff)]
qPCR_ct_diff <- qPCR_ct_diff[-6]
qPCR_ct_diff[, common_status := (candidate_id %in% c("C3","C4","C5","C6"))]
pdf("qPCR_lfc_6Aug2024.pdf")

ggplot(qPCR_ct_diff, aes(x=common_status, y = log10(lfc)))+
    geom_abline(intercept = 0, slope = 0)+ 
    geom_boxplot()+scale_y_continuous(breaks = c(0,1:5))+
    theme_classic()
ggplot(qPCR_ct_diff, aes(x=common_status, y = ct_diff))+
    geom_abline(intercept = 0, slope = 0)+ 
    geom_boxplot()+
    theme_classic()
dev.off()


library(data.table)
library(ggplot2)
perc_values <-c(0.38,0.08,0.47,0.51,0.06,0.19,0.16,1,0.1,0.75,0.18,0.05,0.04,0.03, 0.02, 0.01,0)
xt <- data.table(x = perc_values,
           y = 1-perc_values)
xt[, id := 1:.N]
xt <- melt(xt, id.vars = "id", measure.vars = c("x","y"))

pdf("piecharts_for_validation_diagram.pdf")
ggplot(xt, aes(x = "", y = value, fill = factor(variable, levels =c("y","x"))))+
    geom_bar(width = 1, stat = "identity")+
    coord_polar("y",start = 0)+facet_wrap(~id)+
    theme_classic()
dev.off()

