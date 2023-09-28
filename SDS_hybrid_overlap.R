library(tidyverse)
# Read in SDS
gd <- read.delim('nterm-altTSS-counts-table-SDS.tsv')

# Read in gtf file with exon-level info on exon classification and hybrid status
# gtf contains pairs of transcripts (subsequent rows) from a gene where one transcript uses a hybrid exon as a first exon and one transcript uses a hybrid exon as an internal exon
# Information included along with standard gtf info for each pair includes:
### one, both, or neither are protein coding
### protein code for each transcript
### alignment score
hdf <- read.csv('HFE_DF_withInfo.txt')

# extract transcript names and create column in gtf
tn <- unlist(lapply(strsplit(unlist(lapply(strsplit(hdf$V9, split = ';'), "[[", 6)), split = " "), "[[", 3))
hdf$transcript_name <- tn

m_l <- c(paste(hdf$transcript_name[seq(1, (length(hdf$transcript_name)-1), by=2)], ';',  hdf$transcript_name[seq(2, (length(hdf$transcript_name)), by=2)], sep = ""),
         paste(hdf$transcript_name[seq(2, (length(hdf$transcript_name)), by=2)], ';',  hdf$transcript_name[seq(1, (length(hdf$transcript_name)-1), by=2)], sep = ""))
g_l <- c(paste(gd$anchor, ';', gd$other, sep = ""))
n.gd <- gd[which(g_l %in% m_l),]
overlap <- strsplit(g_l[g_l %in% m_l], split = ';')
annotated_overlap <- data.frame(transcript1 = unlist(lapply(overlap, "[[", 1)), transcript2 = unlist(lapply(overlap, "[[", 2)))
sds_o <- 1-(length(annotated_overlap$transcript1)/length(gd$anchor))
over <- length(annotated_overlap$transcript1)/length(gd$anchor)
(annotated_p1 <- ggplot(data.frame(g = c(1, 1), gt = c('SDS only', 'Annotated + SDS'), value = 
                                     c(over, sds_o)),
                        aes(fill = gt, y = value, x = g)) + geom_bar(stat="identity") + theme_bw() + xlab("") + ylab("Fraction of Number of SDS") +
    scale_fill_manual(values=c("navy", "yellow2", 'darkgrey', 'darkorange3'))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))