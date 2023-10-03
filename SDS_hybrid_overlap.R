library(tidyverse)
# Read in SDS
gd <- read.delim('nterm-altTSS-counts-table-SDS.tsv')


# HFE_DF.csv
# File details every possible annotated hybrid exon usage swap from first to internal

# Pairs of rows contain:
# (top row): details about the transcript which uses the exon as a first followed by 
# (bottom row): details about the transcript which uses the exon as an internal

# Only the transcript names from the pairs are used for this file, looking for overlap with a set of SDS 
hdf <- read.csv('HFE_DF.csv')


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
