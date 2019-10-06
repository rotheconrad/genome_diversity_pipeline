#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
pdf('08_02_gsp_qual.pdf')
layout(1:length(argv))
for(i in argv){
  a <- read.table(i, sep = '\t', header = TRUE, row.names = 1)
  barplot(t(a), beside = TRUE,
    col = c('lightblue', 'darkblue', 'red3', 'darkred'),
    ylab = 'Quality', main = i)
  abline(h = 50, lty = 3)
}
dev.off()

