#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
pdf('08_02_gsp_qual.pdf')
layout(1:length(argv))
for(i in argv){
  a <- read.table(i, sep = '\t', header = TRUE, row.names = 1)
  b <- barplot(t(a), beside = TRUE,
    col = c('lightblue', 'darkblue', 'red3', 'darkred'),
    ylab = 'Quality', main = i)
  rect(1, 0, ceiling(max(b)), 50, border = NA, col = rgb(1,1,1,1/2))
  arrows(x0=1, x1=ceiling(max(b)), y0=c(0,50), y1=c(0,50), lty = c(1,3),
    length = 0)
}
t <- dev.off()

