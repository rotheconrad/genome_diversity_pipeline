#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
pdf('08_02_gsp_qual.pdf')
layout(1:length(argv))
for(i in argv){
  a <- read.table(i, sep = '\t', header = TRUE, row.names = 1)
  b <- barplot(t(a), beside = TRUE,
    col = c('lightblue', 'darkblue', 'red3', 'darkred'),
    ylab = 'Quality', main = i)
  right.end <- ceiling(max(b))
  rect(1, 0, right.end, 50, border = NA, col = rgb(1,1,1,3/4))
  arrows(x0 = 1, x1 = right.end, y0 = c(0,50), y1 = c(0,50), lty = c(1,3),
    length = 0)
  ablines(v = seq(1, right.end, by = nrow(b)+1), col = grey(1/2))
}
t <- dev.off()

