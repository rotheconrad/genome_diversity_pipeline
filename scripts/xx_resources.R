#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
a <- read.table(args[1], sep = '\t', header = TRUE)

plot_res <- function(x, y, k, ...) {
  step_names <- c('trim', 'norm', 'asm', 'maxbin', 'metabat', 'derep')
  col <- rainbow(max(k)-1, s = 3/4, v = 3/4)
  sel <- !is.na(x) & !is.na(y) & !is.na(k)
  plot(x, y, col = col[k-1], pch = 19, cex = 2, las = 1,
    xlab = 'Requested', ylab = 'Used', log = 'xy', ...)
  abline(0, 1, lty = 2)
  legend('topleft', legend = paste(2:max(k), step_names), col = col, pch = 19)
}

pdf(args[2])
plot_res(a$ram_req/1e3, a$ram_use/1e3, a$step, main = 'RAM (Gb)')
plot_res(a$wall_req/24, a$wall_use/24, a$step, main = 'Wall time (h)')
plot_res(a$cpu_req/24, a$cpu_use/24, a$step, main = 'CPU time (h)')
ttt <- dev.off()

