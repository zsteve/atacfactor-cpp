library(ggplot2)

peak_id = 30
bw = 12

t <- read.table('../output', header = T)
t_peak <- t[t$peakid == peak_id, ]

r <- read.table(sprintf("split/%d", peak_id))
names(r) <- c('id', 'chr', 'pos', 'strand', 'insert')

d <- density(r$pos, bw = bw)
d_plus <- density(r[r$strand == '+', ]$pos, bw = bw)
d_minus <- density(r[r$strand == '-', ]$pos, bw = bw)

attach(d)

g <- ggplot() + geom_line(aes(x = x, y = y), linetype = 'dashed') + 
  geom_line(aes(x = d_plus$x, y = d_plus$y), color = 'green') + 
  geom_line(aes(x = d_minus$x, y = d_minus$y), color = 'red') + 
  geom_vline(aes(xintercept = t_peak[, 'pos']), color = 'blue') + 
  ggtitle(sprintf("Peak id = %d; footprints = %d, centre = %d", t_peak$peakid[1], t_peak$n_fp[1], t_peak$centre[1])) + 
  geom_vline(aes(xintercept = t_peak$centre[1]), color = 'blue', linetype = "dotted")
g

detach(d)
