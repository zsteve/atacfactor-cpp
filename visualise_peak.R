library(ggplot2)
library(Rcpp)

# much faster than equivalent R code
cppFunction('
    Rcpp::List kde_diff(NumericVector& x, int bw, NumericVector range){
      NumericVector out_y(512);
      NumericVector out_x(512);
      double stepsize = (range[1] - range[0])/512;
      double min = range[0];
      for(int i = 0; i < 512; i++){
        out_x[i] = min + stepsize*i;
        for(int j = 0; j < x.size(); ++j){
          double t = ((x[j] - out_x[i])/bw);
          out_y[i] += 1.00/(x.size() * bw)*1/sqrt(2*3.14159265359)*exp(-0.5*pow(t, 2));;
        }
      }

      return Rcpp::List::create(Rcpp::Named("x") = out_x, 
                                Rcpp::Named("y") = out_y);
    }
            
')
# 
# kde_diff <- function(x, bw, range){
#   output <- rep(0, 512)
#   output_x <- rep(0, 512)
#   stepsize <- (range$max - range$min)/512
#   min_x <- range$min
#   for(i in 1:512){
#     output_x[i] <- min_x + stepsize*i
#     for(j in 1:length(x)){
#       output[i] <- output[i] + 1/(length(x)*bw)*gaussian((x[j] - (min_x + stepsize*i))/bw)
#     }
#   }
#   return(list(x = output_x, y = output))
# }

peak_id = 2
bw = 12

t <- read.table('../output_fseq', header = T)
t_peak <- t[t$peakid == peak_id, ]

r <- read.table(sprintf("split_fseq/%d", peak_id))
names(r) <- c('id', 'chr', 'pos', 'strand', 'insert')


d_plus <- density(r[r$strand == '+', ]$pos, bw = bw)
d_minus <- density(r[r$strand == '-', ]$pos, bw = bw)
range <- c(min(d_plus$x, d_minus$x), max(d_plus$x, d_minus$x))
d_plus2 <- kde_diff(r[r$strand == '+', ]$pos, bw = bw, range)
d_minus2 <- kde_diff(r[r$strand == '-', ]$pos, bw = bw, range)

diff_y <- d_plus2$y - d_minus2$y

g <- ggplot() + geom_line(aes(x = d_plus2$x, y = diff_y), linetype = 'dashed') + 
  geom_line(aes(x = d_plus$x, y = d_plus$y), color = 'green', alpha = 0.5) + 
  geom_line(aes(x = d_minus$x, y = d_minus$y), color = 'red', alpha = 0.5) + 

  geom_vline(aes(xintercept = t_peak[, 'pos']), color = 'blue') + 
  ggtitle(sprintf("Peak id = %d; footprints = %d, centre = %d", t_peak$peakid[1], t_peak$n_fp[1], t_peak$centre[1])) + 
  geom_vline(aes(xintercept = t_peak$centre[1]), color = 'blue', linetype = "dotted")
g


# pileup of peaks
n_fp = 1

idx <- which(t$n_fp == n_fp)
t_peak_fp <- t[idx, ]
peak_ids <- unique(t_peak_fp$peakid)
indices <- match(peak_ids, t_peak_fp$peakid)
centres <- t_peak_fp$centre[indices]

sum_y_fw = rep(0, 512)
sum_y_rv = rep(0, 512)
for(i in 1:length(peak_ids)){
  r <- read.table(sprintf('split_fseq/%d', peak_ids[i]), header = F)
  nreads <- nrow(r)
  names(r) <- c('id', 'chr', 'pos', 'strand', 'insert')
  r$pos <- r$pos - centres[i]
  r_fw <- r[r$strand == '+', ]
  r_rv <- r[r$strand == '-', ]
  y_new_fw <- kde_diff(r_fw$pos, bw = bw, c(-300, 300))$y
  y_new_rv <- kde_diff(r_rv$pos, bw = bw, c(-300, 300))$y
  #print(ggplot() + geom_line(aes(x = seq(-300, 300, 601/512), y = y_new_fw + y_new_rv)) 
  #      + geom_line(aes(x = seq(-300, 300, 601/512), y=y_new_fw), color = 'green')
  #      + geom_line(aes(x = seq(-300, 300, 601/512), y=y_new_rv), color = 'red') + geom_vline(aes(xintercept = 0)))

  #readline(prompt="Press [enter] to continue")
  sum_y_fw <- sum_y_fw + 1/nreads*y_new_fw
  sum_y_rv <- sum_y_rv + 1/nreads*y_new_rv
  
  if(i %% 100 == 0){
    print(i)
  }
}


print(ggplot() + geom_line(aes(x = seq(-300, 300, 601/512), y = sum_y_fw + sum_y_rv)) 
      + geom_line(aes(x = seq(-300, 300, 601/512), y=sum_y_fw), color = 'green')
      + geom_line(aes(x = seq(-300, 300, 601/512), y=sum_y_rv), color = 'red') + geom_vline(aes(xintercept = 0)))

