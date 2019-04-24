filtstat <- function(df, make_plot = FALSE, filter_len = 20){
  df$year <- round(df$time)
  splt <- split(df, df$year)
  yrd <- sapply(splt, function(x) max(x$reports))
  phlf <- sapply(splt, function(x) mean(x$half == "second"))
  yf <- stats::filter(yrd, filter = filter_len:1, sides = 1,
                      method = "convolution")
  N <- length(yf)
  D <- data.frame(cases = yrd[-1], filt = yf[-N], camp = phlf[-1])
  if (make_plot) {
    g <- ggplot(data = D, aes(x = filt, y = cases, col = camp)) +
      geom_point() + scale_y_continuous(trans = "log1p")
    print(g)
  }
  m <- lm(log(cases + 1)~filt + filt:camp, data = D)
  coef(m)["filt:camp"]
}
