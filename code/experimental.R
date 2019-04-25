


iepc <- function(its) {
  mlevs <- sort(unique(its))
  ieps <- numeric(length(mlevs))
  for(i in seq_along(mlevs)){
    start <- min(which(its >= mlevs[i]))
    stop <- max(which(its >= mlevs[i]))
    enc <- rle(its[start:stop] >= mlevs[i])
    sel <- enc$values == FALSE
    ieps[i] <- mean(enc$lengths[sel])
  }
  cbind(lev = mlevs, miep = ieps)
}

d1 <- iepc(yrd[30:50])
d2 <- iepc(yrd[-c(1:50)])

plot(log(d1))
points(log(d2), col = 2)

abline(lm(log(miep)~log(lev), data = data.frame(d1)))
abline(lm(log(miep)~log(lev), data = data.frame(d2)), col = 2)

plot((d1))
points((d2), col = 2)

abline(lm((miep)~(lev), data = data.frame(d1)))
abline(lm((miep)~(lev), data = data.frame(d2)), col = 2)

x <- 1:10


stats::filter(x, filter = c(3, 2, 1), sides = 1, method = "convolution")



discount_data %>% filter(sim == 17) -> foo

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
      geom_point() + scale_y_log10()
    print(g)
  }
  m <- lm(log(cases)~filt + filt:camp, data = D)
  coef(m)["filt:camp"]
}

simdf <- split(discount_data, discount_data$sim)
stats <- sapply(simdf, filtstat)
mean(stats < 0)

preds <- rbind(cbind(stats, 1), cbind(rep(0, 500), 0))

filter(x, rep(1, 3), sides = 1)
filter(x, rep(1, 3), sides = 1, circular = TRUE)



filtstat <- function(df, make_plot = FALSE, filter_len = 20){
  df$year <- round(df$time)
  splt <- split(df, df$year)
  yrd <- sapply(splt, function(x) max(x$reports))
  phlf <- sapply(splt, function(x) mean(x$half == "second"))
  yf <- stats::filter(yrd, filter = filter_len:1, sides = 1,
                      method = "convolution")
  N <- length(yf)
  D <- data.frame(cases = yrd[-1],
                  filt = scale(yf[-N], center = FALSE),
                  camp = phlf[-1])
  tmp <- cumsum(D$camp)
  D$ccamp <- ifelse(tmp > filter_len, filter_len, tmp)
  #D <- D[D$cases > 100, ]
  if (make_plot) {
    g <- ggplot(data = D, aes(x = filt, y = cases, col = ccamp)) +
      geom_point() + scale_y_continuous(trans = "log1p")
    print(g)
  }

  obj <- function(par) {
    logsz <- par[1]
    logint <- par[2]
    logs1 <- par[3]
    logs2 <- par[4]
    mm <- model.matrix(~filt + filt:camp, data = D)
    mf <- model.frame(cases~filt + filt:camp, data = D)
    eta <- exp(mm %*% c(int, s1, logs2))
    -sum(dnbinom(mf$cases, size = exp(logsz), mu = eta, log = TRUE))
  }
  #init <- c(0, log(mean(D$cases)), 0, 0)

  #browser()
  #ans <- optim(init, obj)
  #ans <- MASS::glm.nb(cases ~ filt + filt:ccamp, link = log, data = D)
  ans <- lm(log(cases + 1) ~ filt + filt:ccamp, data = D)
  #m <- lm(log(cases + 1)~filt + filt:camp, data = D)
  #coef(m)["filt:camp"]
  coef(ans)["filt:ccamp"]
}


filtrecs <- list()
nrecs <- 0
for(do_city in unique(data_for_ews$city)){
  for(do_speed in unique(data_for_ews$vacc_speed)){

    sel_data <- data_for_ews %>% filter(city == do_city) %>%
      filter(vacc_speed == do_speed) %>%
      dplyr::select(time, reports, sim, half)

    simdf <- split(sel_data, sel_data$sim)
    stats <- sapply(simdf, filtstat)
    auc <- mean(stats < 0)
    nrecs <- nrecs + 1
    filtrecs[[nrecs]] <- tibble(vacc_speed = do_speed,
                                    AUC = auc, metric = "self-inhibition",
                                    city = do_city)
  }
}
(filttab <- bind_rows(filtrecs) %>% arrange(city, vacc_speed))

ggplot(filttab, aes(x = vacc_speed, y = AUC, color = city)) + geom_line() + geom_point(aes(shape = city)) + ylim(c(0.5,1)) + theme_classic()

