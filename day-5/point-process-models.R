# Analyse the Global temperature anomaly record with GAM and GAMM
pkgs <- c("mgcv", "spatstat", "ggplot2", "readr", "dplyr", "gratia")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

# load the gorilla data
data(gorillas, package = "spatstat.data")
# loads gorillas and gorillas.extra

# create the point process and quadrature values and set up for a PP as GLM
pp <- data.frame(gorillas,
  lapply(gorillas.extra, function(x) {
    x[gorillas]
  }), pt = 1, wt = 1e-6)

q_xy <- data.frame(gorillas.extra[[1]])[, c("x", "y")] # extract x and y from window

quad <- data.frame(q_xy,
  lapply(gorillas.extra, function(x) {
    x[q_xy]
  }),
  pt = 0,
  wt = area(gorillas$window) / nrow(q_xy))

dat <- merge(pp, quad, all = TRUE, sort = FALSE)

# center and scale covariates
dat <- dat |>
  mutate(
    elevation = scale(elevation),
    slopeangle = scale(slopeangle),
    waterdist = scale(waterdist)
  )

ctrl <- gam.control(trace = FALSE, nthreads = 4)
tic()
m_k200 <- gam(pt / wt ~ elevation + waterdist + slopeangle + heat +
  slopetype + vegetation + s(x, y, bs = "gp", k = 200),
data =  dat, family = poisson(), weights = wt, method = "REML",
control = ctrl)
toc()

tic()
m_k400 <- gam(pt / wt ~ elevation + waterdist + slopeangle + heat +
  slopetype + vegetation + s(x, y, bs = "gp", k = 400),
data = dat, family = poisson(), weights = wt, method = "REML",
control = ctrl)
toc()

# might be broken - used plot.gam(m_k200, all.terms = TRUE)
draw(m_k200, parametric = TRUE, rug = FALSE)

# might be broken - used plot.gam(m_k400, all.terms = TRUE)
draw(m_k400, parametric = TRUE, rug = FALSE)

sum(m_k200$edf) / m_k200$smooth[[1]]$bs.dim
sum(m_k400$edf) / m_k400$smooth[[1]]$bs.dim

# determine the range of between-point distances
dists <- fields::rdist(dat[dat$pt == 1, c("x", "y")])
range_interval <- range(dists[dists != 0])

# set up the function to be minimized
objective_fn <- function(rho) {
  tmp.m <- gam(
    pt / wt ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation +
      s(x, y, bs = "gp", k = 400, m = c(3, rho)),
    data = dat, family = poisson(), weights = wt, method = "REML")
  return(tmp.m$gcv.ubre) # the "method" specific criterion
}

# find the optimized range parameter - super slow
### opt <- optimize(objective_fn, interval = range_interval)
### opt$minimum
opt_min <- 568.5373

m <- gam(pt / wt ~ elevation + waterdist + slopeangle + heat + slopetype +
  vegetation + s(x, y, bs = "gp", k = 400, m = c(3, opt_min)),
data = dat, family = poisson(), weights = wt, method = "REML")

# fit an IPP model to contrast with the fitted LGCP
m_ipp <- gam(pt / wt ~ elevation + waterdist + slopeangle + heat + slopetype +
  vegetation,
data = dat, family = poisson(), weights = wt, method = "REML")

# set the domain data points (in this case the quadrature we used)
domain.grid <- dat[dat$pt == 0, ]

# predict intensity values
domain.grid$z_ipp <- predict(m_ipp, newdata = domain.grid, type = "response")

domain.grid$z <- predict(m, newdata = domain.grid, type = "response")

# create the pixel images (uses the window supplied in the original gorillas data)
pred_ipp.im <- as.im(domain.grid[, c("x", "y", "z_ipp")], W = gorillas$window)

pred.im <- as.im(domain.grid[, c("x", "y", "z")], W = gorillas$window)

# calculate the observed K functions
K_obs_ipp <- Kinhom(gorillas, lambda = pred_ipp.im, correction = "border")
K_obs <- Kinhom(gorillas, lambda = pred.im, correction = "border")

# simulate the envelopes/bounds
K_env_ipp <- envelope(gorillas, fun = Kinhom,
  simulate = expression(rpoispp(lambda = pred_ipp.im)))
K_env <- envelope(gorillas, fun = Kinhom,
  simulate = expression(rpoispp(lambda = pred.im)))

# plotting
layout(mat = matrix(1:2, nrow = 2, ncol = 1, byrow = TRUE), widths = 1,
  heights = c(0.5, 0.5))
op <- par(mar = c(2.1, 3.1, 2.1, 0))
plot(K_env_ipp$r, K_env_ipp$mmean, type = "n",
  ylim = range(c(K_env_ipp$obs, K_env_ipp$hi, K_env_ipp$lo)), ylab = "",
  xlab = "", xaxt = "n", yaxt = "n")
axis(side = 2, at = seq(0, 3e6, by = 1e6),
  labels = c("0", seq(1e6, 3e6, by = 1e6)))
polygon(c(rev(K_env_ipp$r), K_env_ipp$r), c(rev(K_env_ipp$hi), K_env_ipp$lo),
  col = "grey80", border = NA)
lines(K_env_ipp$r, K_env_ipp$mmean, lty = "dashed")
lines(K_obs_ipp$r, K_obs_ipp$border, col = "red")
mtext(text = "A: Poisson Process", side = 3, cex = 1, line = 0.5, adj = 0)
par(mar = c(3.1, 3.1, 1.1, 0))
plot(K_env$r, K_env$mmean, type = "n",
  ylim = range(c(K_env$obs, K_env$hi, K_env$lo)), ylab = "", xlab = "",
  yaxt = "n")
axis(side = 2, at = seq(0, 3e6, by = 1e6),
  labels = c("0", seq(1e6, 3e6, by = 1e6)))
polygon(c(rev(K_env$r), K_env$r), c(rev(K_env$hi), K_env$lo), col = "grey80",
  border = NA)
lines(K_env$r, K_env$mmean, lty = "dashed")
lines(K_obs$r, K_obs$border, col = "red")
mtext(text = "distance (m)", side = 1, srt = 90, cex = 1, line = 2)
mtext(text = "Inhomogeneous K Function", side = 2, srt = 90, cex = 1,
  xpd = TRUE, outer = T, line = -1)
mtext(text = "B: log-Gaussian Cox Process", side = 3, cex = 1, line = 0.5,
  adj = 0)
legend(x = 0, y = 4e6,
legend = c("Observed", "Theoretic", "95% Sim. Bounds"),
  col = c("red", "black", "grey80"),
  lty = c("solid", "dashed", "solid"), cex = 1,
  lwd = c(2, 2, 2), bty = "n")
par(op)
layout(1)

# fit an IPP for comparison
m_ipp_poly <- gam(pt / wt ~ poly(elevation, 2) + poly(waterdist, 2) +
  poly(slopeangle, 2) + heat + slopetype + vegetation,
data = dat, family = poisson(), weights = wt, method = "REML")

# fit the LGCP
m_poly <- gam(pt / wt ~ poly(elevation, 2) + poly(waterdist, 2) +
  poly(slopeangle, 2) + heat + slopetype + vegetation +
  s(x, y, bs = "gp", k = 400, m = c(3, opt_min)),
data = dat, family = poisson(), weights = wt, method = "REML")

# fit an IPP for comparison
m_ipp_sm <- gam(pt / wt ~ s(elevation) + s(waterdist) + s(slopeangle) + heat +
  slopetype + vegetation, data = dat, family = poisson(), weights = wt,
method = "REML")

# fit the LGCP
m_sm <- gam(pt / wt ~ s(elevation) + s(waterdist) + s(slopeangle) + heat +
  slopetype + vegetation + s(x, y, bs = "gp", k = 400, m = c(3, opt_min)),
data = dat, family = poisson(), weights = wt, method = "REML")

info_crit <- function(x){c(‘GAM criterion‘=ifelse(is.null(x$gcv.ubre),NA,x$gcv.ubre),
logLik = logLik(x), AIC = AIC(x), BIC = AIC(x, k = log(sum(dat$pt))))}

data.frame(‘IPP Linear‘=info_crit(m_ipp), ‘IPP Poly‘=info_crit(m_ipp_poly),
‘IPP Smoooth‘=info_crit(m_ipp_sm), ‘LGCP Linear‘ = info_crit(m),
‘LGCP Poly‘ = info_crit(m_poly), ‘LGCP Smoooth‘ = info_crit(m_sm))
