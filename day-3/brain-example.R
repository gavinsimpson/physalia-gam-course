# packages
pkgs <- c("ggplot2", "readr", "mgcv", "gratia", "dplyr", "patchwork")
vapply(pkgs, library, logical(1L), logical.return = TRUE,
       character.only = TRUE)

URL <- "https://bit.ly/gamair-brain"

brain <- read_csv(URL, col_types = "ddddd")
brain

# filter two outliers
brain <- brain |>
  filter(medFPQ > 5e-3)

# plot the response
brain |>
  ggplot(aes(x = medFPQ)) +
    geom_density()

# plot the response
brain |>
  ggplot(aes(x = medFPQ)) +
  geom_histogram()

# plot as a surface
brain |>
  ggplot(aes(x = X, y = Y, fill = medFPQ)) +
  geom_raster() + coord_equal() +
  scale_fill_viridis_c(option = "plasma")

# fit a Gaussian GAM
m_gaus <- gam(medFPQ ~ s(Y, X, k = 100), data = brain, method = "REML")

# plot the smooth
draw(m_gaus, dist = 0.03)

# model diagnostics
appraise(m_gaus, method = "simulate")

# fit a more appropriate GAM
m_gamma <- gam(medFPQ ~ s(Y, X, k = 100), data = brain, method = "REML",
               family = Gamma(link = "log"))

# plot the smooth
draw(m_gamma, dist = 0.03)

# model diagnostics
appraise(m_gamma, method = "simulate")

# Would an additive model be better?
m_add <- gam(medFPQ ~ s(Y, k = 30) + s(X, k = 30),
             data = brain, method = "REML",
             family = Gamma(link = "log"))

# draw model
draw(m_add, dist = 0.03)

# compare
m_gamma_ml <- update(m_gamma, . ~ ., method = "ML")
m_add_ml   <- update(m_add, . ~ ., method = "ML")

AIC(m_add_ml, m_gamma_ml)

# isotropic or anisotropic?
m_gamma_te <- gam(medFPQ ~ te(Y, X, k = 10), data = brain, method = "REML",
                  family = Gamma(link = "log"))

# compare
m_gamma_te_ml <- update(m_gamma_te, . ~ ., method = "ML")
AIC(m_gamma_ml, m_gamma_te_ml)

# ANOVA-like decomposition
m_gamma_ti <- gam(medFPQ ~ s(Y, k = 10, bs = "cr") + s(X, k = 10, bs = "cr") +
                    ti(Y, X, k = 10, bs = c("cr", "cr")),
                  data = brain, method = "REML",
                  family = Gamma(link = "log"))

# compare the fits visually
draw(m_gamma, dist = 0.03) +
  draw(m_gamma_te, dist = 0.03) +
  plot_layout(ncol = 2)

# Are the two halves of the brain symmetric / mirror images
# continuous by variable smooths

brain <- brain %>%
  mutate(Xc = X - 64.5,
         right = as.numeric(X < 64.5))

# !!!!! Note this is Xc
m_sym <- gam(medFPQ ~ s(Y, Xc, k = 100), data = brain, method = "REML",
               family = Gamma(link = "log"))
# !!!!! Note this is Xc
m_asym <- gam(medFPQ ~ s(Y, Xc, k = 100) + s(Y, Xc, k = 100, by = right),
              data = brain, method = "REML",
              family = Gamma(link = "log"))

anova(m_asym)

# this is less well justified
anova(m_sym, m_asym, test = "LRT")

# compare surfaces

# simulate some new data from m_gamma and peturb it
brain_new <- brain
mu <- fitted(m_gamma)
n <- length(mu)
ind <- brain_new$X < 60 & brain_new$Y < 20
mu[ind] <- mu[ind] / 3
set.seed(1)

brain_comb <- brain_new %>%
  mutate(medFPQ = rgamma(rep(1,n), mu / m_gamma$sig2, scale = m_gamma$sig2)) %>%
  bind_rows(brain) %>%
  mutate(sample_a = rep(c(0,1), each = n),
         sample_b = 1 - sample_a,
         sample = factor(rep(c("new", "orig"), each = n),
                         levels = c("orig", "new")))

# fit model assuming no difference
m_same <- gam(medFPQ ~ s(Y, X, k = 100),
              data = brain_comb,
              method = "REML",
              family = Gamma(link = "log"))

# model the difference
m_diff <- gam(medFPQ ~ s(Y, X, k = 100) + s(Y, X, by = sample_a, k = 100),
              data = brain_comb,
              method = "REML", family = Gamma(link = "log"))

anova(m_diff)

draw(m_diff, dist = 0.03)

# model both samples separately with a factor by smooth
m_fac_by <- gam(medFPQ ~ sample + s(Y, X, k = 100, by = sample),
              data = brain_comb,
              method = "REML",
              family = Gamma(link = "log"))

summary(m_fac_by)

draw(m_fac_by, dist = 0.03)

# compute the difference between the two surfaces
diffs <- difference_smooths(m_fac_by, smooth = "s(Y,X)")
diffs

# plot the difference
draw(diffs)

# ordered factor provides another way to do this
# first smooth is for the reference level of the factor and the
# by factor smooth(s) are smooth difference between another level
# reference level
brain_comb <- brain_comb %>%
  mutate(sample_o = ordered(sample))
contrasts(brain_comb$sample_o) <- "contr.treatment"

m_diff2 <- gam(medFPQ ~ sample_o + s(Y, X, k = 100) +
                s(Y, X, by = sample_o, k = 100),
              data = brain_comb,
              method = "REML", family = Gamma(link = "log"))
