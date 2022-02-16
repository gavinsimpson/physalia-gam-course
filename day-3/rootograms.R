# rootograms as a diagnostic plot

# packages
pkgs <- c("ggplot2", "here", "readr", "mgcv", "gratia", "dplyr", "patchwork",
          "tibble")
vapply(pkgs, library, logical(1L), logical.return = TRUE,
       character.only = TRUE)

# install.packages("coenocliner")
library("coenocliner")

## parameters for simulating
set.seed(1)
locs <- runif(100, min = 1, max = 10)     # environmental locations
A0 <- 90                                  # maximal abundance
mu <- 3                                   # position on gradient of optima
alpha <- 1.5                              # parameter of beta response
gamma <- 4                                # parameter of beta response
r <- 6                                    # range on gradient species is present
pars <- list(m = mu, r = r, alpha = alpha, gamma = gamma, A0 = A0)
nb.alpha <- 1.5                           # overdispersion parameter 1/theta
zprobs <- 0.3                             # prob(y == 0) in binomial model

pois <- coenocline(locs, responseModel = "beta", params = pars,
                   countModel = "poisson")
nb   <- coenocline(locs, responseModel = "beta", params = pars,
                   countModel = "negbin", countParams = list(alpha = nb.alpha))
zinb <- coenocline(locs, responseModel = "beta", params = pars,
                   countModel = "ZINB",
                   countParams = list(alpha = nb.alpha, zprobs = zprobs))

df <- setNames(cbind.data.frame(locs, pois, nb, zinb),
               c("x", "yPois", "yNegBin", "yZINB")) %>%
  as_tibble()

m_pois <- gam(yPois ~ s(x), data = df,
              family = poisson(), method = "REML")
m_nb   <- gam(yNegBin ~ s(x), data = df,
              family = poisson(), method = "REML")
m_zinb <- gam(yZINB ~ s(x), data = df,
              family = poisson(), method = "REML")

appraise(m_pois, method = "simulate")
appraise(m_nb, method = "simulate")
appraise(m_zinb, method = "simulate")

rg_pois <- rootogram(m_pois, max_count = 20)
draw(rg_pois) +
  labs(title = "Poisson data - Poisson fit")

rg_nb <- rootogram(m_nb, max_count = 20)
draw(rg_nb) +
  labs(title = "Negative binomial data - Poisson fit")

rg_zinb <- rootogram(m_zinb, max_count = 20)
draw(rg_zinb) +
  labs(title = "ZI Negative binomial data - Poisson fit")

m_nb_nb <- gam(yNegBin ~ s(x), data = df,
               family = nb(), method = "REML")


appraise(m_nb_nb, method = "simulate")

rg_nb_nb <- rootogram(m_nb_nb, max_count = 20)
draw(rg_nb_nb) +
  labs(title = "Negative binomial data - Negative binomial fit")
