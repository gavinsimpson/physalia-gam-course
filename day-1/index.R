## ----setup, include=FALSE, cache=FALSE----------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(cache = TRUE, dev = 'svg', echo = TRUE, message = FALSE, warning = FALSE,
                      fig.height = 6, fig.width = 1.777777 * 6)

library('here')
library('mgcv')
library('gratia')
library('ggplot2')
library('purrr')
library('mvnfast')
library("tibble")
library('patchwork')
library('tidyr')
library("knitr")
library("viridis")
library('readr')
library('dplyr')
library('gganimate')

## plot defaults
theme_set(theme_minimal(base_size = 16, base_family = 'Fira Sans'))

## constants
anim_width <- 1000
anim_height <- anim_width / 1.77777777
anim_dev <- 'png'
anim_res <- 200


## -----------------------------------------------------------------------------
dbinom(x = 7, size = 10, prob = 0.7)


## ----binomial-pdf, echo = FALSE-----------------------------------------------
## Binomial Probability mass function
s <- seq(0, 40, by = 1)
n <- rep(c(20,20,40), each = length(s))
binom.pmf <- data.frame(x = rep(s, 3),
                        n = rep(c(20,20,40), each = length(s)),
                        p = rep(c(0.5, 0.7, 0.5), each = length(s)))
binom.pmf <- transform(binom.pmf,
                       pmf = dbinom(x, size = n, prob = p),
                       params = paste("n=", n, "; p=", p, sep = ""))

plt.binom <- ggplot(binom.pmf, aes(x = x, y = pmf, colour = params)) +
    geom_point() + labs(y = "Probability Mass", x = "No. of Successes")
plt.binom


## ----poisson-pdf, echo = FALSE------------------------------------------------
s <- seq(0, 20, by = 1)
poisson.pmf <- data.frame(x = rep(s, 3),
                          lambda = rep(c(1,4,10), each = length(s)))
poisson.pmf <- transform(poisson.pmf,
                         pmf = dpois(x, lambda = lambda),
                         params = paste("lambda=", lambda, sep = ""))

plt.poisson <- ggplot(poisson.pmf, aes(x = x, y = pmf, colour = params)) +
    geom_point() + labs(y = "Probability Mass", x = "Count")
plt.poisson


## ----load-darl-data, echo = TRUE----------------------------------------------
wasp <- read_csv(here("data", "darlingtonia.csv"), comment = "#",
                 col_types = "dl")


## ----fit-darlingtonia, echo = TRUE--------------------------------------------
m <- glm(visited ~ leafHeight, data = wasp, family = binomial)
m


## ----summary-darlingtonia, echo = TRUE----------------------------------------
summary(m)


## ----predict-darlingtonia, echo = TRUE, eval = FALSE--------------------------
## # data to predict at
## pdat <- with(wasp,
##              tibble(leafHeight = seq(min(leafHeight),
##                                      max(leafHeight),
##                                      length = 100)))
## # predict
## pred <- predict(m, pdat, type = "link", se.fit = TRUE)
## ilink <- family(m)$linkinv # g-1()
## pdat <- pdat %>%
##   bind_cols(data.frame(pred)) %>%
##   mutate(fitted = ilink(fit),
##          upper = ilink(fit + (2 * se.fit)),
##          lower = ilink(fit - (2 * se.fit)))
## # plot
## ggplot(wasp, aes(x = leafHeight,
##                  y = as.numeric(visited))) +
##     geom_point() +
##     geom_ribbon(aes(ymin = lower, ymax = upper,
##                     x = leafHeight), data = pdat,
##                 inherit.aes = FALSE, alpha = 0.2) +
##     geom_line(data = pdat, aes(y = fitted)) +
##     labs(x = "Leaf Height [cm]",
##          y = "Probability of visitation")


## ----predict-darlingtonia, eval = TRUE, echo = FALSE, fig.height = 6, fig.width = 6----
# data to predict at
pdat <- with(wasp,
             tibble(leafHeight = seq(min(leafHeight),
                                     max(leafHeight),
                                     length = 100)))
# predict
pred <- predict(m, pdat, type = "link", se.fit = TRUE)
ilink <- family(m)$linkinv # g-1()
pdat <- pdat %>%
  bind_cols(data.frame(pred)) %>%
  mutate(fitted = ilink(fit),
         upper = ilink(fit + (2 * se.fit)),
         lower = ilink(fit - (2 * se.fit)))
# plot
ggplot(wasp, aes(x = leafHeight,
                 y = as.numeric(visited))) +
    geom_point() +
    geom_ribbon(aes(ymin = lower, ymax = upper,
                    x = leafHeight), data = pdat,
                inherit.aes = FALSE, alpha = 0.2) +
    geom_line(data = pdat, aes(y = fitted)) +
    labs(x = "Leaf Height [cm]",
         y = "Probability of visitation")


## ----coeftab-darlingtonia, results = "asis", echo = FALSE---------------------
knitr::kable(round(summary(m)$coefficients, 4), format = "pipe")


## ----load-maddy, echo = TRUE--------------------------------------------------
maddy <- read_csv(here("data", "maddy-peat.csv"), col_types = "cdddddd")
maddy <- mutate(maddy, midDepth = upperDepth - (0.5 * abs(upperDepth - lowerDepth)),
                calMid = calUpper - (0.5 * abs(calUpper - calLower)))
maddy


## ----plot-maddy, echo = TRUE, eval = FALSE------------------------------------
## ggplot(maddy, aes(x = midDepth, y = calMid)) +
##     geom_point() +
##     labs(y = "Calibrated Age", x = "Depth")


## ----plot-maddy, echo = FALSE, eval = TRUE, fig.height = 6, fig.width = 6-----
ggplot(maddy, aes(x = midDepth, y = calMid)) +
    geom_point() +
    labs(y = "Calibrated Age", x = "Depth")


## ----peat-model, echo = TRUE--------------------------------------------------
m_gamma <- glm(calMid ~ midDepth, data = maddy, family = Gamma(link = "identity"))
summary(m_gamma)


## ----plot-maddy-fitted-gamma, eval=FALSE--------------------------------------
## # data to predict at
## pdat <- with(maddy,
##              tibble(midDepth = seq(min(midDepth),
##                                    max(midDepth),
##                                    length = 100)))
## # predict
## p_gamma <- predict(m_gamma, pdat, type = "link",
##                    se.fit = TRUE)
## ilink <- family(m_gamma)$linkinv
## # confidence interval
## p_gamma <- pdat %>%
##   bind_cols(data.frame(p_gamma)) %>%
##   mutate(fitted = ilink(fit),
##          upper = ilink(fit + (2 * se.fit)),
##          lower = ilink(fit - (2 * se.fit)))
## # plot
## p1 <- ggplot(maddy, aes(x = midDepth, y = calMid)) +
##     geom_ribbon(aes(ymin = lower, ymax = upper,
##                     x = midDepth), data = p_gamma,
##                 inherit.aes = FALSE, alpha = 0.2) +
##     geom_line(data = p_gamma, aes(y = fitted)) +
##     geom_point() +
##     labs(y = "Calibrated Age", x = "Depth",
##          title = "Gamma GLM")
## p1


## ----plot-maddy-fitted-gamma, echo = FALSE, eval = TRUE, fig.height = 6, fig.width = 6----
# data to predict at
pdat <- with(maddy,
             tibble(midDepth = seq(min(midDepth),
                                   max(midDepth),
                                   length = 100)))
# predict
p_gamma <- predict(m_gamma, pdat, type = "link",
                   se.fit = TRUE)
ilink <- family(m_gamma)$linkinv
# confidence interval
p_gamma <- pdat %>%
  bind_cols(data.frame(p_gamma)) %>%
  mutate(fitted = ilink(fit),
         upper = ilink(fit + (2 * se.fit)),
         lower = ilink(fit - (2 * se.fit)))
# plot
p1 <- ggplot(maddy, aes(x = midDepth, y = calMid)) +
    geom_ribbon(aes(ymin = lower, ymax = upper,
                    x = midDepth), data = p_gamma,
                inherit.aes = FALSE, alpha = 0.2) +
    geom_line(data = p_gamma, aes(y = fitted)) +
    geom_point() +
    labs(y = "Calibrated Age", x = "Depth",
         title = "Gamma GLM")
p1


## ----plot-maddy-fitted-gaussian, eval = FALSE---------------------------------
## # fit gaussian GLM
## m_gaus <- glm(calMid ~ midDepth, data = maddy,
##               family = gaussian)
## # predict
## p_gaus <- predict(m_gaus, pdat, type = "link",
##                   se.fit = TRUE)
## ilink <- family(m_gaus)$linkinv
## # prep confidence interval
## p_gaus <- pdat %>%
##   bind_cols(data.frame(p_gaus)) %>%
##   mutate(fitted = ilink(fit),
##          upper = ilink(fit + (2 * se.fit)),
##          lower = ilink(fit - (2 * se.fit)))
## # plot
## p2 <- ggplot(maddy, aes(x = midDepth, y = calMid)) +
##     geom_ribbon(aes(ymin = lower, ymax = upper,
##                     x = midDepth), data = p_gaus,
##                 inherit.aes = FALSE, alpha = 0.2) +
##     geom_line(data = p_gaus, aes(y = fitted)) +
##     geom_point() +
##     labs(y = "Calibrated Age",
##          x = "Depth",
##          title = "Linear Model")
## p2


## ----plot-maddy-fitted-gaussian, echo = FALSE, eval = TRUE, fig.height = 6, fig.width = 6----
# fit gaussian GLM
m_gaus <- glm(calMid ~ midDepth, data = maddy,
              family = gaussian)
# predict
p_gaus <- predict(m_gaus, pdat, type = "link",
                  se.fit = TRUE)
ilink <- family(m_gaus)$linkinv
# prep confidence interval
p_gaus <- pdat %>%
  bind_cols(data.frame(p_gaus)) %>%
  mutate(fitted = ilink(fit),
         upper = ilink(fit + (2 * se.fit)),
         lower = ilink(fit - (2 * se.fit)))
# plot
p2 <- ggplot(maddy, aes(x = midDepth, y = calMid)) +
    geom_ribbon(aes(ymin = lower, ymax = upper,
                    x = midDepth), data = p_gaus,
                inherit.aes = FALSE, alpha = 0.2) +
    geom_line(data = p_gaus, aes(y = fitted)) +
    geom_point() +
    labs(y = "Calibrated Age",
         x = "Depth",
         title = "Linear Model")
p2


## -----------------------------------------------------------------------------
library("patchwork")
p1 + p2


## ----hadcrut-temp-example, echo = FALSE---------------------------------------
URL <- "https://bit.ly/hadcrutv4"
# data are year, median of ensemble runs, certain quantiles in remaining cols
# take only cols 1 and 2
gtemp <- read_table(URL, col_types = 'nnnnnnnnnnnn', col_names = FALSE) %>%
    select(num_range('X', 1:2)) %>% setNames(nm = c('Year', 'Temperature'))

## Plot
gtemp_plt <- ggplot(gtemp, aes(x = Year, y = Temperature)) +
    geom_line() + 
    geom_point() +
    labs(x = 'Year', y = expression(Temeprature ~ degree*C))
gtemp_plt


## ----hadcrut-temp-example, echo = FALSE---------------------------------------
URL <- "https://bit.ly/hadcrutv4"
# data are year, median of ensemble runs, certain quantiles in remaining cols
# take only cols 1 and 2
gtemp <- read_table(URL, col_types = 'nnnnnnnnnnnn', col_names = FALSE) %>%
    select(num_range('X', 1:2)) %>% setNames(nm = c('Year', 'Temperature'))

## Plot
gtemp_plt <- ggplot(gtemp, aes(x = Year, y = Temperature)) +
    geom_line() + 
    geom_point() +
    labs(x = 'Year', y = expression(Temeprature ~ degree*C))
gtemp_plt


## ----hadcrut-temp-polynomial, echo = FALSE------------------------------------
p <- c(1,3,8,15)
N <- 300
newd <- with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N)))
polyFun <- function(i, data = data) {
    lm(Temperature ~ poly(Year, degree = i), data = data)
}
mods <- lapply(p, polyFun, data = gtemp)
pred <- vapply(mods, predict, numeric(N), newdata = newd)
colnames(pred) <- p
newd <- cbind(newd, pred)
polyDat <- gather(newd, Degree, Fitted, - Year)
polyDat <- mutate(polyDat, Degree = ordered(Degree, levels = p))
gtemp_plt + geom_line(data = polyDat, mapping = aes(x = Year, y = Fitted, colour = Degree),
                      size = 1.5, alpha = 0.9) +
    scale_color_brewer(name = "Degree", palette = "PuOr") +
    theme(legend.position = "right")


## ----read-hadcrut, echo = TRUE------------------------------------------------
library('readr')
library('dplyr')
URL <-  "https://bit.ly/hadcrutv4"
gtemp <- read_table(URL, col_types = 'nnnnnnnnnnnn', col_names = FALSE) %>%
    select(num_range('X', 1:2)) %>% setNames(nm = c('Year', 'Temperature'))


## ----show-hadcrut, echo = TRUE, dependson = -1--------------------------------
gtemp


## ----hadcrutemp-fitted-gam, echo = TRUE, results = 'hide'---------------------
library('mgcv')
m <- gam(Temperature ~ s(Year), data = gtemp, method = 'REML')
summary(m)


## ----hadcrutemp-fitted-gam, echo = FALSE--------------------------------------
library('mgcv')
m <- gam(Temperature ~ s(Year), data = gtemp, method = 'REML')
summary(m)


## ----hadcrtemp-plot-gam, echo = FALSE-----------------------------------------
N <- 300
newd <- as_tibble(with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N))))
pred <- as_tibble(as.data.frame(predict(m, newdata = newd, se.fit = TRUE,
                                        unconditional = TRUE)))
pred <- bind_cols(newd, pred) %>%
    mutate(upr = fit + 2 * se.fit, lwr = fit - 2*se.fit)

ggplot(gtemp, aes(x = Year, y = Temperature)) +
    geom_point() +
    geom_ribbon(data = pred,
                mapping = aes(ymin = lwr, ymax = upr, x = Year), alpha = 0.4, inherit.aes = FALSE,
                fill = "#fdb338") +
    geom_line(data = pred,
              mapping = aes(y = fit, x = Year), inherit.aes = FALSE, size = 1, colour = "#025196") +
    labs(x = 'Year', y = expression(Temeprature ~ degree*C))


## ----hadcrtemp-plot-gam, echo = FALSE-----------------------------------------
N <- 300
newd <- as_tibble(with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N))))
pred <- as_tibble(as.data.frame(predict(m, newdata = newd, se.fit = TRUE,
                                        unconditional = TRUE)))
pred <- bind_cols(newd, pred) %>%
    mutate(upr = fit + 2 * se.fit, lwr = fit - 2*se.fit)

ggplot(gtemp, aes(x = Year, y = Temperature)) +
    geom_point() +
    geom_ribbon(data = pred,
                mapping = aes(ymin = lwr, ymax = upr, x = Year), alpha = 0.4, inherit.aes = FALSE,
                fill = "#fdb338") +
    geom_line(data = pred,
              mapping = aes(y = fit, x = Year), inherit.aes = FALSE, size = 1, colour = "#025196") +
    labs(x = 'Year', y = expression(Temeprature ~ degree*C))

