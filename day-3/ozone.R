# Analyze the Ozone data set from LA

#' Ozone in LA in 1976
#'
#' A study the relationship between atmospheric ozone concentration and
#' meteorology in the Los Angeles Basin in 1976.  A number of cases with
#' missing variables have been removed for simplicity.
#'
#'
#' @name ozone
#' @docType data
#' @format A data frame with 330 observations on the following 10 variables.
#' \describe{ \item{O3}{Ozone conc., ppm, at Sandbug AFB.}
#' \item{vh}{a numeric vector} \item{wind}{wind speed}
#' \item{humidity}{a numeric vector} \item{temp}{temperature}
#' \item{ibh}{inversion base height} \item{dpg}{Daggett
#' pressure gradient} \item{ibt}{a numeric vector}
#' \item{vis}{visibility} \item{doy}{day of the year} }
#' @source Breiman, L. and J. H. Friedman (1985). Estimating optimal
#' transformations for multiple regression and correlation. Journal of the
#' American Statistical Association 80, 580-598.
#' @keywords datasets
#' @examples
#'

# load packages
pkgs <- c("mgcv", "here", "readr", "ggplot2", "gratia", "dplyr", "patchwork")
vapply(pkgs, library, logical(1L), character.only = TRUE,
    logical.return = TRUE)

# load the data
ozone <- read_csv("https://bit.ly/gam-ozone-data")
# or with here
#ozone <- read_csv(here("data", "ozone.csv"))

# analysis

# plot the data
o3_lab <- expression(O[3])
temp_lab <- expression(Temperature ~ (degree * F))

ozone |> # %>%
    ggplot(aes(y = O3, x = temp)) +
    geom_point() +
    labs(y = o3_lab, x = temp_lab)

ozone |>
    ggplot(aes(y = O3, x = ibh)) +
    geom_point() +
    labs(y = o3_lab, x = "Inversion base height")

ozone |>
    ggplot(aes(y = O3, x = ibt)) +
    geom_point() +
    labs(y = o3_lab, x = "Inversion base temperature")

# Start with a simple model
# this is basically lm()
lm1 <- gam(O3 ~ temp + ibh + ibt, data = ozone, method = "ML")
summary(lm1)

m1 <- gam(O3 ~ s(temp) + s(ibh) + s(ibt),
    data = ozone, method = "ML")
summary(m1)
draw(m1)

# Is the relationship between ozone and temperature linear?
# We test for this by adding a linear term, but then we have to remove the
# linear basis function from the TPRS of temperature. That is done by specifying
# `m = c(2,0)`, which for a TPRS means we want the usual second order derivative
# penalty but zero (0) null space
b1 <- basis(s(temp), data = ozone)
b0 <- basis(s(temp, m = c(2, 0)), data = ozone)
wrap_plots(draw(b1) + labs(subtitle = "With a null space"),
    draw(b0) + labs(subtitle = "Without a null space"), ncol = 2)

# Is the smooth of temperature adding anything?
m2 <- gam(O3 ~ temp + s(temp, m = c(2, 0)) + s(ibh) + s(ibt),
    data = ozone, method = "ML")
summary(m2)
# seems so, but...
draw(m2)

# ... the credible interval on the smooth is wide, though importantly it does
# not include zero *everywhere*
# We can look at draws from the posterior of the smooth using
samp_s_temp <- smooth_samples(m2, term = "s(temp)", n = 50)

# ww also evaluate the smooth
sm_s_temp <- smooth_estimates(m2, "s(temp)")

# and plot it with the 50 posterior draws
draw(sm_s_temp) +
    geom_line(data = samp_s_temp,
        aes(y = value, x = .x1, group = draw),
    alpha = 0.3, colour = "steelblue", size = 1)

# diagnostics
appraise(m1, method = "simulate")
# ah! non-constant variance

# refit as Tweedie?
m3 <- update(m1, . ~ ., family = tw())
summary(m3)
draw(m3)
appraise(m3)
# that looks better; the bands in the residuals are because of the slight
# integer nature of the data. Alternatively we might (ab)use a count model
# and use the negative binomial distribution for this

# negative binomial model
m4 <- gam(O3 ~ s(temp) + s(ibh) + s(ibt),
    data = ozone, method = "ML", family = nb())
summary(m4)
draw(m4)
appraise(m4, method = "simulate")

# quasipoisson model
m5 <- gam(O3 ~ s(temp) + s(ibh) + s(ibt),
    data = ozone, method = "ML", family = quasipoisson())
summary(m5)
draw(m5)
appraise(m5, method = "simulate")

# the Tweedie has perhaps the best justification so we'll proceed with that
# we could perform model selection using `select = TRUE`
m6 <- gam(O3 ~ s(temp) + s(ibh) + s(ibt),
    data = ozone, method = "ML", family = tw(), select = TRUE)
summary(m6)
draw(m6)
# we see that f(ibt) has been essentially shrunken out of the model

# a full model
m_full <- gam(O3 ~ s(temp) + s(ibh) + s(ibt) + s(humidity) + s(vh) +
        s(wind) + s(dpg) + s(vis) + s(doy),
    data = ozone, method = "REML", family = tw(), select = TRUE)

summary(m_full)
draw(m_full)
appraise(m_full, method = "simulate", n_simulate = 1000)

conc <- concrvity(m_full, pairwise = TRUE)
draw(conc)
