# Analyse the Global temperature anomaly record with GAM and GAMM
pkgs <- c("mgcv", "ggplot2", "readr", "dplyr", "gratia")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

# load temperature record
URL <- "https://bit.ly/hadcrutv4"
# data are year, median of ensemble runs, certain quantiles in remaining cols
# take only cols 1 and 2
gtemp <- read_table(URL, col_types = 'nnnnnnnnnnnn',
                    col_names = FALSE) %>%
    select(num_range('X', 1:2)) %>%
    setNames(nm = c('Year', 'Temperature'))

# Plot
gtemp_plt <- ggplot(gtemp, aes(x = Year, y = Temperature)) +
    #geom_line() +
    geom_point() +
    labs(x = 'Year', y = expression(Temeprature ~ degree * C))
gtemp_plt

# model the global temperature anomaly as a smooth plus AR(1) residuals
m_ar1 <- gamm(Temperature ~ s(Year, k = 20),
              data = gtemp,
              method = "REML",
              correlation = corAR1(form = ~ 1))

# model summary
summary(m_ar1$gam)

# plot
draw(m_ar1$gam)

## is this a better fit over the none AR(1) model?
m_0 <- gamm(Temperature ~ s(Year, k = 20),
            data = gtemp,
            method = "REML")

## GLRT
anova(m_0$lme, m_ar1$lme)

intervals(m_ar1$lme)$corStruct

compare_smooths(m_0, m_ar1, smooths = "s(Year)") %>% draw()

acf(resid(m_0$lme))

## yes!

## what about AR(p) for p > 1?
## For that we need corARMA and we have to refit the AR(1) just to be very sure
## everything is comparable
m_arma_p1 <- gamm(Temperature ~ s(Year, k = 20),
                  data = gtemp,
                  method = "REML",
                  correlation = corARMA(form = ~ 1, p = 1))
m_arma_p2 <- gamm(Temperature ~ s(Year, k = 20),
                  data = gtemp,
                  method = "REML",
                  correlation = corARMA(form = ~ 1, p = 2))
m_arma_p3 <- gamm(Temperature ~ s(Year, k = 20),
                  data = gtemp,
                  method = "REML",
                  correlation = corARMA(form = ~ 1, p = 3))

## GLRT
anova(m_0$lme, m_arma_p1$lme, m_arma_p2$lme, m_arma_p3$lme)

## Do we need AR(p > 1)? Depends which metric you use

## what is/are the estimated values of the AR(p) parameters?
intervals(m_ar1$lme, which = "var-cov")
intervals(m_arma_p1$lme, which = "var-cov")
intervals(m_arma_p2$lme, which = "var-cov")
intervals(m_arma_p3$lme, which = "var-cov")

## What about bam() - we can fit an AR(1), but need to provide rho
m_bam <- bam(Temperature ~ s(Year, k = 20),
    data = gtemp,
    method = "fREML",
    rho = 0.2 # <-- have to state a value of rho
)

acf(m_bam$std.rsd)

draw(m_bam)

## What do you do if the data aren't evenly spaced?
##
## Can use a continuous time AR(1) - CAR(1)
m_car1 <- gamm(Temperature ~ s(Year, k = 20),
               data = gtemp,
               method = "REML",
               correlation = corCAR1(form = ~ 1))
## estimated value of Phi, the CAR(1) term
intervals(m_car1$lme, which = "var-cov")
## should be close to the AR(1) value
