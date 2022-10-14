# braya so example
braya <- read_table(here("data", "braya-so.txt"), skip = 84,
    col_names = FALSE, col_types = "ddddddd") %>%
    setNames(c("Depth", "DepthUpper", "DepthLower", "Year", "YearYoung",
        "YearOld", "UK37")) %>%
    mutate(sampleInterval = YearYoung - YearOld)

braya

braya_lab <- expression(italic(U)[37]^{italic(k)})

ggplot(braya, aes(x = Year, y = UK37)) +
    geom_line(colour = "grey") +
    geom_point() +
    labs(y = braya_lab, x = "Year CE")

## fit the car(1) model --- needs optim as this is not a stable fit!
## also needs k setting lower than default
braya_car1 <- gamm(UK37 ~ s(Year, k = 30), data = braya,
                   correlation = corCAR1(form = ~ Year),
                   method = "REML",
		           control = list(niterEM = 0, optimMethod = "BFGS",
                                  opt = "optim"))

## fit model using GCV
braya_gcv <- gam(UK37 ~ s(Year, k = 30), data = braya)

draw(braya_gcv, n = 200)

## estimate of phi and confidence interval
intervals(braya_car1$lme)$corStruct

# do this properly
braya_m <- gam(UK37 ~ s(Year, k = 40), data = braya,
                  method = "REML",
                  weights = sampleInterval / mean(sampleInterval))
summary(braya_m)

draw(braya_m, n = 200)

braya_m2 <- gam(list(UK37 ~ s(Year, k = 40),
    ~ s(sampleInterval)),
data = braya,
method = "REML",
family = gaulss())

draw(braya_m2, n = 200, overall_uncertainty = FALSE)

