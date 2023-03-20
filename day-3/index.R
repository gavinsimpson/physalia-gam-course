## ----setup, include=FALSE, cache=FALSE----------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(cache = TRUE, dev = "svg", echo = TRUE, message = FALSE, 
  warning = FALSE, fig.height = 6, fig.width = 1.777777 * 6)

library("gridGraphics")
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
theme_set(theme_bw(base_size = 16, base_family = 'Fira Sans'))

## constants
anim_width <- 1000
anim_height <- anim_width / 1.77777777
anim_dev <- 'png'
anim_res <- 200


## ----tprs-2d-basis-setup, echo = FALSE----------------------------------------
test_fun2 <- function(n = 50, sx = 0.3, sz = 0.4, scale = 0.1,
                     seed = NULL) {
    require("tibble")
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        runif(1)
    }
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    }
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    xz <- expand.grid(x = seq(0,1, length = n), z = seq(0,1, length=n))
    #take <- seq_len(n)
    x <- xz[["x"]]
    z <- xz[["z"]]
    f <- 1.2 * exp(-(x - 0.2)^2 / sx^2 - (z - 0.3)^2 / sz^2) +
      08. * exp(-(x - 0.7)^2 / sx^2 - (z - 0.8)^2 / sz^2)
    tibble(y = f + rnorm(n) * scale, x = x, z = z)
}

library("gratia")
library("mgcv")
library("ggplot2")
library("dplyr")


## ----tprs-2d-basis, echo = FALSE, out.width = "90%"---------------------------
test_df2 <- test_fun2(seed = 42)
bfun <- basis(s(x, z), data = test_df2)

bfun %>%
  ggplot(aes(x = x, y = z, fill = value, group = bf)) +
  geom_raster() +
  facet_wrap(~ bf) +
  scale_fill_distiller(palette = "RdBu", type = "div") +
  scale_x_continuous(guide = guide_axis(n.dodge = 2,
                                        check.overlap = TRUE)) +
  theme_bw(base_size = 12, base_family = 'Fira Sans')


## ----include-simons-tprs-2d-image, out.width = "45%", echo = FALSE------------
knitr::include_graphics("resources/wood-2ed-fig-5-12-2-d-tprs-basis-funs.png")


## ----beyond-linearity---------------------------------------------------------
n <- 100
set.seed(2)
df <- tibble(x = runif(n),
             y = x + x^2 * 0.2 + rnorm(n) * 0.1)

model <- gam(y ~ x + s(x, m = c(2, 0)),
             data = df, method = "REML")


## ----beyond-linearity-summary-------------------------------------------------
summary(model)


## ----beyond-linearity-plot, out.width = "95%", fig.align = "center"-----------
draw(model, parametric = TRUE)


## ----draw-mcycle, out.width = "70%", fig.align = "center"---------------------
data(mcycle, package = "MASS")
m <- gam(accel ~ s(times), data = mcycle, method = "REML")
sm_plt <- draw(m, residuals = TRUE)
sm_plt


## ----derivatives-times-smooth-------------------------------------------------
fd <- derivatives(m, type = "central", unconditional = TRUE)
fd


## ----draw-derivatives-times-smooth, fig.align = "center", out.width = "90%"----
fd_plt <- draw(fd) + labs(title = "First derivative s(times)")
sm_plt + fd_plt + plot_layout(ncol = 2)


## ----draw-derivatives-times-smooth-2, fig.align = "center", out.width = "80%"----
fd2 <- derivatives(m, order = 2, eps = 0.1, type = "central", unconditional = TRUE)
fd2_plt <- draw(fd2)
fd_plt + fd2_plt + plot_layout(ncol = 2)


## ----draw-derivatives-times-smooth-2-right, fig.align = "center", out.width = "80%"----
m2 <- gam(accel ~ s(times, m = 3), data = mcycle, method = "REML")
fd2 <- derivatives(m2, order = 2, eps = 0.1, type = "central", unconditional = TRUE)
fd2_plt <- draw(fd2) + labs(title = "Second derivative s(times)")
fd_plt + fd2_plt + plot_layout(ncol = 2)


## ----tprs-vs-tensor-product-setup, echo = FALSE-------------------------------
# following shows how tensor pruduct deals nicely with 
# badly scaled covariates (range of x 5% of range of z )
test1 <- function(x, z, sx = 0.3, sz = 0.4) {
  x <- x * 20
  (pi ** sx * sz) * (1.2 * exp(-(x - 0.2)^2 / sx^2 - (z - 0.3)^2 / sz^2) +
    0.8 * exp(-(x - 0.7)^2 / sx^2 - ( z - 0.8)^2 / sz^2))
}
n <- 500
x <- runif(n) / 20
z <- runif(n)
xs <- seq(0, 1, length = 30) / 20
zs <- seq(0, 1, length = 30)
pr <- tibble(x = rep(xs, 30), z = rep(zs, rep(30, 30)))
truth <- matrix(test1(pr$x, pr$z), 30, 30)
f <- test1(x, z)
y <- f + rnorm(n) * 0.2
df <- tibble(y = y, x = x, z = z)
truth_df <- pr %>% mutate(f = test1(x, z))
m_tprs <- gam(y ~ s(x, z), data = df, method = "REML")
m_te <- gam(y ~ te(x, z), data = df, method = "REML")

truth_plt <- truth_df %>%
  ggplot(aes(x = x, y = z, fill = f)) +
    geom_raster() +
    scale_fill_distiller(palette = "RdBu", type = "div") +
    geom_contour(aes(z = f), colour = "black", bins = 8) +
    labs(title = "f(x,z)")


## ----draw-tprs-vs-tensor-product-truth, echo = FALSE--------------------------
old_par <- par(mar = c(0, 2, 0, 0), bg = NA)
truth_plt +
  wrap_elements(panel = ~ persp(xs, zs, truth), clip = FALSE) +
  plot_layout(ncol = 2)
par(old_par)


## ----tprs-vs-tensor-product---------------------------------------------------
df
m_tprs <- gam(y ~ s(x, z), data = df, method = "REML")
m_te   <- gam(y ~ te(x, z), data = df, method = "REML")


## ----draw-tprs-vs-tensor-product, message = FALSE, out.width = "95%"----------
truth_plt + (draw(m_tprs) + coord_cartesian()) + draw(m_te) + plot_layout(ncol = 3)


## ----plot-tprs-vs-tensor-product-fake, eval = FALSE---------------------------
## layout(matrix(1:3, ncol = 3))
## persp(xs, zs, truth)
## vis.gam(m_tprs)
## vis.gam(m_te)
## layout(1)


## ----plot-tprs-vs-tensor-product, echo = FALSE, fig.width = 6, fig.width = 18----
old_par <- par(mar = c(0, 2, 0, 0), bg = NA)
persp1 <- wrap_elements(panel = ~ persp(xs, zs, truth), clip = FALSE)
persp2 <- wrap_elements(panel = ~ vis.gam(m_tprs), clip = FALSE)
persp3 <- wrap_elements(panel = ~ vis.gam(m_te), clip = FALSE)
plt <- persp1 + labs(title = "Truth") + 
  persp2 + labs(title = "TPRS") +
  persp3 + labs(title = "Tensor Product") +
  plot_layout(ncol = 3)
plt
par(old_par)


## ----echo = FALSE, out.width = "50%"------------------------------------------
knitr::include_graphics("resources/wood-gams-2ed-fig-5-17-tensor-product.svg")


## ----ranefs-------------------------------------------------------------------
m_nlme <- lme(travel ~ 1, data = Rail, ~ 1 | Rail, method = "REML") 

m_gam  <- gam(travel ~ s(Rail, bs = "re"), data = Rail, method = "REML")


## ----misspecify, echo = FALSE, out.width = "95%"------------------------------
set.seed(15)
model_list = c("right model",
               "wrong distribution",
               "heteroskedasticity",
               "dependent data",
               "wrong functional form")
n <- 60
sigma=1
x <- seq(-1,1, length=n)
model_data <- as.data.frame(expand.grid( x=x,model=model_list))
model_data$y <- 5*model_data$x^2 + 2*model_data$x
for(i in model_list){
  if(i == "right model"){
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+ 
      rnorm(n,0, sigma)
  } else if(i == "wrong distribution"){
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+ 
      rt(n,df = 3)*sigma
  } else if(i == "heteroskedasticity"){
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+  
      rnorm(n,0, sigma*10^(model_data[model_data$model==i, "x"]))
  } else if(i == "dependent data"){
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+ 
      arima.sim(model = list(ar=c(.7)), n = n,sd=sigma) 
  } else if(i=="wrong functional form") {
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+ 
      rnorm(n,0, sigma) + ifelse(model_data[model_data$model==i, "x"]>0, 5,-5)
  }
}
ggplot(aes(x,y), data= model_data)+
  geom_point()+
  geom_line(color=ifelse(model_data$model=="dependent data", "black",NA))+
  facet_wrap(~model)+
  geom_smooth(method=gam, formula = y~s(x,k=12),method.args = list(method="REML"))+
  theme(strip.text = element_text(size=16))


## ----sims, include=TRUE,echo=TRUE---------------------------------------------
set.seed(2)
n <- 400
x1 <- rnorm(n)
x2 <- rnorm(n)
y_val <- 1 + 2*cos(pi*x1) + 2/(1+exp(-5*(x2)))
y_norm <- y_val + rnorm(n, 0, 0.5)
y_negbinom <- rnbinom(n, mu = exp(y_val),size=10)
y_binom <- rbinom(n,1,prob = exp(y_val)/(1+exp(y_val)))


## ----sims_plot,fig.width = 11, fig.height = 5.5, echo = FALSE-----------------
p1 <- ggplot(data.frame(x = x1, y = y_norm),
             aes(x = x, y = y)) +
    geom_point() + labs(x = "x1", title = "Gaussian")

p2 <- ggplot(data.frame(x = x2, y = y_norm),
             aes(x = x, y = y)) +
    geom_point() + labs(x = "x2", title = "Gaussian")

p3 <- ggplot(data.frame(x = x1, y = y_negbinom),
             aes(x = x, y = y)) +
    geom_point() + labs(x = "x1", title = "Negative binomial")

p4 <- ggplot(data.frame(x = x2, y = y_negbinom),
             aes(x = x, y = y)) +
    geom_point() + labs(x = "x2", title = "Negative binomial")

p5 <- ggplot(data.frame(x = x1, y = y_binom),
             aes(x = x, y = y)) +
    geom_point() + labs(x = "x1", title = "Binomial")

p6 <- ggplot(data.frame(x = x2, y = y_binom),
             aes(x = x, y = y)) +
    geom_point() + labs(x = "x2", title = "Binomial")

#plot_grid(p1, p3, p5, p2, p4, p6, ncol = 3, align = 'hv', axis = 'lrtb')
wrap_plots(p1, p3, p5, p2, p4, p6, ncol = 3)


## ----gam_check_norm1, fig.keep="none", include=TRUE,echo=TRUE, fig.width=11, fig.height = 5.5, fig.align="center"----
norm_model_1 <- gam(y_norm ~ s(x1, k = 4) + s(x2, k = 4), method = 'REML')
gam.check(norm_model_1)


## ----gam_check_norm2, fig.keep="none", include=TRUE, echo=TRUE, fig.width=15, fig.height = 5.5,fig.align="center"----
norm_model_2 <- gam(y_norm ~ s(x1, k = 12) + s(x2, k = 4), method = 'REML')
gam.check(norm_model_2)


## ----gam_check_norm3, fig.keep="none", include=TRUE, echo=TRUE----------------
norm_model_3 <- gam(y_norm ~ s(x1, k = 12) + s(x2, k = 12),method = 'REML')
gam.check(norm_model_3)


## ----gam_check_norm4, echo = FALSE--------------------------------------------
p1 <- draw(norm_model_1)
p2 <- draw(norm_model_2)
p3 <- draw(norm_model_3)

## plot_grid(p1, p2, p3, nrow = 3, align = 'hv', axis = 'lrtb')
wrap_plots(p1, p2, p3, nrow = 3)


## ----alt-basis-dim-check-1----------------------------------------------------
norm_model_1 <- gam(y_norm ~ s(x1, k = 4) + s(x2, k = 4), method = "REML")
k.check(norm_model_1)


## ----alt-basis-dim-check-2----------------------------------------------------
res <- resid(norm_model_1, type = "deviance")

res_model <- gam(res ~ s(x1, k = 12) + s(x2, k = 12),
  method = "REML",
  family = quasi(link = "identity", variance = "constant"))
edf(res_model)


## ----alt-basis-dim-check-3, fig.align = "center", out.width = "95%"-----------
draw(res_model)


## ----gam_check_plots1, include=TRUE, echo=TRUE, results="hide", out.width = "90%", fig.align = "center"----
norm_model <- gam(y_norm ~ s(x1, k=12) + s(x2, k=12), method="REML")
gam.check(norm_model, rep = 500)


## ----gam_check_plots2, include=T, echo=TRUE, results="hide", out.width = "90%", fig.align = "center"----
pois_model <- gam(y_negbinom ~ s(x1, k=12) + s(x2, k=12), family=poisson, method="REML")
gam.check(pois_model, rep = 500)


## ----gam_check_plots3, include=T,echo=TRUE, results="hide", out.width = "90%", fig.align = "center"----
negbin_model <- gam(y_negbinom ~ s(x1, k=12) + s(x2, k=12), family = nb, method="REML")
gam.check(negbin_model, rep = 500)


## ----appraise-gam-check-example, fig.height = 5.5-----------------------------
appraise(negbin_model, method = 'simulate')


## ----dharma-residuals-image, fig.align = "center", out.width = "50%", echo = FALSE----
knitr::include_graphics("resources/dharma-randomised-residuals.png")


## ----dharma-1-----------------------------------------------------------------
library("mgcViz")
library("DHARMa")

testDispersion(pois_model, plot = FALSE)


## ----dharma-2-----------------------------------------------------------------
testDispersion(negbin_model, plot = FALSE)


## ----dharma-sim-resids--------------------------------------------------------
resids <- simulateResiduals(fittedModel = pois_model, plot = FALSE)


## ----dharma-plots-possion, fig.align = "center", out.width = "90%"------------
plot(resids)


## ----dharma-plots-negbin, fig.align = "center", out.width = "90%"-------------
resids <- simulateResiduals(fittedModel = negbin_model, plot = FALSE)
plot(resids)


## ----setup-shrinkage-example--------------------------------------------------
## an example of automatic model selection via null space penalization
n <- 200
dat <- data_sim("eg1", n = n, scale = .15, dist = "poisson", seed = 3) ## simulate data
set.seed(21)
dat <- dat %>% mutate(x4 = runif(n, 0, 1), x5 = runif(n, 0, 1),
                      f4 = rep(0, n), f5 = rep(0, n))   ## spurious

## ----shrinkage-example-model-fit, echo = TRUE---------------------------------
b <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3) +
             s(x4) + s(x5),
         data = dat, family = poisson, method = 'REML',
         select = TRUE)


## ----shrinkage-example-truth, echo = FALSE, dependson=-1----------------------
p1 <- ggplot(dat, aes(x = x0, y = f0)) + geom_line()
p2 <- ggplot(dat, aes(x = x1, y = f1)) + geom_line()
p3 <- ggplot(dat, aes(x = x2, y = f2)) + geom_line()
p4 <- ggplot(dat, aes(x = x3, y = f3)) + geom_line()
p5 <- ggplot(dat, aes(x = x4, y = f4)) + geom_line()
p6 <- ggplot(dat, aes(x = x5, y = f5)) + geom_line()
#plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, align = 'vh', labels = paste0('x', 1:6))
p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "1", tag_prefix = "x")


## ----shrinkage-example-summary, dependson=-1----------------------------------
summary(b)


## ----shrinkage-example-plot, fig.align = "center", out.width = "95%", dependson=-1----
draw(b, scales = 'fixed')


## ----cross-validated, echo = FALSE--------------------------------------------
knitr::include_graphics("resources/cross-validated.png")


## ----load-galveston-----------------------------------------------------------
galveston <- read_csv("https://bit.ly/gam-galveston") %>%
    mutate(datetime = as.POSIXct(paste(DATE, TIME),
                                 format = '%m/%d/%y %H:%M', tz = "CDT"),
           STATION_ID = factor(STATION_ID),
           DoY = as.numeric(format(datetime, format = '%j')),
           ToD = as.numeric(format(datetime, format = '%H')) +
               (as.numeric(format(datetime, format = '%M')) / 60))
galveston


## ----galveston-full-model-----------------------------------------------------
knots <- list(DoY = c(0.5, 366.5))
m <- bam(MEASUREMENT ~
             s(ToD, k = 10) +
             s(DoY, k = 12, bs = "cc") +
             s(YEAR, k = 30) +
             s(LONGITUDE, LATITUDE, k = 100, bs = "ds", m = c(1, 0.5)) +
             ti(DoY, YEAR, bs = c("cc", "tp"), k = c(12, 15)) +
             ti(LONGITUDE, LATITUDE, ToD, d = c(2,1), bs = c("ds", "tp"),
                m = list(c(1, 0.5), NA), k = c(20, 10)) +
             ti(LONGITUDE, LATITUDE, DoY, d = c(2,1), bs = c("ds", "cc"),
                m = list(c(1, 0.5), NA), k = c(25, 12)) +
             ti(LONGITUDE, LATITUDE, YEAR, d = c(2,1), bs = c("ds", "tp"),
                m = list(c(1, 0.5), NA), k = c(25, 15)),
         data = galveston, method = "fREML", knots = knots,
         nthreads = c(6, 1), discrete = FALSE)


## ----galveston-simple-model---------------------------------------------------
m.sub <- bam(MEASUREMENT ~
             s(ToD, k = 10) +
             s(DoY, k = 12, bs = "cc") +
             s(YEAR, k = 30) +
             s(LONGITUDE, LATITUDE, k = 100, bs = "ds", m = c(1, 0.5)) +
             ti(DoY, YEAR, bs = c("cc", "tp"), k = c(12, 15)),
         data = galveston, method = "fREML", knots = knots,
         nthreads = c(4, 1), discrete = TRUE)


## ----galveston-compare-models-aic---------------------------------------------
AIC(m, m.sub)


## ----galveston-compare-models-anova-------------------------------------------
anova(m, m.sub, test = "F")


## ----galveston-full-model-summary---------------------------------------------
summary(m)


## ----galveston-full-model-plot, fig.height = 5.5------------------------------
plot(m, pages = 1, scheme = 2, shade = TRUE)


## ----galveston-full-model-draw, fig.height = 14, fig.width = 1.777777*14, fig.align = "center", out.width = "90%"----
draw(m, scales = "free", rug = FALSE, n = 50) +  plot_layout(widths = 1) &
  theme(strip.text.x = element_text(size = 8))


## ----galveston-full-predict---------------------------------------------------
pdata <- data_slice(m, ToD = 12, DoY = 180,
                    YEAR = evenly(YEAR, by = 1),
                    LONGITUDE = evenly(LONGITUDE, n = 50),
                    LATITUDE  = evenly(LATITUDE, n = 50))
fv <- fitted_values(m, data = pdata)
# set fitted values to NA for grid points that are too far from the data
ind <- too_far(pdata$LONGITUDE, pdata$LATITUDE,
               galveston$LONGITUDE, galveston$LATITUDE, dist = 0.1)
fv <- fv %>%
  mutate(fitted = if_else(ind, NA_real_, fitted))


## ----galveston-full-predict-plot, fig.show = 'hide', fig.height = 10, fig.width = 1.777777*10----
plt <- ggplot(fv, aes(x = LONGITUDE, y = LATITUDE)) +
    geom_raster(aes(fill = fitted)) + facet_wrap(~ YEAR, ncol = 12) +
    scale_fill_viridis(name = expression(degree*C), option = "plasma",
      na.value = "transparent") +
    coord_quickmap() +
    scale_x_continuous(guide = guide_axis(n.dodge = 2,
                                          check.overlap = TRUE)) +
    theme(legend.position = "top")
plt


## ----galveston-full-predict-plot, echo = FALSE, fig.height = 10, fig.width = 1.777777*10----
plt <- ggplot(fv, aes(x = LONGITUDE, y = LATITUDE)) +
    geom_raster(aes(fill = fitted)) + facet_wrap(~ YEAR, ncol = 12) +
    scale_fill_viridis(name = expression(degree*C), option = "plasma",
      na.value = "transparent") +
    coord_quickmap() +
    scale_x_continuous(guide = guide_axis(n.dodge = 2,
                                          check.overlap = TRUE)) +
    theme(legend.position = "top")
plt


## ----galveston-animation, echo = FALSE, results = 'hide'----------------------
p <- ggplot(fv, aes(x = LONGITUDE, y = LATITUDE, frame = YEAR)) +
    geom_raster(aes(fill = fitted)) +
    scale_fill_viridis(name = expression(degree*C), option = "plasma",
                       na.value = "transparent") +
    coord_quickmap() +
    theme(legend.position = "right") +
    labs(x = "Longitude", y = "Latitude")

anim <- p + transition_time(YEAR) +
    ggtitle("Year {round(frame_time, 0)}")

anim <- animate(anim,
                nframes = 200, height = anim_height, width = anim_width,
                res = 100, dev = anim_dev)

anim_save('./resources/galveston-animation.gif', anim)


## ----galveston-trends-by-month, fig.show = "hide"-----------------------------
ds <- data_slice(m, ToD = 12, DoY = c(1, 90, 180, 270),
  YEAR = evenly(YEAR, n = 250),
  LONGITUDE = -94.8751, LATITUDE  = 29.50866)
fv <- fitted_values(m, data = ds, scale = "response")

plt2 <- ggplot(fv, aes(x = YEAR, y = fitted, group = factor(DoY))) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.2) +
    geom_line() + facet_wrap(~ DoY, scales = "free_y") +
    labs(x = NULL, y = expression(Temperature ~ (degree * C)))
plt2


## ----galveston-trends-by-month, echo = FALSE----------------------------------
ds <- data_slice(m, ToD = 12, DoY = c(1, 90, 180, 270),
  YEAR = evenly(YEAR, n = 250),
  LONGITUDE = -94.8751, LATITUDE  = 29.50866)
fv <- fitted_values(m, data = ds, scale = "response")

plt2 <- ggplot(fv, aes(x = YEAR, y = fitted, group = factor(DoY))) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.2) +
    geom_line() + facet_wrap(~ DoY, scales = "free_y") +
    labs(x = NULL, y = expression(Temperature ~ (degree * C)))
plt2

