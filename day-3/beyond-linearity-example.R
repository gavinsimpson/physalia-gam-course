library("mgcv")
library("ggplot2")

n <- 5000
set.seed(2)
df <- tibble(x = runif(n),
             y = x + x^2 * 0.2 + rnorm(n) * 0.1,
             f = x + x^2 * 0.2)

ggplot(df, aes(x = x, y = y)) + 
  geom_point() + 
  # geom_smooth(method = "gam") +
  geom_line(aes(x = x, y =f))

bfs <- basis(s(x, m = c(2, 0)), data = df)

draw(bfs) + facet_wrap(~ bf)

bfs <- basis(s(x, m = c(2)), data = df)
draw(bfs) + facet_wrap(~ bf)

model <- gam(y ~ x + s(x, m = c(2, 0)),
             data = df, method = "REML")

model0 <- gam(y ~ x,
             data = df, method = "REML")

summary(model)

draw(model)


model2 <- gam(y ~ s(x, m = c(2)),
             data = df, method = "REML")
draw(model2, overall_uncertainty = TRUE)

summary(model2)
