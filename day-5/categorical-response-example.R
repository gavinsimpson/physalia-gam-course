# categorical
pkgs <- c("mgcv", "ggplot2", "dplyr", "gratia")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

# Ordered categorical model ocat()
n_categories <- 4
su_eg1_ocat <- data_sim("eg1", n = 200, dist = "ordered categorical",
                        n_cat = n_categories, seed = 2)
m_ocat <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3),
              family = ocat(R = n_categories), data = su_eg1_ocat, method = "REML")

fitted_values(m_ocat)

draw(m_ocat, scales = "fixed")

ds <- data_slice(m_ocat, x2 = evenly(x2))
fv_ocat <- fitted_values(m_ocat, data = ds)

fv_ocat |>
  ggplot(aes(x = x2, y = fitted, colour = category, group = category)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, x = x2,
                  fill = category, colour = NULL),
              alpha = 0.2) +
  geom_line()

head(predict(m_ocat, type = "response", exclude = "s(individual)"))