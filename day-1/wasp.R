# Darlingtonia example
library("readr")
library("dplyr")
library("ggplot2")
library("here")

# load data
wasp <- read_csv(here("data", "darlingtonia.csv"), comment = "#",
                 col_types = "dl")

m <- glm(visited ~ leafHeight, data = wasp, family = binomial(link = "logit"))
summary(m)

pdat <- with(wasp,
             tibble(leafHeight = seq(min(leafHeight),
                                     max(leafHeight),
                                     length = 100)))
pdat

# predict
pred <- predict(m, pdat, type = "link", se.fit = TRUE)
str(pred)

p_link <- predict(m, pdat, type = "response")
head(p_link)

# inverse link
ilink <- family(m)$linkinv # g-1()

# add prediction data and predicted values
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

# likelihood ratio test
# compare m with a model without leafHeight
m0 <- update(m, . ~ . - leafHeight)

# GLRT
anova(m0, m, test = "LRT")

# what about polynomials of leafHeight?
m2 <- update(m, . ~ . + I(leafHeight^2))
summary(m2)

anova(m, m2, test = "LRT")
