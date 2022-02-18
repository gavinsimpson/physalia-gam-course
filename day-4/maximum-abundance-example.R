# Analyse the simulated species abundance data

# load packages
pkgs <- c("mgcv", "ggplot2", "readr", "dplyr", "gratia")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

# load the data
spp_url <- "https://bit.ly/spp-gradient"
gradient <- read_csv(spp_url, col_types = "dd")
# or if you have the data cloned locally
# gradient <- read_csv(here("data", "simulated-gradient.csv"),
#                      col_types = "dd")
gradient

gradient %>%
  ggplot(aes(x = environment, y = abundance)) +
  geom_point()

## step 1 fit the GAM
max_nb <- gam(abundance ~ s(environment), data = gradient,
              method = "REML",
              family = nb())

draw(max_nb)

max_p <- update(max_nb, . ~ ., family = poisson())

draw(max_p)

AIC(max_nb, max_p)

nb_theta(max_nb)

## step 2
new_df <- with(gradient,
               tibble(environment = seq_min_max(environment, n = 200)))

fs <- fitted_samples(max_p, n = 10000, newdata = new_df,
                     seed = 42)

## a function to find the env value where abundance is maximum
max_loc <- function(row, abund, env) {
  r <- row[which.max(abund)]
  env[r]
}

max_post <- fs %>%
  group_by(draw) %>%
  summarise(env = max_loc(row, fitted, new_df$environment))

summ <- max_post %>%
  summarise(mean = mean(env), median = median(env),
            "q2.5" = quantile(env, prob = 0.025),
            "q97.5" = quantile(env, prob = 0.975))
summ

summ <- summ %>%
  mutate(abundance = 0)

## plot
gradient %>%
  ggplot(aes(x = environment, y = abundance)) +
  geom_point() +
  geom_pointrange(data = summ,
                  aes(y = abundance, x = median,
                      xmax = q97.5, xmin = q2.5), col = "red")
