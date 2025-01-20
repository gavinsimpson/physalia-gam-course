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
new_df <- data_slice(gradient, environment = evenly(environment, n = 200)) |>
  select(environment)

fs <- fitted_samples(max_p, n = 10000, data = new_df, seed = 42)

## a function to find the env value where abundance is maximum
max_loc <- function(row, abund, env) {
  r <- row[which.max(abund)]
  env[r]
}

fv <- fitted_values(max_p, data = new_df)

gradient %>%
  ggplot(aes(x = environment, y = abundance)) +
  geom_point() +
  geom_line(data = fv, mapping = aes(y = .fitted, x = environment),
            colour = "red")

max_loc(fv$.row, fv$.fitted, fv$environment)

fv |> slice_max(order_by = .fitted)

# apply fun to posterior simulations

max_loc(df1$.row, df1$.fitted, new_df$environment)

max_post <- fs %>%
  group_by(.draw) %>%
  summarise(env = max_loc(.row, .fitted, new_df$environment))

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

## using the Metropolis Hasting sampler

fs_mh <- fitted_samples(max_p, n = 10000, data = new_df, seed = 42,
                        method = "mh", burnin = 1000, thin = 2, t_df = 5)

max_post_mh <- fs_mh %>%
  group_by(.draw) %>%
  summarise(env = max_loc(.row, .fitted, new_df$environment))

summ_mh <- max_post_mh %>%
  summarise(mean = mean(env), median = median(env),
            "q2.5" = quantile(env, prob = 0.025),
            "q97.5" = quantile(env, prob = 0.975))
summ_mh

summ_mh <- summ_mh %>%
  mutate(abundance = 0)

## plot
gradient %>%
  ggplot(aes(x = environment, y = abundance)) +
  geom_point() +
  geom_pointrange(data = summ,
    aes(y = abundance, x = median,
      xmax = q97.5, xmin = q2.5), col = "red")
