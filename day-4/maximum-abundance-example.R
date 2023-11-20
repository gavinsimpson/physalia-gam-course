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

## step 2

## a function to find the env value where abundance is maximum

## apply to the posterior

## summarise posterior

## plot