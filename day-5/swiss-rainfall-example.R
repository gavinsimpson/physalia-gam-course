
library("gratia")
library("mgcv")

swiss <- read_csv(here("data", "swiss-rainfall.csv"),
    col_types = "dddccdcdd")

m <- gam(list(exra ~ s(nao) + s(elevation) + climate.region + s(N, E),
                   ~ s(year) + s(elevation) + climate.region + s(N, E),
                   ~ climate.region),
         family = gevlss(),
         data = swiss)

draw(m)

appraise(m)

