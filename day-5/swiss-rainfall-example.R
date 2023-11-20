library("gratia")
library("mgcv")

URL <- "https://bit.ly/gam-swiss-rainfall"

swiss <- read_csv(URL, col_types = "dddccdcdd")

# swiss <- read_csv(here("data", "swiss-rainfall.csv"),
#     col_types = "dddccdcdd")

m <- gam(list(exra ~ s(nao) + s(elevation) + climate.region + s(N, E),
                   ~ s(year) + s(elevation) + climate.region + s(N, E),
                   ~ climate.region),
         family = gevlss(),
         data = swiss)

draw(m, overall_uncertainty = FALSE) + plot_layout(widths = 1)

appraise(m)
