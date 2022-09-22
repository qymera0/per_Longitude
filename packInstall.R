packages <- c("bayesplot", "brms", "broom", "devtools", "flextable", "GGally", "ggmcmc", "ggrepel", "gtools", "loo", "patchwork", "psych", "Rcpp", "remotes", "rstan", "StanHeaders", "survival", "tidybayes", "tidyverse")

install.packages(packages, dependencies = T)

install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

remotes::install_github("stan-dev/posterior")
devtools::install_github("rmcelreath/rethinking")