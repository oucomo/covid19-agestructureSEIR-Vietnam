nSim <- 3000
nDaySim <- 366
POP <- (sum(wuhanpop$popage))*wuhanpop$propage 
initialI <- c(rep(0,4), rep(1, 10), 0, 0) 
rho <- rep(0.8, nDaySim)
sim_dat <-
    list(
        nDaySim = nDaySim,
        mean_DurInf = 10,
        mean_DurLat = 7,
        mean_R0 = 2,
        s_R0 = .4,
        mean_R0postoutbreak = 1.5,
        s_R0postoutbreak = .2,
        nAgeGroups = 16,
        POP = POP,
        initialI = initialI,
        contact_matrix = contacts[-5], #delete "all"
        rho = rho,
        tStartIntenseIntervention = 22,
        nIntenseStages = 4,
        IntenseStageWeeks = c(4,2,4,4),
        pWorkOpen = c(.1, .25, .5, .75),
        tCloseSchool = 1,
        tReopenSchool = 6*7
    )

library(rstan)
options(mc.cores = parallel::detectCores()-1)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
# SEIR_sim <- stan('codes/stan/model_SEIR.stan', data = sim_dat, iter=nSim, chains=1, seed=200408)
eigCalc <- tools::file_path_as_absolute('codes/stan/eigenCalc.hpp')
SEIR_sim <- stan_model('codes/stan/model_SEIR.stan',
    allow_undefined = TRUE,
    includes = paste0('\n#include "', eigCalc, '"\n'))
SEIR_res <- sampling(SEIR_sim, data = sim_dat, iter=nSim, chains=4, warmup=1000, seed=296, control=list(adapt_delta=0.8))
SEIR_ext <- extract_SEIR(SEIR_res)

