nSim <- 5000
POP <- (sum(wuhanpop$popage))*wuhanpop$propage 
initialI <- c(rep(0,4), rep(10^-6, 10), 0, 0) 
rho <- rep(0.8, 430)
sim_dat <-
    list(
        nDaySim = 430,
        mean_DurInf = 10,
        # s_DurInf = 5,
        mean_DurLat = 7,
        # s_DurLat = 7,
        mean_R0 = 2,
        s_R0 = .5,
        mean_R0postoutbreak = 1.5,
        s_R0postoutbreak = .1,
        nAgeGroups = 16,
        POP = POP,
        initialI = initialI,
        contact_matrix = contacts[-5], #delete "all"
        rho = rho,
        tStartIntenseIntervention = 22,
        nIntenseStages = 4,
        IntenseStageWeeks = c(2,2,4,4),
        pWorkOpen = c(.1, .25, .5, .9),
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
SEIR_res <- sampling(SEIR_sim, data = sim_dat, iter=nSim, chains=1, warmup=1000, seed=200408, control=list(adapt_delta=0.8))
SEIR_ext <- extract(SEIR_res)

S=sapply(1:430, function(i) sapply(1:16, function(j) mean(SEIR_ext$SEIR[,1,i,j])))
E=sapply(1:430, function(i) sapply(1:16, function(j) mean(SEIR_ext$SEIR[,2,i,j])))
I=sapply(1:430, function(i) sapply(1:16, function(j) mean(SEIR_ext$SEIR[,3,i,j])))
R=sapply(1:430, function(i) sapply(1:16, function(j) mean(SEIR_ext$SEIR[,4,i,j])))
