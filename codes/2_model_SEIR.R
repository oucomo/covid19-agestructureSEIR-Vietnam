nSim <- 5000 # number of simulations
nDaySim <- 366 # number of days simulated
POP <- (sum(wuhanpop$popage))*wuhanpop$propage # get the population for each age group
initialI <- c(rep(0,4), rep(1, 10), 0, 0) # number of infected a t0, current all adults
rho <- rep(0.8, nDaySim) # pct of reported cases reported = incidence * rho
sim_dat <-
    list(
        nDaySim = nDaySim, # number of day simulated
        mean_DurInf = 10, # days in group I
        mean_DurLat = 7, # day in group E
        mean_R0 = 2, # mean of R0
        s_R0 = .4, # variance of log(R0)
        mean_R0postoutbreak = 1.5, # mean of R0 post lockdown
        s_R0postoutbreak = .2, #variance of log(R0 post lockdowm)
        nAgeGroups = 16, # number of age groups
        POP = POP, 
        initialI = initialI,
        contact_matrix = contacts[-5], #delete "all", statified by places
        rho = rho, 
        tStartIntenseIntervention = 99999, #time of start lockdown, if > nDaySim => no lockdown
        nIntenseStages = 2, # number of lockdown stages, differentiated by the prob of Work open, must be >1 due to limitation in stan
        IntenseStageWeeks = c(0,0), # number of weeks in each stage
        pWorkOpen = c(1,1), # prop of work open in each stage
        tCloseSchool = 0, # relative day of school closure
        tReopenSchool = 0 # realative day of school reopen
    )

library(rstan)
options(mc.cores = 3)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

eigCalc <- tools::file_path_as_absolute('codes/stan/eigenCalc.hpp') # function to get the real part of pseudoeigenvalue
SEIR_sim <- stan_model('codes/stan/model_SEIR.stan',
    allow_undefined = TRUE,
    includes = paste0('\n#include "', eigCalc, '"\n'))
SEIR_res <- sampling(SEIR_sim, data = sim_dat, iter=nSim, chains=3, warmup=1000, seed=296, control=list(adapt_delta=0.8))
SEIR_ext <- extract_SEIR(SEIR_res)

saveRDS(SEIR_ext, file='outputs/SEIR.RDS')
