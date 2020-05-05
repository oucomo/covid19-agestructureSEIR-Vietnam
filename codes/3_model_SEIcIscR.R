nSim <- 5000 #number of simulations
nDaySim <- 366 #number of days simulated
POP <- (sum(wuhanpop$popage))*wuhanpop$propage # population per age group 
initialI <- c(rep(0,4), rep(1, 10), 0, 0)  # initial infectious people, currently all adults has one
rho <- c(rep(0.2, 4), rep(.5, 12)) # probability of being symptomatic per age group
sim_dat <-
  list(
    nDaySim = nDaySim,
    mean_DurInf = 10, #days in I for symptomatic cases
    mean_DurInfsc = 10, #days in I for asymptomatic/subclinical cases
    mean_DurLat = 5, #days in E
    mean_lnR0 = log(2), # mean of R0
    s_lnR0 = 1, # variance of log(R0)
    mean_lnR0postoutbreak = log(1.5), #mean of R0 postlockdown
    s_lnR0postoutbreak = 1, #variance of log(R0postlockdown)
    nAgeGroups = 16, #number of age groups
    POP = POP, 
    initialI = initialI,
    contact_matrix = contacts[-5], #delete "all", contact matrices per places:: HOME, WOKR, SCHOOL, OTHER
    mu_rho = rho,
    tStartIntenseIntervention = 9999,  #time of start lockdown, if > nDaySim => no lockdown
    # nIntenseStages = 4, # number of lockdown stages, differentiated by the prob of Work open, must be >1 due to limitation in stan
    # IntenseStageWeeks = c(4,2,2,4), # number of weeks in each stage
    # pWorkOpen = c(.1, .25, .5, .9), # prop of work open in each stag
    nIntenseStages=2, IntenseStageWeeks=c(0,0), pWorkOpen=c(1,1),
    tCloseSchool = 0,  # relative day of school closure
    tReopenSchool = 0 ## realative day of school reopen
  )

library(rstan)
options(mc.cores = 3)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
eigCalc <- tools::file_path_as_absolute('codes/stan/eigenCalc.hpp') # function to get the real part of pseudoeigenvalue
SEIcIscR_sim <- stan_model('codes/stan/model_SEIcIscR.stan',
                       allow_undefined = TRUE,
                       includes = paste0('\n#include "', eigCalc, '"\n'))
SEIcIscR_res <- sampling(SEIcIscR_sim, data = sim_dat, iter=nSim, chains=3, warmup=1000, seed=296, control=list(adapt_delta=0.7))
SEIcIscR_ext <- extract_SEIR(SEIcIscR_res)

saveRDS(SEIcIscR_ext, file='outputs/SEIcIscR.RDS')
