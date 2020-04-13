nsim = 200

set.seed(202004)
r0postCrI = r0posterior
# hist(r0postCrI)
# summary(r0postCrI)
R0est = sample(x = r0postCrI,size = nsim)
# print(R0est)

## To simulate n_simSEIR outbreaks: duration of infection = 3 days, initial infected  n=~54 infected
epi_doNothingDurInf3 = vector('list',nsim)
epi_baseDurInf3 = vector('list',nsim)
epi_marchDurInf3 = vector('list',nsim)
epi_aprilDurInf3 = vector('list',nsim)
start = Sys.time()
durInfSim = 3
initialI = c(rep(0,4), rep(10^-6, 10), 0, 0) #assume 10 case before 10/3

for(sim in 1:nsim)
{
  epi_doNothingDurInf3[[sim]] = simulateOutbreakSEIR(R0t =R0est[sim] ,rho = rep(0.5,3660),dateStartSchoolClosure = as.Date('2020-03-10'),
                                                 dateStartIntenseIntervention = as.Date('2020-03-10'), dateEndIntenseIntervention = as.Date('2020-03-10'),
                                                 pWorkOpen = c(1,1,1,1),numWeekStagger = c(0,0,0),pInfected=initialI,durInf = durInfSim)
  epi_baseDurInf3[[sim]] = simulateOutbreakSEIR(R0t =R0est[sim] ,rho = rep(0.5,3660),dateStartSchoolClosure = as.Date('2020-03-10'), 
                                                dateStartIntenseIntervention  = as.Date('2020-04-01'),
                                                dateEndIntenseIntervention = as.Date('2020-04-15'),pWorkOpen = c(0.1, 0.25, 0.75, 0.9),
                                                numWeekStagger = c(2,2,4),pInfected=initialI,durInf = durInfSim)
  if(sim%%10==0) cat(paste0('Done with simulation ',sim), '\n')
}
end = Sys.time()
print(end-start)

covid_SDurInf3 = list() 
covid_SDurInf3[[1]] = summariseSimulations(VAR = 'S',CI = 95,SIMS = epi_doNothingDurInf3)
covid_SDurInf3[[2]] = summariseSimulations(VAR = 'S',CI = 95,SIMS = epi_baseDurInf3)

covid_IDurInf3 = list() 
covid_IDurInf3[[1]] = summariseSimulations(VAR = 'incidence',CI = 95,SIMS = epi_doNothingDurInf3)
covid_IDurInf3[[2]] = summariseSimulations(VAR = 'incidence',CI = 95,SIMS = epi_baseDurInf3)

peaktime_DurInf3 = list()
peaktime_DurInf3[[1]] = summarisePeakTimePeakSize(SIMS = epi_doNothingDurInf3)
peaktime_DurInf3[[2]] = summarisePeakTimePeakSize(SIMS = epi_baseDurInf3)

covid_DurInf3 = list() 
covid_DurInf3[[1]] = summariseSimulations_mid(CI = 95,SIMS = epi_doNothingDurInf3)
covid_DurInf3[[2]] = summariseSimulations_mid(CI = 95,SIMS = epi_baseDurInf3)

AGEcovid_IDurInf3 = list()
AGEcovid_IDurInf3[[1]] = summariseSimulationsAGE(VAR = 'incidence',CI = 95,SIMS = epi_doNothingDurInf3)
AGEcovid_IDurInf3[[2]] = summariseSimulationsAGE(VAR = 'incidence',CI = 95,SIMS = epi_baseDurInf3)

epiFirstSimDurInf3 = list(epi_doNothingDurInf3 = epi_doNothingDurInf3[[1]],
                          epi_baseDurInf3= epi_baseDurInf3[[1]])

# SEIR_dt = data.frame(time=covid_SDurInf3[[1]]$Sim1$time, S=covid_SDurInf3[[1]]$Sim1$S, E=covid_SDurInf3[[1]]$Sim1$E, I=covid_SDurInf3[[1]]$I,R=covid_SDurInf3[[1]]$Sim1$R)
# plot(covid_IDurInf3[[1]]$summary$mean)
# plot(covid_IDurInf3[[2]]$summary$mean)

covid_IIDurInf3 = list() 
covid_IIDurInf3[[1]] = summariseSimulations(VAR = 'I',CI = 95,SIMS = epi_doNothingDurInf3)
covid_IIDurInf3[[2]] = summariseSimulations(VAR = 'I',CI = 95,SIMS = epi_baseDurInf3)
# plot(covid_SDurInf3[[1]]$summary$mean)
# plot(covid_SDurInf3[[2]]$summary$mean)