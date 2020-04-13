## To simulate n_simSEIcIscR outbreaks

nsim = 200

set.seed(123)
r0postCrI = r0posterior
# hist(r0postCrI)
# summary(r0postCrI)
R0est = sample(x = r0postCrI,size = nsim)
# print(R0est)

## To simulate n_sim SEIcIscR outbreaks: duration of infection = 3 days, initial infected  n=~200 infected
epi_doNothingDurInf3 = vector('list',nsim)
epi_baseDurInf3 = vector('list',nsim)
epi_marchDurInf3 = vector('list',nsim)
epi_aprilDurInf3 = vector('list',nsim)
start = Sys.time()
durInfSim = 3
initialI = c(rep(0,4), rep(10^-6, 10), 0, 0) #assume 10 case before 10/3
for(sim in 1:nsim)
{
  epi_doNothingDurInf3[[sim]] = simulateOutbreakSEIcIscR(R0t =R0est[sim] ,dateStartSchoolClosure = as.Date('2020-03-10'),
                                                 dateStartIntenseIntervention = as.Date('2020-03-10'), dateEndIntenseIntervention = as.Date('2020-03-10'),
                                                 pWorkOpen = c(1,1,1,1),numWeekStagger = c(0,0,0),pInfected=initialI,durInf = durInfSim)
  epi_baseDurInf3[[sim]] = simulateOutbreakSEIcIscR(R0t =R0est[sim], dateStartSchoolClosure = as.Date('2020-03-10'),
                                                    dateStartIntenseIntervention = as.Date('2020-04-01'),dateEndIntenseIntervention = as.Date('2020-04-15'),
                                                    pWorkOpen = c(.1,.25,.75,.9),
                                                    numWeekStagger = c(2,2,4),pInfected=initialI,durInf = durInfSim)
  if(sim%%10==0) cat(paste0('Done with simulation ',sim), '\n')
}
end = Sys.time()
print(end-start)

covid_SDurInf3sc = list() 
covid_SDurInf3sc[[1]] = summariseSimulations(VAR = 'S',CI = 95,SIMS = epi_doNothingDurInf3)
covid_SDurInf3sc[[2]] = summariseSimulations(VAR = 'S',CI = 95,SIMS = epi_baseDurInf3)


covid_IDurInf3sc = list() 
covid_IDurInf3sc[[1]] = summariseSimulations(VAR = 'incidence',CI = 95,SIMS = epi_doNothingDurInf3)
covid_IDurInf3sc[[2]] = summariseSimulations(VAR = 'incidence',CI = 95,SIMS = epi_baseDurInf3)

peaktime_DurInf3sc = list()
peaktime_DurInf3sc[[1]] = summarisePeakTimePeakSize(SIMS = epi_doNothingDurInf3)
peaktime_DurInf3sc[[2]] = summarisePeakTimePeakSize(SIMS = epi_baseDurInf3)

covid_DurInf3sc = list() 
covid_DurInf3sc[[1]] = summariseSimulations_mid(CI = 95,SIMS = epi_doNothingDurInf3)
covid_DurInf3sc[[2]] = summariseSimulations_mid(CI = 95,SIMS = epi_baseDurInf3)


AGEcovid_IDurInf3sc = list()
AGEcovid_IDurInf3sc[[1]] = summariseSimulationsAGE(VAR = 'incidence',CI = 95,SIMS = epi_doNothingDurInf3)
AGEcovid_IDurInf3sc[[2]] = summariseSimulationsAGE(VAR = 'incidence',CI = 95,SIMS = epi_baseDurInf3)

epiFirstSimDurInf3sc = list(epi_doNothingDurInf3 = epi_doNothingDurInf3[[1]],
                          epi_baseDurInf3= epi_baseDurInf3[[1]])

# library(ggplot2)
# library(patchwork)
# p1 = ggplot(mapping=aes(x=1:428, y=covid_IDurInf3sc[[1]]$summary$mean)) + geom_density(stat='identity') + geom_ribbon(aes(ymin=covid_IDurInf3sc[[1]]$summary$lci, ymax=covid_IDurInf3sc[[1]]$summary$uci), alpha=0.2)
# p2 = ggplot(mapping=aes(x=1:428, y=covid_IDurInf3sc[[2]]$summary$mean)) + geom_density(stat='identity') + geom_ribbon(aes(ymin=covid_IDurInf3sc[[2]]$summary$lci, ymax=covid_IDurInf3sc[[2]]$summary$uci), alpha=0.2)
# p1/p2

# covid_IIscDurInf3sc = list() 
# covid_IIscDurInf3sc[[1]] = summariseSimulations(VAR = 'Isc',CI = 50,SIMS = epi_doNothingDurInf3)
# covid_IIscDurInf3sc[[2]] = summariseSimulations(VAR = 'Isc',CI = 50,SIMS = epi_baseDurInf3)