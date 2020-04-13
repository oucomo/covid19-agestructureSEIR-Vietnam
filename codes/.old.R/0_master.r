# load relevant the data files
source('codes/1_loadData.r')

# source the age-structured SEIcIscR model functions 
source('codes/function_modelSEIcIscR.r')

# source the age-structured SEIcIscR model functions 
source('codes/function_postprocessing.r')

# simulate N oubtreaks
source('codes/2_simOutbreak_ncov_SEIR_vn.r')
source('codes/2_simOutbreak_ncov_SEIcIscR_vn.r')


library(ggplot2)
library(patchwork)
p1 = ggplot(mapping=aes(x=1:428, y=covid_IDurInf3[[1]]$summary$mean)) + geom_density(stat='identity') + geom_ribbon(aes(ymin=covid_IDurInf3[[1]]$summary$lci, ymax=covid_IDurInf3[[1]]$summary$uci), alpha=0.2) + xlab('times') + ylab('I')
# p2 = ggplot(mapping=aes(x=1:428, y=covid_IDurInf3sc[[1]]$summary$mean)) + geom_density(stat='identity') + geom_ribbon(aes(ymin=covid_IDurInf3sc[[1]]$summary$lci, ymax=covid_IDurInf3sc[[1]]$summary$uci), alpha=0.2) + xlab('times') + ylab('I')
# print(p1|p2)

p2+geom_density(mapping=aes(y=covid_IDurInf3sc[[1]]$summary$mean), color='#963939', stat='identity') +
 geom_ribbon(aes(ymin=covid_IDurInf3sc[[1]]$summary$lci, ymax=covid_IDurInf3sc[[1]]$summary$uci),fill='#963939', alpha=0.2)


p1+geom_density(mapping=aes(x=1:428, y=covid_IDurInf3[[2]]$summary$mean), stat='identity') + geom_ribbon(aes(ymin=covid_IDurInf3[[2]]$summary$lci, ymax=covid_IDurInf3[[2]]$summary$uci), alpha=0.2) + xlab('times') + ylab('I')

ggplot(mapping=aes(x=1:428))+geom_density(mapping=aes(y=covid_IDurInf3sc[[2]]$summary$mean), color='#963939', stat='identity') +
 geom_ribbon(aes(ymin=covid_IDurInf3sc[[2]]$summary$lci, ymax=covid_IDurInf3sc[[2]]$summary$uci),fill='#963939', alpha=0.2)+geom_density(mapping=aes(y=covid_IDurInf3sc[[1]]$summary$mean), color='#963939', stat='identity') +
 geom_ribbon(aes(ymin=covid_IDurInf3sc[[1]]$summary$lci, ymax=covid_IDurInf3sc[[1]]$summary$uci),fill='#963939', alpha=0.2)


p1
