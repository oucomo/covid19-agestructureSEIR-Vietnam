library('socialmixr')
library(tidyverse)
library(magrittr)
library(stringr)
library(Hmisc)
setwd("C:\\Users\\nhatlth\\Google Drive\\COvid\\")

### Estimate the population of Vietnam
library(readr)
vietpop<-as.data.frame(read_csv("Viet Nam-2007.csv"))
vietpop %<>% mutate(agegroup=Age,
                    popage=M+F)
total<-sum(vietpop$popage)
vietpop %<>% mutate(propage=popage/total) %>% select(agegroup,popage,propage)
vietpop$agegroup.new<-ifelse(vietpop$agegroup%in%c("5-9","10-14"),"5-14",
                             ifelse(vietpop$agegroup%in%c("15-19","20-24"),"15-24",
                                    ifelse(vietpop$agegroup%in%c("25-29","30-34"),"25-34",
                                           ifelse(vietpop$agegroup%in%c("35-39","40-44","45-49"),"35-49",
                                                  ifelse(vietpop$agegroup%in%c("50-54","55-59","60-64"),"50-64",
                                                         ifelse(vietpop$agegroup%in%c("65-69","70-74","75-79","80-84","85-89","90-94","95-99","100+"),"65+","0-4"))))))

vietpop %<>%group_by(agegroup.new) %>% dplyr::summarise(popage=sum(popage),propage=sum(propage)) %>% ungroup() %>% as.data.frame()
#write.csv(vietpop,file="vietpop.csv")

# Numextract <- function(string){
#   unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
# }
# list_surveys()
# 
# Vietnam_survey <- get_survey("https://doi.org/10.5281/zenodo.1289474")
# saveRDS(Vietnam_survey, "Vietnam.rds")
Vietnam_survey <- readRDS("Vietnam.rds")
contacts<-Vietnam_survey$contacts %>% 
  mutate(cnt_age_group=as.character(cnt_age_group),
         cnt_age_group=ifelse(cnt_age_group=="65+","65-100",cnt_age_group))
contact.age<-as.integer(unlist(strsplit(contacts$cnt_age_group,"-")))
odd.index<-seq(1,length(contact.age),by=2)
even.index<-seq(2,length(contact.age),by=2)
contacts %<>%mutate(cnt_age_est_min=contact.age[odd.index],
                    cnt_age_est_max=contact.age[even.index],
                    cnt_age_exact=round((cnt_age_est_min+cnt_age_est_max)/2)) 
participants<-Vietnam_survey$participants %>% mutate(year=2007,country=as.factor(country))

Vietnam_survey<-survey(participants=participants,
                       contacts=contacts,
                       reference=Vietnam_survey$reference)
check(Vietnam_survey)
clean(Vietnam_survey)
library(geepack)
contacts <-Vietnam_survey$contacts%>%mutate(flag=1) %>% 
            group_by(part_id,cnt_age_group,cnt_home,cnt_work,cnt_school,cnt_transport,cnt_leisure) %>% 
            dplyr::summarise(counts=sum(flag,na.rm = T),frequency_multi=mean(frequency_multi,na.rm = T),
                             phys_contact=mean(phys_contact,na.rm = T),duration_multi=mean(duration_multi,na.rm = T)) %>% ungroup() %>% 
            as.data.frame() %>% 
            mutate(contact_age_group=factor(cnt_age_group,levels = c("0-5","6-15","16-25","26-34","35-49","50-64","65-100"),
                                  labels=c("0-4","5-14","15-24","25-34","35-49","50-64","65+"))) 


participants<-Vietnam_survey$participants %>% 
            mutate(age.group=cut(part_age,breaks=c(-1e-10,5,15,25,35,50,65,100),
                                 labels=c("0-4","5-14","15-24","25-34","35-49","50-64","65+")),
                   part_id=factor(part_id))

dat<-merge(contacts,participants,by="part_id")

# 
# dat2<-as.data.frame(table())
# 
# dat.wide<-


#fit.weight<-glm(phys_contact~cnt_age_group+cnt_home+cnt_work+cnt_school+cnt_transport+)
dat$dayofweek<-as.integer(dat$dayofweek)
dat$dayofweek<-ifelse(is.na(dat$dayofweek),0,dat$dayofweek)
x<-model.matrix(~age.group:contact_age_group-1,data=dat)


# 
# # Method 1
# X.syms.weighted<-as.data.frame(cbind(dat[,c("counts","part_id","hh_id")],"a1_1"=x[,"age.group0-4:contact_age_group0-4"],
#                             "a1_2"=(x[,"age.group5-14:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="5-14"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+x[,"age.group0-4:contact_age_group5-14"]),
#                             "a1_3"=(x[,"age.group15-24:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="15-24"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+x[,"age.group0-4:contact_age_group15-24"]),
#                             "a1_4"=(x[,"age.group25-34:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="25-34"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+x[,"age.group0-4:contact_age_group25-34"]),
#                             "a1_5"=(x[,"age.group35-49:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="35-49"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+x[,"age.group0-4:contact_age_group35-49"]),
#                             "a1_6"=(x[,"age.group50-64:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+x[,"age.group0-4:contact_age_group50-64"]),
#                             "a1_7"=(x[,"age.group65+:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+x[,"age.group0-4:contact_age_group65+"]),
#               
#                             "a2_2"=x[,"age.group5-14:contact_age_group5-14"],
#                             "a2_3"=(x[,"age.group15-24:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="15-24"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+x[,"age.group5-14:contact_age_group15-24"]),
#                             "a2_4"=(x[,"age.group25-34:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="25-34"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+x[,"age.group5-14:contact_age_group25-34"]),
#                             "a2_5"=(x[,"age.group35-49:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="35-49"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+x[,"age.group5-14:contact_age_group35-49"]),
#                             "a2_6"=(x[,"age.group50-64:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+x[,"age.group5-14:contact_age_group50-64"]),
#                             "a2_7"=(x[,"age.group65+:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+x[,"age.group5-14:contact_age_group65+"]),
#               
#                             "a3_3"=x[,"age.group15-24:contact_age_group15-24"],
#                             "a3_4"=(x[,"age.group25-34:contact_age_group15-24"]*vietpop$propage[vietpop$agegroup.new=="25-34"]/vietpop$propage[vietpop$agegroup.new=="15-24"]+x[,"age.group15-24:contact_age_group25-34"]),
#                             "a3_5"=(x[,"age.group35-49:contact_age_group15-24"]*vietpop$propage[vietpop$agegroup.new=="35-49"]/vietpop$propage[vietpop$agegroup.new=="15-24"]+x[,"age.group15-24:contact_age_group35-49"]),
#                             "a3_6"=(x[,"age.group50-64:contact_age_group15-24"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="15-24"]+x[,"age.group15-24:contact_age_group50-64"]),
#                             "a3_7"=(x[,"age.group65+:contact_age_group15-24"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="15-24"]+x[,"age.group15-24:contact_age_group65+"]),
#               
#               
#                             "a4_4"=x[,"age.group25-34:contact_age_group25-34"],
#                             "a4_5"=(x[,"age.group35-49:contact_age_group25-34"]*vietpop$propage[vietpop$agegroup.new=="35-49"]/vietpop$propage[vietpop$agegroup.new=="25-34"]+x[,"age.group25-34:contact_age_group35-49"]),
#                             "a4_6"=(x[,"age.group50-64:contact_age_group25-34"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="25-34"]+x[,"age.group25-34:contact_age_group50-64"]),
#                             "a4_7"=(x[,"age.group65+:contact_age_group25-34"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="25-34"]+x[,"age.group25-34:contact_age_group65+"]),
#               
#               
#                             "a5_5"=x[,"age.group35-49:contact_age_group35-49"],
#                             "a5_6"=(x[,"age.group50-64:contact_age_group35-49"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="35-49"]+x[,"age.group35-49:contact_age_group50-64"]),
#                             "a5_7"=(x[,"age.group65+:contact_age_group35-49"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="35-49"]+x[,"age.group35-49:contact_age_group65+"]),
#               
#                             "a6_6"=x[,"age.group50-64:contact_age_group50-64"],
#                             "a6_7"=(x[,"age.group65+:contact_age_group50-64"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="50-64"]+x[,"age.group50-64:contact_age_group65+"]),
#               
#                             "a7_7"=x[,"age.group65+:contact_age_group65+"],
#                             "weight_hh_size"=1/x[,"hh_size"],
#                             "hh_size"=dat[,"hh_size"],
#                             "type_of_day"=ifelse(dat[,"dayofweek"]>4,"working_day","non_working")
#               ))
# fit<-geeglm(counts~a1_1+a1_2+a1_3+a1_4+a1_5+a1_6+a1_7+
#                    a2_2+a2_3+a2_4+a2_5+a2_6+a2_7+
#                    a3_3+a3_4+a3_5+a3_6+a3_7+
#                    a4_4+a4_5+a4_6+a4_7+
#                    a5_5+a5_6+a5_7+
#                    a6_6+a6_7+a7_7-1,id=part_id,family=poisson(link = "log"),weights = weight_hh_size,data=X.syms.weighted,corstr = "exchangeable")
# 
# 
# 
# contact_matrix<- matrix(0, 7, 7)
# contact_matrix[upper.tri(contact_matrix, diag=TRUE)]<-exp(unname(fit$coefficients)[1:28])
# contact_matrix <- contact_matrix + t(contact_matrix) - diag(diag(contact_matrix))
# 

# Method 2
vietpop$age.group<-vietpop$agegroup.new
dat<-merge(dat,vietpop,by="age.group")
dat<-arrange(dat,hh_id,part_id)
fit<-geeglm(counts~age.group:contact_age_group-1,id=hh_id,family=poisson(link = "log"),weights = popage,data=dat,corstr = "exchangeable")

contact_matrix<- matrix(exp(unname(fit$coefficients)), 7, 7)
vector<-c("a1_1"=fit$coefficients["age.group0-4:contact_age_group0-4"],
                            "a1_2"=0.5*(fit$coefficients["age.group5-14:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="5-14"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+fit$coefficients["age.group0-4:contact_age_group5-14"]),
                            "a1_3"=0.5*(fit$coefficients["age.group15-24:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="15-24"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+fit$coefficients["age.group0-4:contact_age_group15-24"]),
                            "a1_4"=0.5*(fit$coefficients["age.group25-34:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="25-34"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+fit$coefficients["age.group0-4:contact_age_group25-34"]),
                            "a1_5"=0.5*(fit$coefficients["age.group35-49:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="35-49"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+fit$coefficients["age.group0-4:contact_age_group35-49"]),
                            "a1_6"=0.5*(fit$coefficients["age.group50-64:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+fit$coefficients["age.group0-4:contact_age_group50-64"]),
                            "a1_7"=0.5*(fit$coefficients["age.group65+:contact_age_group0-4"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="0-4"]+fit$coefficients["age.group0-4:contact_age_group65+"]),

                            "a2_2"=fit$coefficients["age.group5-14:contact_age_group5-14"],
                            "a2_3"=0.5*(fit$coefficients["age.group15-24:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="15-24"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+fit$coefficients["age.group5-14:contact_age_group15-24"]),
                            "a2_4"=0.5*(fit$coefficients["age.group25-34:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="25-34"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+fit$coefficients["age.group5-14:contact_age_group25-34"]),
                            "a2_5"=0.5*(fit$coefficients["age.group35-49:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="35-49"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+fit$coefficients["age.group5-14:contact_age_group35-49"]),
                            "a2_6"=0.5*(fit$coefficients["age.group50-64:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+fit$coefficients["age.group5-14:contact_age_group50-64"]),
                            "a2_7"=0.5*(fit$coefficients["age.group65+:contact_age_group5-14"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="5-14"]+fit$coefficients["age.group5-14:contact_age_group65+"]),

                            "a3_3"=fit$coefficients["age.group15-24:contact_age_group15-24"],
                            "a3_4"=0.5*(fit$coefficients["age.group25-34:contact_age_group15-24"]*vietpop$propage[vietpop$agegroup.new=="25-34"]/vietpop$propage[vietpop$agegroup.new=="15-24"]+fit$coefficients["age.group15-24:contact_age_group25-34"]),
                            "a3_5"=0.5*(fit$coefficients["age.group35-49:contact_age_group15-24"]*vietpop$propage[vietpop$agegroup.new=="35-49"]/vietpop$propage[vietpop$agegroup.new=="15-24"]+fit$coefficients["age.group15-24:contact_age_group35-49"]),
                            "a3_6"=0.5*(fit$coefficients["age.group50-64:contact_age_group15-24"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="15-24"]+fit$coefficients["age.group15-24:contact_age_group50-64"]),
                            "a3_7"=0.5*(fit$coefficients["age.group65+:contact_age_group15-24"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="15-24"]+fit$coefficients["age.group15-24:contact_age_group65+"]),


                            "a4_4"=fit$coefficients["age.group25-34:contact_age_group25-34"],
                            "a4_5"=0.5*(fit$coefficients["age.group35-49:contact_age_group25-34"]*vietpop$propage[vietpop$agegroup.new=="35-49"]/vietpop$propage[vietpop$agegroup.new=="25-34"]+fit$coefficients["age.group25-34:contact_age_group35-49"]),
                            "a4_6"=0.5*(fit$coefficients["age.group50-64:contact_age_group25-34"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="25-34"]+fit$coefficients["age.group25-34:contact_age_group50-64"]),
                            "a4_7"=0.5*(fit$coefficients["age.group65+:contact_age_group25-34"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="25-34"]+fit$coefficients["age.group25-34:contact_age_group65+"]),


                            "a5_5"=fit$coefficients["age.group35-49:contact_age_group35-49"],
                            "a5_6"=0.5*(fit$coefficients["age.group50-64:contact_age_group35-49"]*vietpop$propage[vietpop$agegroup.new=="50-64"]/vietpop$propage[vietpop$agegroup.new=="35-49"]+fit$coefficients["age.group35-49:contact_age_group50-64"]),
                            "a5_7"=0.5*(fit$coefficients["age.group65+:contact_age_group35-49"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="35-49"]+fit$coefficients["age.group35-49:contact_age_group65+"]),

                            "a6_6"=fit$coefficients["age.group50-64:contact_age_group50-64"],
                            "a6_7"=0.5*(fit$coefficients["age.group65+:contact_age_group50-64"]*vietpop$propage[vietpop$agegroup.new=="65+"]/vietpop$propage[vietpop$agegroup.new=="50-64"]+fit$coefficients["age.group50-64:contact_age_group65+"]),

                            "a7_7"=fit$coefficients["age.group65+:contact_age_group65+"])

contact_matrix_sym<- matrix(0, 7, 7)
contact_matrix_sym[lower.tri(contact_matrix_sym, diag=TRUE)]<-exp(vector)
contact_matrix_sym <- contact_matrix_sym + t(contact_matrix_sym) - diag(diag(contact_matrix_sym))

library(ggplot2)
library(dplyr)
library(tidyr)
#df1<-(contact_matrix)
df1<-(contact_matrix_sym)
colnames(df1)<-c("0-4","5-14","15-24","25-34","35-49","50-64","65+")
rownames(df1)<-c("0-4","5-14","15-24","25-34","35-49","50-64","65+")
plotDat <- as.data.frame(reshape2::melt(df1))
colnames(plotDat)<-c("contact_group","participant_group","intensity")
plotDat$intensity<-round(plotDat$intensity,2)
ggplot(plotDat, aes(participant_group,contact_group, col = intensity, fill = intensity, label = intensity)) +
  geom_tile() +
  geom_text(col = "black") +
  theme_minimal() +
  scale_fill_gradient(low = "darkblue", high ="greenyellow") +
  scale_color_gradient(low = "darkblue", high ="greenyellow")

write.csv(contact_matrix,file="Data/contact_matrix.csv")
