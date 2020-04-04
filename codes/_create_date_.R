date <- as.Date(as.Date('2020-3-10'):(as.Date('2020-3-10')+100), origin='1970-1-1')
#26 cases before 10
write.csv(date,'data/vn_date.csv')