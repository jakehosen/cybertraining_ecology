install.packages("tidyr")
install.packages("dplyr")# clean data
install.packages("tidyr")# clean data
install.packages("rLakeAnalyzer")# lake analyses
install.packages("lubridate")# working with time
install.packages("LakeMetabolizer")# lake analyses


library(tidyr)
library(dplyr)# clean data
library(tidyr)# clean data
library(rLakeAnalyzer)# lake analyses
library(lubridate)# working with time
library(LakeMetabolizer)# lake analyses
library(curl)
library(ggplot2)
library(readr)


bass_do<-readRDS("~/metabolism_modeling/BASS_do_compiled_20210729-20220712.rds")

bass_temp<-readRDS("~/BASS_temperature_compiled_20210729-20220712.rds")



crooked_temp<-readRDS("~/metabolism_modeling/CROK_temperature_compiled_20220606-20230808.rds")

crooked_do<-readRDS("~/metabolism_modeling/CROK_do_compiled_20220606-20230808.rds")



failing_temp<-readRDS("~/metabolism_modeling/FAIL_temperature_compiled_20220412-20230703.rds")

failing_do<-readRDS("~/metabolism_modeling/FAIL_do_compiled_20220412-20230703.rds")


citizen_temp<-readRDS("~/metabolism_modeling/CITZ_temperature_compiled_20220615-20221119.rds")
citizen_do<-readRDS("~/metabolism_modeling/CITZ_do_compiled_20220615-20221119.rds")


wtr.heat.map(wtr=as.data.frame(bass_temp))
wtr.plot.temp(wtr=as.data.frame(bass_temp),na.rm=T)

wtr.heat.map(wtr=as.data.frame(crooked_temp))
wtr.plot.temp(wtr=as.data.frame(crooked_temp),na.rm=T)

wtr.heat.map(wtr=as.data.frame(failing_temp))
wtr.plot.temp(wtr=as.data.frame(failing_temp),na.rm=T)

wtr.heat.map(wtr=as.data.frame(citizen_temp))
wtr.plot.temp(wtr=as.data.frame(citizen_do),na.rm=T)



