


'4.结果'
'4.1.MODS'
'这是MODS的代码,供参考
end_bc_sofa <- sofa_bc %>%  
  group_by(stay_id) %>%  
  summarise(end_sofa=max(sofa_24hours), mods_endtime=sofa_charttime[which.max(sofa_24hours)])

end_nobc_sofa <- sofa_nobc %>%  
  group_by(stay_id) %>%  
  summarise(end_sofa=max(sofa_24hours), mods_endtime=sofa_charttime[which.max(sofa_24hours)])'
#划定亚组
'统一一下,sub1都是大于,sub2都是小于'
'注:修改前的newdata_sub1在newdata_sub1_save中,在外面也有储存'
write.csv(newdata,file = 'newdata_亚组前保存_2月29日.csv')
densityplot(newdata_sub1$bc_base,xlim = c(5,30))
newdata = subset(newdata_亚组前保存_2月29日, glucose_base >= 28)#做完了打个叉,20242月26 23:33开始做,做的就把保险栓先拉掉
newdata = subset(newdata_亚组前保存_2月29日, potassium_base > 5)
newdata_sub1
densityplot(newdata_sub1$sofar)

densityplot(newdata_sub1$glucose_base,xlim = c(0,500))
newdata_sub1_sub1 = subset(newdata_sub1, glucose >2)
newdata_sub1_sub1 = subset(newdata_sub1, glucose >2)
densityplot(newdata_sub1$bc_base)

densityplot(newdata_sub1$)

densityplot(newdata_sub1$chloride_base,xlim = c(90,120))

densityplot(newdata_sub1$ph_base,xlim = c(7.1,7.5))

densityplot(newdata_sub1$calcium_base,xlim = c(6,11))

densityplot(newdata_sub1$sirs)


newdata_sub1$dmods=newdata_sub1$end_sofa-newdata_sub1$sofat
newdata_sub1$timetomods=difftime(newdata_sub1$mods_endtime,newdata_sub1$admittime,units = "days")
newdata_sub1$mods= ifelse(newdata_sub1$dmods >= 2, 1,0)
newdata_sub1$timetomods=ifelse(newdata_sub1$mods == 0,"",newdata_sub1$timetomods)
newdata_sub1$timetomods=as.numeric(newdata_sub1$timetomods)#列的类型转化
newdata_sub1$mods=as.numeric(newdata_sub1$mods)#列的类型转化
surv_mods_sub1=survfit(Surv(newdata_sub1$timetomods,newdata_sub1$mods)~bc_use,data = newdata_sub1)
'图片需要美化'
ggsurvplot(surv_mods_sub1,data = newdata_sub1,conf.int = T,pval=T,
           xlim = c(0,7),
           ggtheme = theme_bw(base_size=7),
           pval.coord = c(5.6, 0.1),
           pval.size = 3,
           legend.title = '',
           legend = 'top',
           legend.labs = c("Control group", "Bicarbonate group"),
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           xlab = 'Time (days)',
           font.main = "Arial",
           ylab = 'SOFA exceeding threshold (%)',
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           palette =c("#E7B800", "#2E9FDF"),
           break.x.by = 2,
           fun = 'event',
)
cox_mods=coxph(Surv(newdata_sub1$timetomods,newdata_sub1$mods)~bc_use,data = newdata_sub1)
autoReg(cox_mods)
cox_mods=coxph(Surv(newdata_sub1$timetomods,newdata_sub1$mods)~bc_use+#分组变量
                 gender+admission_age+type+sirs+ag_base+ph_base+glucose_base+sodium_base+potassium_base+chloride_base+calcium_base+
                 bc_base+bun_base+creatinine_base+gcs+sofat+sofar+kdigo#协变量
               ,data = newdata_sub1)
autoReg(cox_mods)
''

mods_check = subset(newdata_sub1, mods == 1)




'4.2.肌酐,肾脏指标1'



#肌酐上升,阈值和老师讨论一下

newdata_sub1$timetocreat=difftime(newdata_sub1$creat_endtime,newdata_sub1$admittime,units = 'days')
newdata_sub1$dcreat=newdata_sub1$end_creat-newdata_sub1$creatinine_base
newdata_sub1$creat= ifelse(newdata_sub1$dcreat >= 0, 1,0)
newdata_sub1$timetocreat=ifelse(newdata_sub1$creat == 0,"",newdata_sub1$timetocreat)
newdata_sub1$timetocreat=as.numeric(newdata_sub1$timetocreat)#列的类型转化
newdata_sub1$creatinine_base =as.numeric(newdata_sub1$creatinine_base)#列的类型转化
surv_creat=survfit(Surv(newdata_sub1$timetocreat,newdata_sub1$creat)~bc_use,data = newdata_sub1)
ggsurvplot(surv_creat,data = newdata_sub1,conf.int = T,pval=T,
           xlim = c(0,7),
           ggtheme = theme_bw(),
           pval.coord = c(5.3, 0.7),
           legend.title = '',
           legend = 'top',
           legend.labs = c("Control group", "Bicarbonate group"),
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           xlab = 'Time (days)',
           font.main = "Arial",
           ylab = 'Cumulative rate (percentage)',
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           palette = c("#8BA583", "#CDA59E"),
           break.x.by = 2,
           fun = event
)
cox_creat=coxph(Surv(newdata_sub1$timetocreat,newdata_sub1$creat)~bc_use+gender+admission_age+dialysis+oasis+electivesurgery+sapsii+sirs+vasopressin+
                  vent_use+ag_base+glucose_base+sodium_base+potassium_base+chloride_base+
                  bc_base+bun_base+creatinine_base+gcs+sofat+sofar+sofac+uo+kdigo,data = newdata_sub1,weights = newdata_sub1$ipw)
summary(cox_creat)
cox_creat=coxph(Surv(newdata_sub1$timetocreat,newdata_sub1$creat)~bc_use,data = newdata_sub1)
autoReg(cox_creat)
ggsurvplot(surv_creat,data = newdata_sub1_sub1,title='creat',conf.int=T,pval=T,xlim=c(0,10000))
'4.3肾脏SOFA,这是指标2'
newdata_sub1$drenal=newdata_sub1$end_sofar-newdata_sub1$sofar
newdata_sub1$timetorenal=difftime(newdata_sub1$sofa_endtime,newdata_sub1$admittime,units = 'days')
newdata_sub1$renal= ifelse(newdata_sub1$drenal >= 1, 1,0)
newdata_sub1$timetorenal=ifelse(newdata_sub1$renal == 0,"",newdata_sub1$timetorenal)
newdata_sub1$timetorenal=as.numeric(newdata_sub1$timetorenal)#列的类型转化
newdata_sub1$renal=as.numeric(newdata_sub1$renal)#列的类型转化
surv_renal_sub1=survfit(Surv(newdata_sub1$timetorenal,newdata_sub1$renal)~bc_use,data = newdata_sub1)
ggsurvplot(surv_renal_sub1,data = newdata_sub1,title='sofa_renal',conf.int=T,xlim=c(0,7),pval=T)

cox_renal=coxph(Surv(newdata_sub1$timetorenal,newdata_sub1$renal)~bc_use+gender+admission_age+dialysis+oasis+electivesurgery+sapsii+sirs+vasopressin+
                  vent_use+ag_base+glucose_base+sodium_base+potassium_base+chloride_base+
                  bc_base+bun_base+creatinine_base+gcs+sofat+sofar+sofac+uo+kdigo,data = newdata_sub1,weights = newdata_sub1$distance)
cox_renal_sub1=coxph(Surv(newdata_sub1$timetorenal,newdata_sub1$renal)~bc_use,data = newdata_sub1)
autoReg(cox_renal_sub1)


'4.4.神经功能指标'
newdata_sub1$dgcs=newdata_sub1$end_gcs-newdata_sub1$gcs
newdata_sub1$timetogcs=difftime(newdata_sub1$gcs_endtime,newdata_sub1$admittime,units = 'days')
newdata_sub1$egcs= ifelse(newdata_sub1$dgcs <= -1, 1,0)
newdata_sub1$timetogcs=ifelse(newdata_sub1$egcs == 0,"",newdata_sub1$timetogcs)
newdata_sub1$timetogcs=as.numeric(newdata_sub1$timetogcs)#列的类型转化
newdata_sub1$egcs=as.numeric(newdata_sub1$egcs)#列的类型转化
gcs_check = subset(newdata_sub1, egcs == 1)
surv_gcs=survfit(Surv(newdata_sub1$timetogcs,newdata_sub1$egcs)~bc_use,data = newdata_sub1)
cox_gcs=coxph(Surv(newdata_sub1$timetogcs,newdata_sub1$egcs)~bc_use+
                ag_base+glucose_base+sodium_base+potassium_base+chloride_base+
                bc_base+bun_base+creatinine_base+gcs+sofat+sofar+sofac,data = newdata_sub1,weights=newdata_sub1$ipw)
summary(cox_gcs)
ggsurvplot(surv_gcs,data = newdata_sub1,title='gcs',conf.int=T,xlim=c(0,7),pval=T)
density_plot <- density(newdata_sub1$glucose_base)
plot(density_plot,xlim = c(0,600))
newdata_sub1_sub1 = subset(newdata_sub1,glucose_base > 3)
surv_gcs=survfit(Surv(newdata_sub1_sub1$timetogcs,newdata_sub1_sub1$egcs)~bc_use,data = newdata_sub1_sub1)
ggsurvplot(surv_gcs,data = newdata_sub1_sub1,title='gcs',conf.int=T,xlim=c(0,10000),pval=T)
cox_gcs_sub1=coxph(Surv(newdata_sub1$timetogcs,newdata_sub1$egcs)~bc_use,data = newdata_sub1)
autoReg(cox_gcs)
#KDIGO评分

'4.4.二次入院'
#二次入院,做logistic回归
sec = newdata_sub1[,c('stay_id','admittime')]
second_admission = second_admission %>%
  distinct()
newdata_sub1_sec = newdata_sub1 %>%
  left_join(second_admission,by = 'stay_id')
newdata_sub1_sec$first_icu_stay = ifelse(newdata_sub1_sec$first_icu_stay == 'f', 1, 0)
newdata_sub1_sec$first_icu_stay = ifelse(is.na(newdata_sub1_sec$first_icu_stay),0,1)
sec_model = glm(first_icu_stay~bc_use,data = newdata_sub1_sec)
autoReg(sec_model)

newdata_sub1_sec$secadmit = ifelse(is.na(newdata_sub1$secadmit),0,1)
second = glm(secadmit ~ bc_use,data = newdata_sub1)
summary(second)
'4.5.死亡'
#死亡,死亡没有统计学差异,因此可以放在附件当中
mortality = df_merge[,c('stay_id','hospital_expire_flag','dod')]
newdata_sub1_death = newdata_sub1 %>%
  inner_join(mortality,by = 'stay_id')
newdata_sub1$death=newdata_sub1_death$hospital_expire_flag
newdata_sub1$timetodeath=newdata_sub1_death$dod - newdata_sub1_death$admittime
newdata_sub1$timetodeath=ifelse(newdata_sub1$death == 0,"",newdata_sub1$timetodeath)
newdata_sub1$timetodeath=as.numeric(newdata_sub1$timetodeath)#列的类型转化
newdata_sub1$death=as.numeric(newdata_sub1$death)#列的类型转化
surv_death_sub1=survfit(Surv(newdata_sub1$timetodeath,newdata_sub1$death)~bc_use,data = newdata_sub1)
cox_death_sub1=coxph(Surv(newdata_sub1$timetodeath,newdata_sub1$death)~bc_use,data = newdata_sub1)
death_model_sub1 = glm(death~bc_use,data = newdata_sub1)
autoReg(death_model_sub1)
ggsurvplot(surv_death,data = newdata_sub1,title='death',conf.int=T,pval = T)

'4.6钙盐钾盐补液次数分析'
potassium_model_sub1 = glm(n_potassium~bc_use,data = newdata_sub1)
autoReg(potassium_model_sub1)
calcium_model_sub1 = glm(n_calcium~bc_use,data = newdata_sub1)
autoReg(calcium_model_sub1)
autoReg(cox_mods)
'MSCM'

#**********************正式稿件开始
'先更新一下list'
newdata_sub1_list=newdata_sub1[,c("stay_id","admittime")]
names(newdata_sub1_list)=c("stay_id","admittime")
'这是正式定稿'

'钾离子'
potassium_mscm_sub1 = potassium %>%
  arrange(stay_id,potassium_time) %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(potassium_time > admittime)

potassium_t0 = potassium_mscm_sub1 %>%
  group_by(stay_id) %>%
  summarise(to = min(potassium_time)) %>%
  distinct()

potassium_mscm_sub1 = potassium_mscm_sub1 %>%
  left_join(potassium_t0,by = 'stay_id') %>%
  distinct()
potassium_mscm_sub1$tstart = difftime(potassium_mscm_sub1$potassium_time,potassium_mscm_sub1$to,units = "days")

potassium_mscm_sub1$tstart = floor(potassium_mscm_sub1$tstart)

potassium_mscm_sub1 = potassium_mscm_sub1 %>%
  group_by(stay_id) %>%
  mutate(tstop = tstart+1) %>%
  filter(!is.na(tstop)) #这个代码可以一键出tstart和tstop

potassium_mscm_sub1_day = potassium_mscm_sub1 %>%
  group_by(stay_id,tstart) %>%
  summarise(max_potassium = max(potassium)) %>%
  distinct()

'血糖'
glucose_mscm_sub1 = glucose %>%
  arrange(stay_id,glucose_time) %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(glucose_time > admittime)

glucose_t0 = glucose_mscm_sub1 %>%
  group_by(stay_id) %>%
  summarise(to = min(glucose_time)) %>%
  distinct()

glucose_mscm_sub1 = glucose_mscm_sub1 %>%
  left_join(glucose_t0,by = 'stay_id') %>%
  distinct()
glucose_mscm_sub1$tstart = difftime(glucose_mscm_sub1$glucose_time,glucose_mscm_sub1$to,units = "days")

glucose_mscm_sub1$tstart = floor(glucose_mscm_sub1$tstart)

glucose_mscm_sub1 = glucose_mscm_sub1 %>%
  group_by(stay_id) %>%
  mutate(tstop = tstart+1) %>%
  filter(!is.na(tstop)) #这个代码可以一键出tstart和tstop

glucose_mscm_sub1_day = glucose_mscm_sub1 %>%
  group_by(stay_id,tstart) %>%
  summarise(max_glucose = max(glucose)) %>%
  distinct()



'AG'
ag_mscm_sub1 = ag %>%
  arrange(stay_id,ag_time) %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(ag_time > admittime)

ag_t0 = ag_mscm_sub1 %>%
  group_by(stay_id) %>%
  summarise(to = min(ag_time)) %>%
  distinct()

ag_mscm_sub1 = ag_mscm_sub1 %>%
  left_join(ag_t0,by = 'stay_id') %>%
  distinct()
ag_mscm_sub1$tstart = difftime(ag_mscm_sub1$ag_time,ag_mscm_sub1$to,units = "days")

ag_mscm_sub1$tstart = floor(ag_mscm_sub1$tstart)

ag_mscm_sub1 = ag_mscm_sub1 %>%
  group_by(stay_id) %>%
  mutate(tstop = tstart+1) %>%
  filter(!is.na(tstop)) #这个代码可以一键出tstart和tstop

ag_mscm_sub1_day = ag_mscm_sub1 %>%
  group_by(stay_id,tstart) %>%
  summarise(max_ag = max(ag)) %>%
  distinct()

'pH'
ph_mscm_sub1 = ph %>%
  arrange(stay_id,ph_time) %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(ph_time > admittime)

ph_t0 = ph_mscm_sub1 %>%
  group_by(stay_id) %>%
  summarise(to = min(ph_time)) %>%
  distinct()

ph_mscm_sub1 = ph_mscm_sub1 %>%
  left_join(ph_t0,by = 'stay_id') %>%
  distinct()
ph_mscm_sub1$tstart = difftime(ph_mscm_sub1$ph_time,ph_mscm_sub1$to,units = "days")

ph_mscm_sub1$tstart = floor(ph_mscm_sub1$tstart)

ph_mscm_sub1 = ph_mscm_sub1 %>%
  group_by(stay_id) %>%
  mutate(tstop = tstart+1) %>%
  filter(!is.na(tstop)) #这个代码可以一键出tstart和tstop

ph_mscm_sub1_day = ph_mscm_sub1 %>%
  group_by(stay_id,tstart) %>%
  summarise(min_ph = min(ph)) %>%
  distinct()

'氯离子'
chloride_mscm_sub1 = chloride %>%
  arrange(stay_id,chloride_time) %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(chloride_time > admittime)

chloride_t0 = chloride_mscm_sub1 %>%
  group_by(stay_id) %>%
  summarise(to = min(chloride_time)) %>%
  distinct()

chloride_mscm_sub1 = chloride_mscm_sub1 %>%
  left_join(chloride_t0,by = 'stay_id') %>%
  distinct()
chloride_mscm_sub1$tstart = difftime(chloride_mscm_sub1$chloride_time,chloride_mscm_sub1$to,units = "days")

chloride_mscm_sub1$tstart = floor(chloride_mscm_sub1$tstart)

chloride_mscm_sub1 = chloride_mscm_sub1 %>%
  group_by(stay_id) %>%
  mutate(tstop = tstart+1) %>%
  filter(!is.na(tstop)) #这个代码可以一键出tstart和tstop

chloride_mscm_sub1_day = chloride_mscm_sub1 %>%
  group_by(stay_id,tstart) %>%
  summarise(max_chloride = max(chloride)) %>%
  distinct()

'钙离子'
calcium_mscm_sub1 = calcium %>%
  arrange(stay_id,calcium_time) %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(calcium_time > admittime)

calcium_t0 = calcium_mscm_sub1 %>%
  group_by(stay_id) %>%
  summarise(to = min(calcium_time)) %>%
  distinct()

calcium_mscm_sub1 = calcium_mscm_sub1 %>%
  left_join(calcium_t0,by = 'stay_id') %>%
  distinct()
calcium_mscm_sub1$tstart = difftime(calcium_mscm_sub1$calcium_time,calcium_mscm_sub1$to,units = "days")

calcium_mscm_sub1$tstart = floor(calcium_mscm_sub1$tstart)

calcium_mscm_sub1 = calcium_mscm_sub1 %>%
  group_by(stay_id) %>%
  mutate(tstop = tstart+1) %>%
  filter(!is.na(tstop)) #这个代码可以一键出tstart和tstop

calcium_mscm_sub1_day = calcium_mscm_sub1 %>%
  group_by(stay_id,tstart) %>%
  summarise(min_calcium = min(calcium)) %>%
  distinct()

'碳酸氢根浓度'
bc_mscm_sub1 = bc %>%
  arrange(stay_id,bc_time) %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(bc_time > admittime)

bc_t0 = bc_mscm_sub1 %>%
  group_by(stay_id) %>%
  summarise(to = min(bc_time)) %>%
  distinct()

bc_mscm_sub1 = bc_mscm_sub1 %>%
  left_join(bc_t0,by = 'stay_id') %>%
  distinct()
bc_mscm_sub1$tstart = difftime(bc_mscm_sub1$bc_time,bc_mscm_sub1$to,units = "days")

bc_mscm_sub1$tstart = floor(bc_mscm_sub1$tstart)

bc_mscm_sub1 = bc_mscm_sub1 %>%
  group_by(stay_id) %>%
  mutate(tstop = tstart+1) %>%
  filter(!is.na(tstop)) #这个代码可以一键出tstart和tstop

bc_mscm_sub1_day = bc_mscm_sub1 %>%
  group_by(stay_id,tstart) %>%
  summarise(min_bc = min(bc)) %>%
  distinct()



'creat'
creat_mscm_sub1 = creatinine %>%
  arrange(stay_id,creat_time) %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(creat_time > admittime)

creat_t0 = creat_mscm_sub1 %>%
  group_by(stay_id) %>%
  summarise(to = min(creat_time)) %>%
  distinct()

creat_mscm_sub1 = creat_mscm_sub1 %>%
  left_join(creat_t0,by = 'stay_id') %>%
  distinct()
creat_mscm_sub1$tstart = difftime(creat_mscm_sub1$creat_time,creat_mscm_sub1$to,units = "days")

creat_mscm_sub1$tstart = floor(creat_mscm_sub1$tstart)

creat_mscm_sub1 = creat_mscm_sub1 %>%
  group_by(stay_id) %>%
  mutate(tstop = tstart+1) %>%
  filter(!is.na(tstop)) #这个代码可以一键出tstart和tstop

creat_mscm_sub1_day = creat_mscm_sub1 %>%
  group_by(stay_id,tstart) %>%
  summarise(max_creat = max(creat)) %>%
  distinct()



'bc使用'
bc_merge=bc_time_vary %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(charttime > admittime)
bc_infusion=bc_merge[,c("stay_id","charttime")]
names(bc_infusion)=c("stay_id","charttime")
bc_infusion$bc_use=1

infusion_mscm_sub1 = bc_time_vary %>%
  arrange(stay_id,charttime) %>%
  left_join (newdata_sub1_list, by='stay_id') %>%
  filter(charttime > admittime)
infusion_t0 = infusion_mscm_sub1 %>%
  group_by(stay_id) %>%
  summarise(to = min(charttime)) %>%
  distinct()
infusion_mscm_sub1 = infusion_mscm_sub1 %>%
  left_join(infusion_t0,by = 'stay_id') %>%
  distinct()
infusion_mscm_sub1$tstart = difftime(infusion_mscm_sub1$charttime,infusion_mscm_sub1$to,units = "days")
infusion_mscm_sub1$tstart = floor(infusion_mscm_sub1$tstart)
infusion_mscm_sub1 = unique(infusion_mscm_sub1)
infusion_mscm_sub1$bc_use=1
bc_infusion = infusion_mscm_sub1[,c("stay_id","tstart","bc_use")]

'以上合并了时依变量'
''

'合并变量'
mscm_time = Reduce(function(x, y) merge(x, y, by = c("stay_id","tstart"),all = T),
                   list(bc_mscm_sub1_day,ag_mscm_sub1_day,ph_mscm_sub1_day,chloride_mscm_sub1_day,calcium_mscm_sub1_day,creat_mscm_sub1_day,
                        glucose_mscm_sub1_day,potassium_mscm_sub1_day,bc_infusion))
mscm_time = mscm_time %>%
  distinct()
'合并'
mscm_fix = newdata_sub1[,c("stay_id","sapsii","admission_age","vasopressin","sirs","dialysis","electivesurgery","gender")]
mscm = Reduce(function(x, y) merge(x, y, by = c("stay_id"),all = T),
              list(mscm_time,mscm_fix))

#***********参照这个代码来安排tstart和tstop
mscm <- mscm %>%
  group_by(stay_id) %>%
  mutate(tstop = tstart+1) %>%
  filter(!is.na(tstop)) 
'拆分成若干个主题,然后一个一个做'


#************
'creat的MSCM结局,如果creat定义改掉了,就从这里开始往下就行了'
termin_creat = newdata_sub1[,c('stay_id','timetocreat','creat')]
names(termin_creat) = c('stay_id',"tstop","end_creat")
termin_creat$tstop = ceiling(termin_creat$tstop)
termin_creat = termin_creat %>%
  filter(end_creat == 1)
mscm_creat = mscm
mscm_creat$bc_label = 0
mscm_creat$fail = 0

for (i in 1:nrow(termin_creat)) {
  for (j in 1:nrow(mscm_creat)){
    if(termin_creat$stay_id[i] == mscm_creat$stay_id[j] & termin_creat$tstop[i] == mscm_creat$tstop[j])
    {mscm_creat$fail[j] = 1}
  }
}

"上面代码合并碳酸氢钠使用和结局时间"  

while(TRUE) {
  namelist_creat = unique(mscm$stay_id)
  data_filtered = mscm_creat[(0:0),]
  for (i in namelist_creat) {
    totide = mscm_creat %>%
      filter(stay_id == i)
    index = which(totide$fail == 1)[1]
    if (!is.na(index)){
      cleaned = totide[seq_len(index),]
      data_filtered = rbind(data_filtered,cleaned)}
    else{
      data_filtered = rbind(data_filtered,totide)}
  }
  rm(totide)
  rm(cleaned)
  break
}
data_filtered = data_filtered %>%
  mutate_at(vars("min_bc":"max_potassium"),as.numeric)
data_filtered
set.seed(666)
'mice插补'
df_input = kNN(data_filtered[,c('min_bc','max_ag','min_ph','max_glucose','min_calcium','max_chloride','max_potassium','max_creat')])
newdata_sub1_mscm_sub1_creat = cbind(complete(df_input),data_filtered)

'生存分析'
creat_final = newdata_sub1_mscm_sub1_creat[,-c(19:26)]
creat_final$bc_label = ifelse(is.na(creat_final$bc_use) == T, 0, 1)
head(creat_final)
temp_creat = ipwtm(exposure = bc_label,family = "survival",
                   numerator = ~sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender,
                   denominator = ~ sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender+
                     max_ag+max_glucose+min_ph+min_bc+max_chloride+min_calcium+max_potassium+max_creat,
                   id = stay_id, tstart = tstart, timevar = tstop, type = "first", data = creat_final)
temp_gcs_sub1 = ipwtm(exposure = bc_label,family = "survival",
                      denominator = ~ ag+glucose,
                      id = stay_id, tstart = start, timevar = stop, type = "first", data = gcs_final_sub1)
surv_creat_mscm_sub1=survfit(Surv(creat_final$tstart,creat_final$fail)~bc_label,data = creat_final,weights = temp_creat$ipw.weights)
ggsurvplot(surv_creat_mscm_sub1,data = creat_final,title='creat_mscm_sub1',conf.int = T,xlim = c(0,20))
rm(cox_gcs_mscm_sub1)
cox_creat_mscm_sub1 = coxph(Surv(tstart,tstop,fail)~bc_label+min_calcium+
                         max_ag+max_glucose+min_bc+max_chloride+max_potassium+max_creat,data = creat_final,weights = temp_creat$ipw.weights)
summary(cox_creat_mscm_sub1)
cox_creat_mscm_sub1 = coxph(Surv(tstart,tstop,fail)~bc_label,data = creat_final,weights = temp_creat$ipw.weights)
cox_creat_mscm_sub1 = coxph(Surv(tstart,tstop,fail)~bc_label,data = creat_final,weights = temp_creat$ipw.weights)
autoReg(cox_creat_mscm_sub1)
write.csv(autoReg(cox_creat_mscm_sub1),file = 'mscm_creat.csv')







#************
'MODS的MSCM结局'
termin_mods = newdata_sub1[,c('stay_id','timetomods','mods')]
names(termin_mods) = c('stay_id',"tstop","end_mods")
termin_mods$tstop = ceiling(termin_mods$tstop)
termin_mods = termin_mods %>%
  filter(end_mods == 1)
mscm_mods = mscm
mscm_mods$bc_label = 0
mscm_mods$fail = 0

for (i in 1:nrow(termin_mods)) {
  for (j in 1:nrow(mscm_mods)){
    if(termin_mods$stay_id[i] == mscm_mods$stay_id[j] & termin_mods$tstop[i] == mscm_mods$tstop[j])
    {mscm_mods$fail[j] = 1}
  }
}

"上面代码合并碳酸氢钠使用和结局时间"  

while(TRUE) {
  namelist_mods = unique(mscm$stay_id)
  data_filtered = mscm_mods[(0:0),]
  for (i in namelist_mods) {
    totide = mscm_mods %>%
      filter(stay_id == i)
    index = which(totide$fail == 1)[1]
    if (!is.na(index)){
      cleaned = totide[seq_len(index),]
      data_filtered = rbind(data_filtered,cleaned)}
    else{
      data_filtered = rbind(data_filtered,totide)}
  }
  rm(totide)
  rm(cleaned)
  break
}
data_filtered = data_filtered %>%
  mutate_at(vars("min_bc":"max_potassium"),as.numeric)
data_filtered
set.seed(666)
'mice插补'
df_input = kNN(data_filtered[,c('min_bc','max_ag','min_ph','max_glucose','max_chloride','min_calcium','max_potassium','max_creat')])
newdata_sub1_mscm_sub1_mods = cbind(complete(df_input),data_filtered)

'生存分析'
mods_final = newdata_sub1_mscm_sub1_mods[,-c(19:26)]
mods_final$bc_label = ifelse(is.na(mods_final$bc_use) == T, 0, 1)
head(mods_final)
?ipwtm
temp_mods = ipwtm(exposure = bc_label,family = "survival",
                  numerator = ~sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender,
                  denominator = ~ sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender+
                    max_ag+max_glucose+min_ph+min_bc+max_chloride+min_calcium+max_potassium+max_creat,
                  id = stay_id, tstart = tstart, timevar = tstop, type = "first", data = mods_final)
mods_final_sub1 = subset(mods_final,)
temp_mods_sub1 = ipwtm(exposure = bc_label,family = "survival",
                       numerator = ~sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender,
                       denominator = ~ sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender+
                         max_ag+max_glucose+min_bc+max_chloride+max_potassium+max_creat,
                       id = stay_id, tstart = tstart, timevar = tstop, type = "first", data = mods_final_sub1)
surv_mods_mscm_sub1=survfit(Surv(mods_final$tstart,mods_final$fail)~bc_label,data = mods_final,weights = temp_mods$ipw.weights)
mods_landmark = jskm()
ggsurvplot(surv_mods_mscm_sub1,data = mods_final,conf.int = T,pval=T,
           xlim = c(0,7),
           ggtheme = theme_bw(base_size=7),
           pval.coord = c(5.6, 0.1),
           pval.size = 3,
           legend.title = '',
           legend = 'top',
           legend.labs = c("Control group", "Bicarbonate group"),
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           xlab = 'Time (days)',
           font.main = "Arial",
           ylab = 'SOFA exceeding threshold (%)',
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           palette =c("#E7B800", "#2E9FDF"),
           break.x.by = 2,
           fun = 'event',
)
rm(cox_gcs_mscm_sub1)
cox_mods_mscm_sub1 = coxph(Surv(tstart,tstop,fail)~bc_label,data = mods_final,weights = temp_mods$ipw.weights)
cox_mods_mscm_sub1 = coxph(Surv(tstart,tstop,fail)~bc_label+ max_ag+max_glucose+min_ph+min_bc+
                        max_chloride+min_calcium+max_potassium+max_creat,data = mods_final,weights = temp_mods$ipw.weights)

autoReg(cox_mods_mscm_sub1)

#************
'gcs的MSCM结局'
termin_gcs = newdata_sub1[,c('stay_id','timetogcs','egcs')]
names(termin_gcs) = c('stay_id',"tstop","end_gcs")
termin_gcs$tstop = ceiling(termin_gcs$tstop)
termin_gcs = termin_gcs %>%
  filter(end_gcs == 1)
mscm_gcs = mscm
mscm_gcs$bc_label = 0
mscm_gcs$fail = 0

for (i in 1:nrow(termin_gcs)) {
  for (j in 1:nrow(mscm_gcs)){
    if(termin_gcs$stay_id[i] == mscm_gcs$stay_id[j] & termin_gcs$tstop[i] == mscm_gcs$tstop[j])
    {mscm_gcs$fail[j] = 1}
  }
}

"上面代码合并碳酸氢钠使用和结局时间"  

while(TRUE) {
  namelist_gcs = unique(mscm$stay_id)
  data_filtered = mscm_gcs[(0:0),]
  for (i in namelist_gcs) {
    totide = mscm_gcs %>%
      filter(stay_id == i)
    index = which(totide$fail == 1)[1]
    if (!is.na(index)){
      cleaned = totide[seq_len(index),]
      data_filtered = rbind(data_filtered,cleaned)}
    else{
      data_filtered = rbind(data_filtered,totide)}
  }
  rm(totide)
  rm(cleaned)
  break
}
data_filtered = data_filtered %>%
  mutate_at(vars("min_bc":"max_potassium"),as.numeric)
data_filtered
set.seed(666)
'mice插补'
df_input = kNN(data_filtered[,c('min_bc','max_ag','min_ph','max_glucose','max_chloride','min_calcium','max_potassium','max_creat')])
newdata_sub1_mscm_sub1_gcs = cbind(complete(df_input),data_filtered)

'生存分析'
gcs_final = newdata_sub1_mscm_sub1_gcs[,-c(19:26)]
gcs_final$bc_label = ifelse(is.na(gcs_final$bc_use) == T, 0, 1)
head(gcs_final)
?ipwtm
temp_gcs = ipwtm(exposure = bc_label,family = "survival",
                 numerator = ~sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender,
                 denominator = ~ sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender+
                   max_ag+max_glucose+min_ph+min_bc+max_chloride+min_calcium+max_potassium+max_creat,
                 id = stay_id, tstart = tstart, timevar = tstop, type = "first", data = gcs_final)
gcs_final_sub1 = subset(gcs_final,)
temp_gcs_sub1 = ipwtm(exposure = bc_label,family = "survival",
                      numerator = ~sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender+
                        max_ag+max_glucose+min_ph+min_bc+max_chloride+min_calcium+max_potassium+max_creat,
                      id = stay_id, tstart = tstart, timevar = tstop, type = "first", data = gcs_final_sub1)
surv_gcs_mscm_sub1=survfit(Surv(gcs_final$tstart,gcs_final$fail)~bc_label,data = gcs_final,weights = temp_gcs$ipw.weights)
ggsurvplot(surv_gcs_mscm_sub1,data = gcs_final,title='gcs_mscm_sub1',conf.int = T,xlim = c(0,15),fun = 'event')
rm(cox_gcs_mscm_sub1)
cox_gcs_mscm_sub1 = coxph(Surv(tstart,tstop,fail)~bc_label,data = gcs_final,weights = temp_gcs$ipw.weights)
cox_gcs_mscm_sub1 = coxph(Surv(tstart,tstop,fail)~bc_label+max_ag+max_glucose+min_ph+min_bc+max_chloride+
                       min_calcium+max_potassium+max_creat,data = gcs_final,weights = temp_gcs$ipw.weights)
autoReg(cox_gcs_mscm_sub1)

#************
'sofar的MSCM结局'
termin_sofar = newdata_sub1[,c('stay_id','timetorenal','renal')]
names(termin_sofar) = c('stay_id',"tstop","end_sofar")
termin_sofar$tstop = ceiling(termin_sofar$tstop)
termin_sofar = termin_sofar %>%
  filter(end_sofar == 1)
mscm_sofar = mscm
mscm_sofar$bc_label = 0
mscm_sofar$fail = 0

for (i in 1:nrow(termin_sofar)) {
  for (j in 1:nrow(mscm_sofar)){
    if(termin_sofar$stay_id[i] == mscm_sofar$stay_id[j] & termin_sofar$tstop[i] == mscm_sofar$tstop[j])
    {mscm_sofar$fail[j] = 1}
  }
}

"上面代码合并碳酸氢钠使用和结局时间"  

while(TRUE) {
  namelist_sofar = unique(mscm$stay_id)
  data_filtered = mscm_sofar[(0:0),]
  for (i in namelist_sofar) {
    totide = mscm_sofar %>%
      filter(stay_id == i)
    index = which(totide$fail == 1)[1]
    if (!is.na(index)){
      cleaned = totide[seq_len(index),]
      data_filtered = rbind(data_filtered,cleaned)}
    else{
      data_filtered = rbind(data_filtered,totide)}
  }
  rm(totide)
  rm(cleaned)
  break
}
data_filtered = data_filtered %>%
  mutate_at(vars("min_bc":"max_potassium"),as.numeric)
data_filtered
set.seed(666)
'mice插补'
df_input = kNN(data_filtered[,c('min_bc','max_ag','min_ph','max_glucose','max_chloride','min_calcium','max_potassium','max_creat')])
newdata_sub1_mscm_sub1_sofar = cbind(complete(df_input),data_filtered)

'生存分析'
sofar_final = newdata_sub1_mscm_sub1_sofar[,-c(19:26)]
sofar_final$bc_label = ifelse(is.na(sofar_final$bc_use) == T, 0, 1)
head(sofar_final)
?ipwtm
temp_sofar = ipwtm(exposure = bc_label,family = "survival",
                   numerator = ~sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender,
                   denominator = ~ sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender+
                     max_ag+max_glucose+min_ph+min_bc+max_chloride+min_calcium+max_potassium+max_creat,
                   id = stay_id, tstart = tstart, timevar = tstop, type = "first", data = sofar_final)
sofar_final_sub1 = subset(sofar_final,)
temp_sofar_sub1 = ipwtm(exposure = bc_label,family = "survival",
                        numerator = ~sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender,
                        denominator = ~ sapsii+admission_age+vasopressin+sirs+dialysis+electivesurgery+gender+
                          max_ag+max_glucose+min_bc+max_chloride+max_potassium+max_creat,
                        id = stay_id, tstart = tstart, timevar = tstop, type = "first", data = sofar_final_sub1)
surv_sofar_mscm_sub1=survfit(Surv(sofar_final$tstart,sofar_final$fail)~bc_label,data = sofar_final,weights = temp_sofar$ipw.weights)
ggsurvplot(surv_sofar_mscm_sub1,data = sofar_final,title='sofar_mscm_sub1',conf.int = T,xlim = c(0,15),fun = 'event')
rm(cox_gcs_mscm_sub1)
cox_sofar_mscm_sub1 = coxph(Surv(tstart,tstop,fail)~bc_label,data = sofar_final,weights = temp_sofar$ipw.weights)
cox_sofar_mscm_sub1 = coxph(Surv(tstart,tstop,fail)~bc_label+max_ag+max_glucose+min_ph+min_bc+
                         max_chloride+min_calcium+max_potassium+max_creat,data = sofar_final,weights = temp_sofar$ipw.weights)

autoReg(cox_sofar_mscm_sub1)

