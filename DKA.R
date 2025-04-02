'0-包的导入'
library("MatchIt")
library("dplyr")
library("survival")
library("survminer")
library("ipw")
library("contTimeCausal")
library("tableone")
library("WeightIt")
library("mice")
library("corrplot")
library("imputeTS")
library(VIM)
library(autoReg)
library(glmnet)
library(ggplot2)
library(gridExtra)

'1.基线数据清洗'
'1.1基线数据清洗-筛选入院之后的数据'






ph_select=ph %>%
  left_join (namelist, by='stay_id') %>%
  filter(ph_time > admittime)
ag_select=ag %>%
  left_join (namelist, by='stay_id') %>%
  filter(ag_time > admittime)
glucose_select=glucose %>%
  left_join (namelist, by='stay_id') %>%
  filter(glucose_time > admittime)
sodium_select=sodium %>%
  left_join (namelist, by='stay_id') %>%
  filter(sodium_time > admittime)
potassium_select=potassium %>%
  left_join (namelist, by='stay_id') %>%
  filter(potassium_time > admittime)
chloride_select=chloride %>%
  left_join (namelist, by='stay_id') %>%
  filter(chloride_time > admittime)
bc_select=bc %>%
  left_join (namelist, by='stay_id') %>%
  filter(bc_time > admittime)
bun_select=bun %>%
  left_join (namelist, by='stay_id') %>%
  filter(bun_time > admittime)
creat_select=creatinine %>%
  left_join (namelist, by='stay_id') %>%
  filter(creat_time > admittime)
gcs_select=gcs %>%
  left_join (namelist, by='stay_id') %>%
  filter(gcs_time > admittime)
sofa_select=sofa %>%
  left_join (namelist, by='stay_id') %>%
  filter(sofa_time > admittime)
calcium_select=calcium %>%
  left_join (namelist, by='stay_id') %>%
  filter(calcium_time > admittime)
kdigo_select= kdigo %>%
  left_join (namelist, by='stay_id') %>%
  filter(kdigo_time > admittime)


'1.2基线内在时依变量的时间'
base_ag <- ag_select %>%
  group_by(stay_id) %>%
  summarise(ag_base=ag[which.min(ag_time)])  

base_ph <- ph_select %>%
  group_by(stay_id) %>%
  summarise(ph_base=ph[which.min(ph_time)])  

base_sodium <- sodium_select %>%  
  group_by(stay_id) %>%  
  summarise(sodium_base=sodium[which.min(sodium_time)])  

base_potassium <- potassium_select %>%  
  group_by(stay_id) %>%  
  summarise(potassium_base=potassium[which.min(potassium_time)])

base_chloride <- chloride_select %>%  
  group_by(stay_id) %>%  
  summarise(chloride_base=chloride[which.min(chloride_time)])

base_bc <- bc_select %>%  
  group_by(stay_id) %>%  
  summarise(bc_base=bc[which.min(bc_time)])

base_calcium <- calcium_select %>%  
  group_by(stay_id) %>%  
  summarise(calcium_base=calcium[which.min(calcium_time)])

base_bun <- bun_select %>%  
  group_by(stay_id) %>%  
  summarise(bun_base=bun[which.min(bun_time)])

base_creatinine <- creat_select %>%  
  group_by(stay_id) %>%  
  summarise(creatinine_base=creat[which.min(creat_time)])

base_kdigo <- kdigo_select %>%  
  group_by(stay_id) %>%  
  summarise(kdigo=kdigo[which.min(kdigo_time)])

base_sofa <- sofa_select %>%  
  group_by(stay_id) %>%  
  summarise(sofat=sofa_24hours[which.min(sofa_time)],sofar=renal_24hours[which.min(sofa_time)],sofac=cardiovascular_24hours[which.min(sofa_time)])

base_gcs <- gcs_select %>%  
  group_by(stay_id) %>%  
  summarise(gcs=gcs_first[which.min(gcs_time)])

'--------旧结局时间定义,这个不用了!!!!!!'
end_bc_gcs <- gcs_bc %>%  
  group_by(stay_id) %>%  
  summarise(end_gcs=min(gcs),gcs_endtime=gcs_charttime[which.min(gcs)])

end_nobc_gcs <- gcs_nobc %>%  
  group_by(stay_id) %>%  
  summarise(end_gcs=min(gcs),gcs_endtime=gcs_charttime[which.min(gcs)])

end_bc_creat<- creat_bc %>%  
  group_by(stay_id) %>%  
  summarise(end_creat=max(creat),creat_endtime=creat_charttime[which.max(creat)])

end_nobc_creat<- creat_nobc %>%  
  group_by(stay_id) %>%  
  summarise(end_creat=max(creat),creat_endtime=creat_charttime[which.max(creat)])

end_bc_sofar <- sofa_bc %>%  
  group_by(stay_id) %>%  
  summarise(end_sofar=max(renal_24hours), sofa_endtime=sofa_charttime[which.max(renal_24hours)])

end_nobc_sofar <- sofa_nobc %>%  
  group_by(stay_id) %>%  
  summarise(end_sofar=max(renal_24hours), sofa_endtime=sofa_charttime[which.max(renal_24hours)])

end_bc_sofa <- sofa_bc %>%  
  group_by(stay_id) %>%  
  summarise(end_sofa=max(sofa_24hours), mods_endtime=sofa_charttime[which.max(sofa_24hours)])

end_nobc_sofa <- sofa_nobc %>%  
  group_by(stay_id) %>%  
  summarise(end_sofa=max(sofa_24hours), mods_endtime=sofa_charttime[which.max(sofa_24hours)])

'----------旧终点时间,这个已经不用了!!!!!!!'



'1.3.1.nobc备用代码'
nobc$nobc_id=nobc[-which(total$ %in% nobc$stay_id),]
'1.3.2数据保存备用'
save_df_initial_20240201=df_initial#保存的,数据搞乱了就找这一版
table(df_initial$bc_use)
'1.3.3 二次入院,停留时间'

second_admission=unique(second_admission)
second_admission$secadmit=1
df_initial = df_initial %>%
  left_join(second_admission, by = 'subject_id')
'1.4.糖尿病分类'
group1 = namelist %>%
  left_join(group_1_diabetes,by = 'stay_id')
names(group_2_diabetes) = c('stay_id','icd_code_2')
group1and2 = group1 %>%
  left_join(group_2_diabetes,by = 'stay_id')
group1and2$type = ifelse(!is.na(group1and2$icd_code) & !is.na(group1and2$icd_code_2), 'others',
                         ifelse(!is.na(group1and2$icd_code), 1,
                                ifelse(!is.na(group1and2$icd_code_2), 2, 'others')))
diabetes_type = group1and2[,c('stay_id','type')]

'2.基线表格初步绘制成型'
'------从此处开始运行代码--------'
df_initial <- Reduce(function(x, y) merge(x, y, by = "stay_id", all = TRUE), 
                   list(patients_all,diabetes_type,charlson_baseline,dialysis_baseline,oasis_baseline,sapsii,sirs_baseline,vasopressin,vent_baseline,weight_baseline,
                        #上面是直接抠的
                        base_ag,base_ph,base_glucose,base_sodium,base_potassium,base_chloride,base_calcium,base_bc,base_bun,base_creatinine,base_gcs,base_sofa,
                        #这是电解质基线，以及sofa和gcs的基线
                        base_kdigo
                        #这是基线KDIGO
                        #死亡结局已经包含在大表当中
                        ))

修剪后重新导入数据

write_xlsx(df_exclusion,"beta_data.xlsx")
rm(end_bc_sofa,end_nobc_sofa)
#加上脑肾结局
end_creat=rbind(end_bc_creat,end_nobc_creat)
end_gcs=rbind(end_bc_gcs,end_nobc_gcs)
end_sofar=rbind(end_bc_sofar,end_nobc_sofar)
end_sofa=rbind(end_bc_sofa,end_nobc_sofa)

write_xlsx(df_merge,"beta_data.xlsx")
head(df_exclusion)
head(df_merge)
df_merge=Reduce(function(x, y) merge(x, y, by = "stay_id", all = TRUE), 
                list(df_initial,end_creat,end_gcs,end_sofar,end_sofa))
#加上糖尿病病型，兼有多种的，应该归于others or complex

#先做一个肾脏的
####我们暂且把这个叫做trymimic表
bc_select=bc_time_vary %>%
  left_join (namelist, by='stay_id') %>%
  filter(charttime > admittime)

bc_use_namelist = unique(bc_select[,c('stay_id')])

View(bc_use_namelist)
bc_use_namelist$bc_label = 1

df_merge = df_merge %>%
  left_join(bc_use_namelist,by = 'stay_id')
View(df_merge)

df_merge$bc_label = ifelse(is.na(df_merge$bc_label),0,1)
table(df_merge$bc_label)
write.csv(df_merge,file = '暑假重置版.csv')
df_merge = 暑假重置版
names(second_admission) = c('subject_id',"hospital_round","2ndtime")
df_try=select(df_merge,stay_id,admittime,bc_label,gender,type,dialysis,admission_age,
              oasis,electivesurgery,myocardial_infarct,congestive_heart_failure,chronic_pulmonary_disease,mild_liver_disease,
              severe_liver_disease,renal_disease,malignant_cancer,
              sapsii,sirs,vasopressin,vent_use,ag_base,ph_base,glucose_base,sodium_base,potassium_base,chloride_base,
              calcium_base,bc_base,bun_base,
              creatinine_base,gcs,sofat,sofar,kdigo,los_hospital,los_icu,
              end_sofa,end_sofar,end_gcs,end_creat,creat_endtime,gcs_endtime,sofa_endtime,mods_endtime)

#注意，时间变量无法插补，一旦有空白请全部删除
#这个filter函数是筛选出来的部分，不是筛选的部分
View(df_try)
df_try=df_try %>%
  filter(!is.na(end_creat) & !is.na(creat_endtime) & !is.na(end_gcs) & !is.na(gcs_endtime) & !is.na(end_sofar) & !is.na(sofa_endtime) & !is.na(end_sofa) & !is.na(mods_endtime))
numberofmissing = colSums(is.na(df_try))
write.csv(df_try,file = 'df_try_before_imput.csv')


View(df_try)
'使用序贯KNN进行填补,终点缺失的不填,只填基线缺失,mscm的同理'
cols_to_impute <- 17:34

# 选择不需要填补的列
cols_to_keep <- setdiff(1:ncol(df_try), cols_to_impute)

# 执行序贯KNN填补
imputed_data <- kNN(df_try[, cols_to_impute])

# 将填补后的数据和原始数据合并
final_data <- cbind(df_try[, cols_to_keep], imputed_data)
df_try = final_data
write.csv(df_try,file = 'df_try_after_imput.csv')


'----这一段暂时不需要------'
'钾离子补充'
#把bc_use单独存起来
bc_use_namelist = bc_use
bc_use = Reduce(function(x, y) merge(x, y, by = "stay_id", all = F), 
                list(df_try,bc_use))

nobc_use = df_try[,c('stay_id','bc_use','admittime')]
nobc_use = df_try %>%
  subset(bc_use != 1)

potassium_sup_nobc = potassium_add %>%
  left_join (nobc_use, by='stay_id') %>%
  filter(charttime > admittime)
potassium_sup_bc = potassium_add %>%
  left_join (bc_use, by='stay_id') %>%
  filter(charttime > first_bc_charttime)
potassium_sup_bc = potassium_sup_bc[,c('stay_id')]
potassium_sup_nobc = potassium_sup_nobc[,c('stay_id')]
potassium_use = rbind(potassium_sup_bc,potassium_sup_nobc)
potassium_use = potassium_use%>%
  group_by(stay_id) %>%           # 按照 ID 列分组
  mutate(n_potassium = n()) %>% 
  distinct()%>%  
  ungroup()   
'钙离子补充'
calcium_sup_nobc = calcium_add %>%
  left_join (nobc_use, by='stay_id') %>%
  filter(charttime > admittime)
calcium_sup_bc = calcium_add %>%
  left_join (bc_use, by='stay_id') %>%
  filter(charttime > first_bc_charttime)
calcium_sup_bc = calcium_sup_bc[,c('stay_id')]
calcium_sup_nobc = calcium_sup_nobc[,c('stay_id')]
calcium_use = rbind(calcium_sup_bc,calcium_sup_nobc)
calcium_use = calcium_use%>%
  group_by(stay_id) %>%           # 按照 ID 列分组
  mutate(n_calcium = n()) %>% 
  distinct()%>%  
  ungroup()   
df_try = Reduce(function(x, y) merge(x, y, by = "stay_id", all = T), 
                list(df_try,potassium_use,calcium_use))
df_try$n_potassium= ifelse(is.na(df_try$n_potassium),0,df_try$n_potassium)
df_try$n_calcium= ifelse(is.na(df_try$n_calcium),0,df_try$n_calcium)
'---这一段暂时不需要------'



'3.匹配,基线表绘制以及相关矩阵'
'3.1.1.首先转换待match的类型'
"'stay_id','admittime','bc_use','gender','type','dialysis','admission_age','
oasis','electivesurgery','myocardial_infarct','congestive_heart_failure','chronic_pulmonary_disease','mild_liver_disease','
severe_liver_disease','renal_disease','malignant_cancer','
sapsii','sirs','vasopressin','vent_use','ag_base','ph_base','glucose_base','sodium_base','potassium_base','chloride_base','bc_base','bun_base','
creatinine_base','gcs','sofat','sofar','kdigo','end_creat','creat_endtime','end_gcs','gcs_endtime','
end_sofar','sofa_endtime','end_sofa','mods_endtime','los_hospital','los_icu'"
df_initial = beta_data
df_try = df_try %>%
  mutate_at(vars("admission_age":"malignant_cancer"),as.numeric)#列的类型转化
df_try = df_try %>%
  mutate_at(vars("ag_base":"kdigo"),as.numeric)#列的类型转化
df_try = df_try %>%
  mutate_at(vars("dialysis",'gender','vent_use','electivesurgery','vasopressin','myocardial_infarct','congestive_heart_failure',
                 'chronic_pulmonary_disease','mild_liver_disease',
                 'severe_liver_disease','renal_disease','malignant_cancer'),as.factor)
'3.1.匹配'
newdata=df_try
set.seed(666)
newdata$bc_label = ifelse(is.na(newdata$bc_label), 0, 1)
matching_data=matchit(bc_label ~ sirs+sapsii+oasis+ag_base+ph_base+calcium_base+bc_base,data=newdata,method='nearest',ratio=1)
summary(matching_data)
matched_data = match.data(matching_data)
newdata = matched_data
write.csv(matched_data,'matched_data.2.23.csv')
density_plot <- ggplot(matched_data, aes(x = bun_base, fill = bc_use,group = bc_use)) +
  geom_density(alpha = 0.5) +
  labs(title = "Kernel Density Estimation Plot (PSM Matched)",
       x = 'variables',
       y = "Density",
       fill = "bc_use") +
  theme_minimal()
print(density_plot)
densityplot(creat_baseline)

#基线表
tables1 = CreateTableOne(vars= c('admission_age','gender','type','dialysis',
                                'oasis','electivesurgery','myocardial_infarct','congestive_heart_failure','chronic_pulmonary_disease',
'mild_liver_disease','severe_liver_disease','renal_disease','malignant_cancer',
'sapsii','sirs','vasopressin','vent_use','ag_base','ph_base','glucose_base','sodium_base','potassium_base',
'chloride_base','calcium_base','bc_base','bun_base','creatinine_base','gcs','sofat','sofar','kdigo','los_hospital','los_icu'),
               strata = "bc_label",data = newdata,test = T,smd = T)
table_output = print(tables1,smd = T,nonparametric = T)
write.csv(table_output,file = 'tables1.csv')
library(forestplot)
SMD_long <- tidyr::pivot_longer(SMD, -Variables, names_to = "Group", values_to = "SMD")

# 添加行号用于连接
SMD_long$row <- 1:nrow(SMD_long)
ggplot(SMD_long, aes(x = Variables, y = SMD, group = Group, color = Group)) +
  geom_point() +  # 添加散点
  geom_line(aes(group = row)) +  # 添加散点的连接线
  labs(title = "Standardized Mean Difference (SMD) before and after PSM",
       x = "Variable",
       y = "Standardized Mean Difference (SMD)") +  # 添加标题和轴标签
  theme_bw() +    # 设置主题为白色背景
  theme(legend.position = "bottom") +  # 设置图例位置
  scale_color_manual(values = c("Before" = "blue", "After" = "red"))  # 设置颜色
#作差。做出所有的差值作为其结局
#首先把它们都转化为数,然后再将其转化为是非关系
#时间标准化。将时间全部作差，按照天数标准化
typeof(newdata$end_creat)
typeof(newdata$creatinine_base)

'做一下相关矩阵'
'bc_use'
'定指标'
var_fixed = newdata[,c('type','gender','admission_age','dialysis','oasis','electivesurgery','sapsii','sirs','vasopressin',
                       'vent_use')]
var_fixed$gender = ifelse(var_fixed$gender == 'M',0,1)
var_fixed = var_fixed %>%
  mutate_at(vars("gender":"vent_use"),as.numeric)
cor_matrix_fix <- cor(var_fixed)
corrplot(cor_matrix_fix, method = "color", addCoef.col = "black",tl.col = "black",order = "FPC")
'时依协变量-动指标'
var_time = newdata[,c('ag_base','ph_base','glucose_base','calcium_base','sodium_base','potassium_base','chloride_base',
                      'bc_base','bun_base','creatinine_base','gcs','sofat','sofar','kdigo')]
names(var_time) = c('AG','PH','Glucose','Calcium','Sodium','Potassium','Chloride','Bicarbonate','Bun','Creatinine',
                 'GCS','SOFA','SOFA (renal)','KDIGO')
cor_matrix_time <- cor(var_time)

var_time = var_time %>%
  mutate_at(vars("ag_base":"kdigo"),as.numeric)
corrplot(cor_matrix_time, method = "color", addCoef.col = "black",tl.col = "black",order = "FPC")
is.na(newdata$ag_base)



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
'注:修改前的newdata在newdata_save中,在外面也有储存'
write.csv(newdata,file = 'newdata.2.23.csv')
newdata
densityplot(newdata$bc_base,xlim = c(5,30))
newdata_sub1 = subset(newdata, bc_base> 17)#做完了打个叉
newdata_sub1 = subset(newdata, bc_base<= 17)
densityplot(newdata$sofar)
newdata_sub1 = subset(newdata, creatinine_base >= 1.5)
newdata_sub1 = subset(newdata, sofar <=2)
densityplot(newdata$glucose_base,xlim = c(0,500))
newdata_sub1 = subset(newdata, glucose >2)
newdata_sub1 = subset(newdata, glucose >2)
densityplot(newdata$ag_base)

densityplot(newdata$creatinine_base)

densityplot(newdata$chloride_base,xlim = c(90,120))

densityplot(newdata$ph_base,xlim = c(7.1,7.5))

densityplot(newdata$calcium_base,xlim = c(6,11))

densityplot(newdata$sirs)


newdata = newdata_亚组前保存_2月29日


'mods'
grid.arrange(cmods$plot, ccreat$plot,crenal$plot,cgcs$plot
             )
plot_grid(cmods$plot, ccreat$plot,crenal$plot,cgcs$plot)
labels <- c("A", "B", "C", "D", "E", "F", "G", "H")
label_grobs <- lapply(labels, function(label) {
  textGrob(label, x = 0, y = 1, just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
})
grid.arrange(cmods$plot, mmods$plot,ccreat$plot,mcreat$plot,
             crenal$plot,mrenal$plot,cgcs$plot,mgcs$plot,ncol = 4)
grid.arrange()
View(newdata)
newdata$dmods=newdata$end_sofa-newdata$sofat
newdata$timetomods=difftime(newdata$mods_endtime,newdata$admittime,units = "days")
newdata$mods= ifelse(newdata$dmods >= 2, 1,0)
newdata$timetomods=ifelse(newdata$mods == 0,"",newdata$timetomods)
newdata$timetomods=as.numeric(newdata$timetomods)#列的类型转化
newdata$mods=as.numeric(newdata$mods)#列的类型转化
surv_mods=survfit(Surv(newdata$timetomods,newdata$mods)~bc_use,data = newdata)
'图片需要美化'
cmods = ggsurvplot(surv_mods,data = newdata,conf.int = F,pval=T,
           xlim = c(0,7),
           ylim = c(0,1),
           ggtheme = theme_classic(base_size=7),
           axes.offset=FALSE,
           pval.coord = c(5.1, 0.1),
           pval.size = 3,
           lwd = 0.5,
           legend.title = '',
           legend = c('top'),
           legend.labs = c("Control group", "Bicarbonate group"),
           xlab = 'Days since admission',
           font.main = "Helvetica",
           ylab = 'SOFA exceeding',
           surv.median.line = "hv", # Specify median survival
           palette = 'lancet' ,
           break.x.by = 1,
           expand=c(0,0),
           fun = 'event'
)
install.packages("gridExtra")
library(gridExtra)
grid.arrange(plot$plot, plot$plot,ncol = 2)
cox_mods=coxph(Surv(newdata$timetomods,newdata$mods)~bc_use,data = newdata)
autoReg(cox_mods)
cox_mods=coxph(Surv(newdata$timetomods,newdata$mods)~bc_use+#分组变量c("#E7B800", "#2E9FDF")
                 gender+admission_age+type+sirs+ag_base+ph_base+glucose_base+sodium_base+potassium_base+chloride_base+calcium_base+
                 bc_base+bun_base+creatinine_base+gcs+sofat+sofar+kdigo#协变量
               ,data = newdata)
#LASSO死亡
autoReg(cox_mods)




mods_check = subset(newdata, mods == 1)




'4.2.肌酐,肾脏指标1'



#肌酐上升,阈值和老师讨论一下

newdata$timetocreat=difftime(newdata$creat_endtime,newdata$admittime,units = 'days')
newdata$dcreat=newdata$end_creat-newdata$creatinine_base
newdata$creat= ifelse(newdata$dcreat >= 0.3, 1,0)
newdata$timetocreat=ifelse(newdata$creat == 0,"",newdata$timetocreat)
newdata$timetocreat=as.numeric(newdata$timetocreat)#列的类型转化
newdata$creatinine_base =as.numeric(newdata$creatinine_base)#列的类型转化
surv_creat=survfit(Surv(newdata$timetocreat,newdata$creat)~bc_label,data = newdata)
ccreat = ggsurvplot(surv_creat,data = newdata,conf.int = F,pval=T,
                    xlim = c(0,7),
                    ylim = c(0,1),
                    ggtheme = theme_classic(base_size=7),
                    axes.offset=FALSE,
                    pval.coord = c(5.1, 0.1),
                    pval.size = 3,
                    lwd = 0.5,
                    legend.title = '',
                    legend = c('top'),
                    legend.labs = c("Control group", "Bicarbonate group"),
                    xlab = 'Days since admission',
                    font.main = "Helvetica",
                    ylab = 'Cr exceeding',
                    surv.median.line = "hv", # Specify median survival
                    palette = 'lancet' ,
                    break.x.by = 1,
                    expand=c(0,0),
                    fun = 'event'
)
plot(surv_creat)
cox_creat=coxph(Surv(newdata$timetocreat,newdata$creat)~bc_use+gender+admission_age+dialysis+oasis+electivesurgery+sapsii+sirs+vasopressin+
                  vent_use+ag_base+glucose_base+sodium_base+potassium_base+chloride_base+
                  bc_base+bun_base+creatinine_base+gcs+sofat+kdigo,data = newdata)
summary(cox_creat)
creat_check = subset(newdata, creat == 1)
cox_creat=coxph(Surv(newdata$timetocreat,newdata$creat)~bc_use,data = newdata)
autoReg(cox_creat)
ggsurvplot(surv_creat,data = newdata_sub1,title='creat',conf.int=T,pval=T,xlim=c(0,10000))
'4.3肾脏SOFA,这是指标2'
newdata$drenal=newdata$end_sofar-newdata$sofar
newdata$timetorenal=difftime(newdata$sofa_endtime,newdata$admittime,units = 'days')
newdata$renal= ifelse(newdata$drenal >= 1, 1,0)
newdata$timetorenal=ifelse(newdata$renal == 0,"",newdata$timetorenal)
newdata$timetorenal=as.numeric(newdata$timetorenal)#列的类型转化
newdata$renal=as.numeric(newdata$renal)#列的类型转化
surv_renal=survfit(Surv(newdata$timetorenal,newdata$renal)~bc_label,data = newdata)
crenal=ggsurvplot(surv_renal,data = newdata,conf.int = F,pval=T,
             xlim = c(0,7),
             ylim = c(0,1),
             ggtheme = theme_classic(base_size=7),
             axes.offset=FALSE,
             pval.coord = c(5.1, 0.1),
             pval.size = 3,
             lwd = 0.5,
             legend.title = '',
             legend = c('top'),
             legend.labs = c("Control group", "Bicarbonate group"),
             xlab = 'Days since admission',
             font.main = "Helvetica",
             ylab = 'SOFA (renal) exceeding',
             surv.median.line = "hv", # Specify median survival
             palette = 'lancet' ,
             break.x.by = 1,
             expand=c(0,0),
             fun = 'event'
)

cox_renal=coxph(Surv(newdata$timetorenal,newdata$renal)~bc_use+ag_base+glucose_base,data = newdata,weights = newdata$distance)
cox_renal=coxph(Surv(newdata$timetorenal,newdata$renal)~bc_use,data = newdata)
autoReg(cox_renal)


'4.4.神经功能指标'
newdata$dgcs=newdata$gcs-newdata$end_gcs
newdata$timetogcs=difftime(newdata$gcs_endtime,newdata$admittime,units = 'days')
newdata$egcs= ifelse(newdata$dgcs >= 1, 1,0)
newdata$timetogcs=ifelse(newdata$egcs == 0,ifelse(newdata$los_hospital <=3,newdata$los_hospital,3),newdata$timetogcs)
newdata$egcs  = ifelse(newdata$timetogcs > 3,0,newdata$egcs)
newdata$timetogcs=ifelse(newdata$timetogcs > 3, 3, newdata$timetogcs)
newdata$timetogcs=as.numeric(newdata$timetogcs)#列的类型转化
newdata$egcs=as.numeric(newdata$egcs)#列的类型转化
gcs_check = subset(newdata, egcs == 1)
surv_gcs=survfit(Surv(newdata$timetogcs,newdata$egcs)~bc_label,data = newdata)
cox_gcs=coxph(Surv(newdata$timetogcs,newdata$egcs)~bc_label,data = newdata)
summary(cox_gcs)
cgcs=ggsurvplot(surv_gcs,data = newdata,conf.int = F,pval=T,
                xlim = c(0,3),
                ylim = c(0,1),
                ggtheme = theme_classic(base_size=7),
                axes.offset=FALSE,
                pval.coord = c(5.1, 0.1),
                pval.size = 3,
                lwd = 0.5,
                legend.title = '',
                legend = c('top'),
                legend.labs = c("Control group", "Bicarbonate group"),
                xlab = 'Days since admission',
                font.main = "Helvetica",
                ylab = 'GCS worsening',
                surv.median.line = "hv", # Specify median survival
                palette = 'lancet' ,
                break.x.by = 1,
                expand=c(0,0),
                fun = 'event'
)
density_plot <- density(newdata$glucose_base)
plot(density_plot,xlim = c(0,600))
newdata_sub1 = subset(newdata,glucose_base > 3)
surv_gcs=survfit(Surv(newdata_sub1$timetogcs,newdata_sub1$egcs)~bc_use,data = newdata_sub1)
ggsurvplot(surv_gcs,data = newdata_sub1,title='gcs',conf.int=T,xlim=c(0,10000),pval=T)
cox_gcs=coxph(Surv(newdata$timetogcs,newdata$egcs)~bc_use,data = newdata)

autoReg(cox_gcs)
#设置其哑变量的名称,等级变量不作为哑变量,作为常规变量纳入
#其他变量全部设置为哑变量
#这个过程是自动的,应该指出
colnames(newdata)

'这是所有变量的名单'
lasso_x<-model.matrix(~gender+type+dialysis+electivesurgery+myocardial_infarct+
                         congestive_heart_failure+chronic_pulmonary_disease+mild_liver_disease+
                         severe_liver_disease+renal_disease+malignant_cancer+vasopressin+vent_use+#以上是哑变量
                         admission_age+oasis+sirs+sapsii+ag_base+ph_base+
                         glucose_base+sodium_base+potassium_base+chloride_base+
                         calcium_base+bc_base+bun_base+creatinine_base+gcs+sofat+kdigo,data = newdata)
'请加上随访时间'

newdata$timetorenal = ifelse(is.na(newdata$timetorenal),ifelse(newdata$los_hospital <=7,newdata$los_hospital,7),newdata$timetorenal)
newdata$timetocreat = ifelse(is.na(newdata$timetocreat),ifelse(newdata$los_hospital <=7,newdata$los_hospital,7),newdata$timetocreat)
newdata$timetomods = ifelse(is.na(newdata$timetomods),ifelse(newdata$los_hospital <=7,newdata$los_hospital,7),newdata$timetomods)
newdata$timetodeath = ifelse(is.na(newdata$timetodeath),ifelse(newdata$los_hospital <=7,newdata$los_hospital,7),newdata$timetodeath)
x<-as.matrix(lasso_x)#指定自变量为矩阵,自变量就公用一套即可.
gcs_obj = Surv(newdata$timetogcs,newdata$egcs) #指定因变量为矩阵
cv.lasso_gcs <- cv.glmnet(x, gcs_obj, alpha = 1,family = 'cox')
best_lambda_gcs <- cv.lasso_gcs$lambda.min
lasso_gcs_model <- glmnet(x, gcs_obj, alpha = 1,family = 'cox',lambda = best_lambda_gcs)
lasso_gcs_coefficients <- coef(lasso_gcs_model)
lasso_gcs_coefficients
creat_obj = Surv(newdata$timetocreat,newdata$creat) #指定因变量为矩阵
cv.lasso_creat <- cv.glmnet(x, creat_obj, alpha = 1,family = 'cox')
best_lambda_creat <- cv.lasso_creat$lambda.min
lasso_creat_model <- glmnet(x, creat_obj, alpha = 1,family = 'cox',lambda = best_lambda_creat)
lasso_creat_coefficients <- coef(lasso_creat_model)
lasso_creat_coefficients
mods_obj = Surv(newdata$timetomods,newdata$emods) #指定因变量为矩阵
cv.lasso_mods <- cv.glmnet(x, mods_obj, alpha = 1,family = 'cox')
best_lambda_mods <- cv.lasso_mods$lambda.min
lasso_mods_model <- glmnet(x, mods_obj, alpha = 1,family = 'cox',lambda = best_lambda_mods)
lasso_mods_coefficients <- coef(lasso_mods_model)
lasso_mods_coefficients
#KDIGO评分

'4.4.二次入院'
#二次入院,做logistic回归
sec = newdata[,c('stay_id','admittime')]
second_admission = second_admission %>%
  distinct()
newdata_sec = newdata %>%
  left_join(second_admission,by = 'stay_id')
newdata_sec$first_icu_stay = ifelse(newdata_sec$first_icu_stay == 'f', 1, 0)
newdata_sec$first_icu_stay = ifelse(is.na(newdata_sec$first_icu_stay),0,1)
sec_model = glm(first_icu_stay~bc_use,data = newdata_sec)
autoReg(sec_model)

newdata_sec$secadmit = ifelse(is.na(newdata$secadmit),0,1)
second = glm(secadmit ~ bc_use,data = newdata)
summary(second)
'4.5.死亡'
#死亡,死亡没有统计学差异,因此可以放在附件当中
mortality = df_merge[,c('stay_id','hospital_expire_flag','dod')]
newdata_death = newdata %>%
  inner_join(mortality,by = 'stay_id')
newdata$death=newdata_death$hospital_expire_flag
newdata$timetodeath=newdata_death$dod - newdata_death$admittime
newdata$timetodeath=difftime(newdata_death$dod,newdata_death$admittime,units = 'days')
newdata$timetodeath=ifelse(newdata$death == 0,"",newdata$timetodeath)
newdata$timetodeath=as.numeric(newdata$timetodeath)#列的类型转化
newdata$death=as.numeric(newdata$death)#列的类型转化
surv_death=survfit(Surv(newdata$timetodeath,newdata$death)~bc_use,data = newdata)
cox_death=coxph(Surv(newdata$timetodeath,newdata$death)~bc_use,data = newdata)
autoReg(cox_death)
death_model = glm(death~bc_use,data = newdata)
autoReg(death_model)
ggsurvplot(surv_death,data = newdata,conf.int = F,pval=T,
           xlim = c(0,7),
           ylim = c(0,0.25),
           ggtheme = theme_classic(base_size=7),
           axes.offset=FALSE,
           pval.coord = c(5.1, 0.1),
           pval.size = 3,
           lwd = 0.5,
           legend.title = '',
           legend = c('top'),
           legend.labs = c("Control group", "Bicarbonate group"),
           xlab = 'Days since admission',
           font.main = "Helvetica",
           ylab = 'mortality',
           surv.median.line = "hv", # Specify median survival
           palette = 'lancet' ,
           break.x.by = 1,
           expand=c(0,0),
           fun = 'event'
)
'4.6钙盐钾盐补液次数分析'
potassium_model = glm(n_potassium~bc_use,data = newdata)
autoReg(potassium_model)
calcium_model = glm(n_calcium~bc_use,data = newdata)
autoReg(calcium_model)
autoReg(cox_mods)


'4.7 WIN RATIO'









