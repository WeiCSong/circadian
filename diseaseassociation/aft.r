library(data.table)
library(mgcv)
library(lubridate)
library(dplyr)
library(survival)
library(ggplot2)
library(CMAverse)
#library(ggrepel)
time=fread("bloodtime_participant.csv",data.table=F)
pred=fread("bart/allpred.csv",data.table=F)
pred=data.frame(pred[,1],pred[,1],pred[,3]-pred[,4],pred[,6])
colnames(pred)=c("ID","IID","shift","disturbance")
x=pred[,4]
x=qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
pred[,4]=x
x=pred[,3]
x=qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
pred[,3]=x
file=system("ls ICD/*.csv",intern=T)
covar=fread("democovar.csv",data.table=F)
bcovar=fread("lifestylecovar.csv",data.table=F)
end_date <- ymd("2025-01-01")
blood_collection <- ymd_hms(time[,2]) 
lp=match(time[,1],pred$ID)
shift=pred[lp,3]


chro=fread("/lustre/home/acct-bioxsyy/bioxsyy-user2/ukb/file/common.regenie.pheno",data.table=F)
chro=chro[match(time[,1],chro[,1]),c(1,7)]

res=fread("allbart2dis_zph.csv",data.table=F)
res=res[which(res[,4]<0.05),1]
resname=unique(res)

for (path in file){
  diag=fread(path,data.table=F)
  for(dis in colnames(diag)[-1]){   
    if(! dis %in% resname){next} 
    diagnosis_date <- ymd(diag[match(time[,1],diag[,1]),dis])  
    status <- ifelse(!is.na(diagnosis_date) & diagnosis_date > blood_collection, 1, 0)
    dtime <- ifelse(status == 1,
                 as.numeric(difftime(diagnosis_date, blood_collection, units = "days")),
                 as.numeric(difftime(end_date, blood_collection, units = "days")))
    formerdiag=which(!is.na(diagnosis_date) & diagnosis_date <= blood_collection)
    dtime[formerdiag] <- NA
    status[formerdiag] <- NA
    
    df=data.frame(status=status,time=dtime,
                  covar[match(time[,1],covar[,1]),2:16],
                  bcovar[match(time[,1],bcovar[,1]),-1],
                  shift=pred[lp,"shift"],
                  disturbance=pred[lp,"disturbance"],
                  chronotype=chro[,2])
    df=df[-which(is.na(df[,1]) | is.na(df[,2]) | is.na(df$shift) | is.na(df$chronotype)),]
    df[is.na(df)]=0
    aft_model <- survreg(Surv(time, status) ~ ., data = df,dist = "weibull")
    coef=summary(aft_model)$table
    int1=coef["disturbance",]
        
    res1=data.frame(dis=dis,metric="disturbance",int1[1],int1[2],int1[3],int1[4])
    fwrite(res1,"allbart2dis_aft.csv",append=T)
  }
}

pred=fread("bart/g1pred.csv",data.table=F)
pred=data.frame(pred[,1],pred[,1],pred[,3]-pred[,4],pred[,6])
colnames(pred)=c("ID","IID","shift","disturbance")
x=pred[,4]
x=qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
pred[,4]=x
x=pred[,3]
x=qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
pred[,3]=x

lp=match(time[,1],pred$ID)
shift=pred[lp,3]

res=fread("g1bart2dis_zph.csv",data.table=F)
res=res[which(res[,4]<0.05),1]
resname=unique(res)

for (path in file){
  diag=fread(path,data.table=F)
  for(dis in colnames(diag)[-1]){   
    if(! dis %in% resname){next} 
    diagnosis_date <- ymd(diag[match(time[,1],diag[,1]),dis])  
    status <- ifelse(!is.na(diagnosis_date) & diagnosis_date > blood_collection, 1, 0)
    dtime <- ifelse(status == 1,
                 as.numeric(difftime(diagnosis_date, blood_collection, units = "days")),
                 as.numeric(difftime(end_date, blood_collection, units = "days")))
    formerdiag=which(!is.na(diagnosis_date) & diagnosis_date <= blood_collection)
    dtime[formerdiag] <- NA
    status[formerdiag] <- NA
    
    df=data.frame(status=status,time=dtime,
                  covar[match(time[,1],covar[,1]),2:16],
                  bcovar[match(time[,1],bcovar[,1]),-1],
                  shift=pred[lp,"shift"],
                  disturbance=pred[lp,"disturbance"],
                  chronotype=chro[,2])
    df=df[-which(is.na(df[,1]) | is.na(df[,2]) | is.na(df$shift) | is.na(df$chronotype)),]
    df[is.na(df)]=0
    aft_model <- survreg(Surv(time, status) ~ ., data = df,dist = "weibull")
    coef=summary(aft_model)$table
    int1=coef["disturbance",]
        
    res1=data.frame(dis=dis,metric="disturbance",int1[1],int1[2],int1[3],int1[4])
    fwrite(res1,"g1bart2dis_aft.csv",append=T)
  }
}

pred=fread("bart/g2pred.csv",data.table=F)
pred=data.frame(pred[,1],pred[,1],pred[,3]-pred[,4],pred[,6])
colnames(pred)=c("ID","IID","shift","disturbance")
x=pred[,4]
x=qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
pred[,4]=x
x=pred[,3]
x=qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
pred[,3]=x

lp=match(time[,1],pred$ID)
shift=pred[lp,3]

res=fread("g2bart2dis_zph.csv",data.table=F)
res=res[which(res[,4]<0.05),1]
resname=unique(res)

for (path in file){
  diag=fread(path,data.table=F)
  for(dis in colnames(diag)[-1]){   
    if(! dis %in% resname){next} 
    diagnosis_date <- ymd(diag[match(time[,1],diag[,1]),dis])  
    status <- ifelse(!is.na(diagnosis_date) & diagnosis_date > blood_collection, 1, 0)
    dtime <- ifelse(status == 1,
                 as.numeric(difftime(diagnosis_date, blood_collection, units = "days")),
                 as.numeric(difftime(end_date, blood_collection, units = "days")))
    formerdiag=which(!is.na(diagnosis_date) & diagnosis_date <= blood_collection)
    dtime[formerdiag] <- NA
    status[formerdiag] <- NA
    
    df=data.frame(status=status,time=dtime,
                  covar[match(time[,1],covar[,1]),2:16],
                  bcovar[match(time[,1],bcovar[,1]),-1],
                  shift=pred[lp,"shift"],
                  disturbance=pred[lp,"disturbance"],
                  chronotype=chro[,2])
    df=df[-which(is.na(df[,1]) | is.na(df[,2]) | is.na(df$shift) | is.na(df$chronotype)),]
    df[is.na(df)]=0
    aft_model <- survreg(Surv(time, status) ~ ., data = df,dist = "weibull")
    coef=summary(aft_model)$table
    int1=coef["disturbance",]
        
    res1=data.frame(dis=dis,metric="disturbance",int1[1],int1[2],int1[3],int1[4])
    fwrite(res1,"g2bart2dis_aft.csv",append=T)
  }
}



