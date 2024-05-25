library("devtools")
library(installr)
install.Rtools()
library(`HapEstXXR`)
options(digits=2)
options(scipen=999)
is.odd <- function(x) x %% 2 != 0

Diverstiypanel=read.csv("combinedLSmean.csv",header=T)

Geno_Diverstiypanel=read.csv("diversity panel geno.csv",header=T,na.strings="")

convertSNP=function(SNP){
  majorallele=tail(names(sort(table(SNP[!SNP=="H"&!is.na(SNP)]))),2)[2] 
  minorallele=tail(names(sort(table(SNP[!SNP=="H"&!is.na(SNP)]))),2)[1]
  new_SNP=gsub(majorallele,0,SNP)
  new_SNP=gsub(minorallele,1,new_SNP)
  new_SNP=gsub("H",2,new_SNP)
  return(as.numeric(new_SNP))
}
new_geno_Diverstiypanel=as.data.frame(cbind(as.character(Geno_Diverstiypanel$ID),apply(Geno_Diverstiypanel[,-1],2,convertSNP)))

pheno=Diverstiypanel[match(new_geno_Diverstiypanel$V1, Diverstiypanel$X),]
rownames(pheno)=NULL
#colnames(pheno)[3:ncol(pheno)]=c("GW_1k","GM2","gy","SPI","HI","FE","TL","IN2wt","iNthreewt","INL2","INL3" ,"IN2PI",     
             #                    "IN3PI","IN2IN3PI","AverageNDVi")
Ttestgeno=new_geno_Diverstiypanel
Ttestgeno[Ttestgeno==2]=NA

t.test(pheno[,7]~Ttestgeno[,17])

selectedgeno=c("S4A_609462220",
"S7A_106952238",
"S2B_644162831",
"S3B_598200325",
"S5A_438411966",
"S4B_13299920",
"S4B_12533949",
"S6B_157236831",
"S2D_76639743",
"S2D_72780208",
"S2D_68734015",
"S2D_69502623",
"S2D_76639743",
"S5B_479718074")

selecedlines=cbind(pheno,Ttestgeno[,colnames(Ttestgeno)%in%selectedgeno])



write.csv(selecedlines,"selectedlines.csv")
selectedlines=read.csv(file.choose(),header=T)
sum(selectedlines[selectedlines[,3]==0,2],na.rm=T)/length(which(selectedlines[,3]==0))
sum(selectedlines[selectedlines[,3]==1,2],na.rm=T)/length(which(selectedlines[,3]==1))
selectedlines[which(selectedlines[,3]==1),3]="+"
selectedlines[which(selectedlines[,3]==0),3]="-"

whichgeno=29
whichpheno=27
sum(selectedlines[selectedlines[,whichgeno]==1,whichpheno],na.rm=T)/length(which(selectedlines[,whichgeno]==1))
sum(selectedlines[selectedlines[,whichgeno]==0,whichpheno],na.rm=T)/length(which(selectedlines[,whichgeno]==0))
selectedlines[which(selectedlines[,whichgeno]==1),whichgeno]="-"
selectedlines[which(selectedlines[,whichgeno]==0),whichgeno]="+"
write.csv(selectedlines,"selectedlines_organized.csv")


Pvaltable=matrix(NA,nrow=ncol(pheno[,2:ncol(pheno)]),ncol=ncol(Ttestgeno[,2:ncol(Ttestgeno)]))
meantable=matrix(NA,nrow=ncol(pheno[,2:ncol(pheno)]),ncol=2*ncol(Ttestgeno[,2:ncol(Ttestgeno)]))
for (i in 2:ncol(pheno)){
  for (j in 2:ncol(Ttestgeno)){
    temp=t.test(pheno[,i]~Ttestgeno[,j])
    Pvaltable[i-1,j-1]=temp$p.value
    meantable[i-1,(2*(j-1)-1)]=temp$estimate[[1]]
    meantable[i-1,2*(j-1)]=temp$estimate[[2]]
  }
}

rownames(Pvaltable)=colnames(pheno[,2:ncol(pheno)])
rownames(meantable)=colnames(pheno[,2:ncol(pheno)])
colnames(Pvaltable)=colnames(Ttestgeno[,2:ncol(Ttestgeno)])
colnames(meantable)[is.odd(1:ncol(meantable))]=colnames(Ttestgeno[,2:ncol(Ttestgeno)])

Pvaltable[Pvaltable>0.05]=NA

diverstiy_Pval=Pvaltable
diverstiy_mean=meantable
write.csv(diverstiy_Pval,"diverstiy_Pval.csv")
write.csv(diverstiy_mean,"diverstiy_mean.csv")



Pvaltable_organized=matrix(NA,nrow=0,ncol=3)
for (k in 1:nrow(Pvaltable)){
  index=!is.na(Pvaltable[k,])
  table_temp=cbind(rownames(Pvaltable)[k],cbind(colnames(Pvaltable)[index],
                                                as.vector(Pvaltable[k,index])))
         Pvaltable_organized=rbind(Pvaltable_organized,table_temp)          
 
}
write.csv(Pvaltable_organized,"Pvaltable_organized.csv")


#####################################################Stepwise mm###########################
library(MASS)
library(leaps)
library(caret)

trait=colnames(pheno)[15]
# Fit the full model 
stepgeno=Ttestgeno[,-1]
stepgeno=stepgeno[,is.na(Pvaltable[rownames(Pvaltable)==trait,])==F]
stepdata=cbind(pheno[,colnames(pheno)==trait],stepgeno)
colnames(stepdata)[1]=trait
stepdata[stepdata==2]=NA
stepdata<-stepdata[complete.cases(stepdata),]
f <- as.formula(
  paste(trait,
        paste(colnames(stepdata[,2:ncol(stepdata)]), collapse = " + "),
        sep = " ~ "))
full.model <- lm(f, data =stepdata)
# Stepwise regression model
step.model <- stepAIC(full.model, direction = "both", 
                      trace = FALSE)
summary(step.model)
#####################################################ML mm###########################
set.seed(123)
# Set up repeated k-fold cross-validation
train.control <- trainControl(method = "cv", number = 10)
# Train the model
step.model <- train(f, data = stepdata,method = "leapBackward",
                    tuneGrid = data.frame(nvmax = 1:2),
                    trControl = train.control)
step.model$results

step.model$bestTune
summary(step.model$finalModel)
coef(step.model$finalModel, 4)

