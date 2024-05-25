library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
Diverstiypanel=read.csv("combinedLSmean.csv",header=T)
linelist=read.csv("master.csv")
spring_lines=linelist[linelist$source%in%c("University of Idaho","MN"),]
winter_lines=linelist[!linelist$source%in%c("University of Idaho","MN"),]
mypheno=Diverstiypanel[Diverstiypanel$entry%in%winter_lines$ID,]

Geno_Diverstiypanel=read.csv("diversity panel geno.csv",header=T,na.strings="")
mygeno=Geno_Diverstiypanel[Geno_Diverstiypanel$ID%in%winter_lines$ID,]


mypheno=mypheno[mypheno$entry%in%mygeno$ID,]
mygeno=mygeno[mygeno$ID%in%mypheno$entry,]

mygeno[mygeno=="H"]=NA
mypheno=mypheno[match(mygeno$ID,mypheno$entry),]

rownames(mygeno)=NULL
rownames(mypheno)=NULL
colnames(mypheno)[2]="DTH"
##########################################################################
plot_data <- left_join(mypheno, mygeno, by = join_by(entry == ID))
plot_data_long <- pivot_longer(plot_data, cols = grep('_', colnames(plot_data), value = T), names_to = 'SNPname', values_to = 'GT')
plot_data_long2 <- pivot_longer(plot_data_long, cols = grep('_|GT|SNPname|entry', colnames(plot_data_long), value = T, invert = T),
  names_to = 'trait', values_to = 'BLUE')





plot_data_long_f <- plot_data_long %>% 
  filter(SNPname = 'S5A_438411966' & trait %in% c('DTH', 'GW_1k', 'FE', 'INL2', 'IN2PI'),
         SNPname = 'S4B_13299920' & trait %in% c('INL3'),
         SNPname = 'S7A_681517593' & trait %in% c('INL2wt', 'iNthreewt', 'FE', 'INL2', 'IN2PI'),
         SNPname = 'S5A_438411966' & trait %in% c('DTH', 'GW_1k', 'FE', 'INL2', 'IN2PI'),
         SNPname = 'S5A_438411966' & trait %in% c('DTH', 'GW_1k', 'FE', 'INL2', 'IN2PI'))



ggplot(data=subset(plot_data_long,!is.na(geno)),aes(x=value, y=pheno,color=DTH))+
  geom_boxplot(fill=c("orange","cornflowerblue"),outlier.size = -1,lwd=3)+
  theme(axis.text=element_text(size=40,face="bold"),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 3))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  ggtitle(mytitle)+
  theme(plot.title = element_text(color="blue", size=50, face="bold",hjust=0.5))












mytrait=mypheno$InternodebothPI
mytitle=c("IN2_3_PI")
mymarker="S5B_479718074"
plot_data=data.frame(geno=mygeno[,colnames(mygeno)==mymarker],pheno=mytrait,DTH=mypheno$DTH)
str(plot_data)

png(file=paste(mytitle,"_",mymarker,".png",sep=""),width=4,height=10,units="in",res=300)
ggplot(data=subset(plot_data,!is.na(geno)),aes(x=geno,y=pheno,color=DTH))+
  geom_boxplot(fill=c("orange","cornflowerblue"),outlier.size = -1,lwd=3)+
  theme(axis.text=element_text(size=40,face="bold"),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 3))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  ggtitle(mytitle)+
  theme(plot.title = element_text(color="blue", size=50, face="bold",hjust=0.5))
dev.off()

##########################################################################
mytrait=mypheno$yield
mytitle=c("GY")
mymarker1="S1B_636796017"
mymarker2="S3A_727396932"
plot_data=data.frame(geno=mygeno[,colnames(mygeno)%in%c(mymarker1,mymarker2)],pheno=mytrait,DTH=mypheno$DTH)
str(plot_data)

plot_data=plot_data[!(is.na(plot_data$geno.S1B_636796017)|is.na(plot_data$geno.S3A_727396932)),]
rownames(plot_data)=NULL
plot_data$new_geno=paste(plot_data$geno.S1B_636796017,plot_data$geno.S3A_727396932,sep="+")

k=compare_means(pheno ~ new_geno,  data = plot_data)
mygroup=cbind(c(k$group1[1]),c(k$group2)[1])
my_comparisons <-split(mygroup, seq(nrow(mygroup)))

png(file=paste(mytitle,"_",mymarker1,"_",mymarker2,".png",sep=""),width=8,height=10,units="in",res=300)
ggplot(data=subset(plot_data,!is.na(new_geno)),aes(x=new_geno,y=pheno))+
  geom_boxplot(fill=c("blueviolet","black","red","chartreuse2"),outlier.size = -1,lwd=3)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",hide.ns=T,vjust=0.6,size=25,bracket.size = 2.5)+
  theme(axis.text=element_text(size=40,face="bold"),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 3))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  ggtitle(mytitle)+
  theme(plot.title = element_text(color="blue", size=50, face="bold",hjust=0.5))
dev.off()

##########################################################################
mytrait=mypheno$HI
mytitle=c("HI")
mymarker1="S4A_609462220"
mymarker2="S7A_106952238"
plot_data=data.frame(geno=mygeno[,colnames(mygeno)%in%c(mymarker1,mymarker2)],pheno=mytrait,DTH=mypheno$DTH)
str(plot_data)

plot_data=plot_data[!(is.na(plot_data$geno.S4A_609462220)|is.na(plot_data$geno.S7A_106952238)),]
rownames(plot_data)=NULL
plot_data$new_geno=paste(plot_data$geno.S4A_609462220,plot_data$geno.S7A_106952238,sep="+")

k=compare_means(pheno ~ new_geno,  data = plot_data)
mygroup=cbind(c(k$group1[c(2,3,4,6)]),c(k$group2[c(2,3,4,6)]))
my_comparisons <-split(mygroup, seq(nrow(mygroup)))

png(file=paste(mytitle,"_",mymarker1,"_",mymarker2,".png",sep=""),width=8,height=10,units="in",res=300)
ggplot(data=subset(plot_data,!is.na(new_geno)),aes(x=new_geno,y=pheno))+
  geom_boxplot(fill=c("blueviolet","black","red","chartreuse2"),outlier.size = -1,lwd=3)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",hide.ns=T,vjust=0.6,size=25,bracket.size = 2.5)+
  theme(axis.text=element_text(size=40,face="bold"),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 3))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  ggtitle(mytitle)+
  theme(plot.title = element_text(color="blue", size=50, face="bold",hjust=0.5))
dev.off()

##########################################################################
mytrait=mypheno$InternodebothPI
mytitle=c("IN2_3_PI")
mymarker1="S2D_63397906"
mymarker2="S2D_63397962"
mymarker3="S5B_479718074"
plot_data=data.frame(geno=mygeno[,colnames(mygeno)%in%c(mymarker1,mymarker2,mymarker3)],pheno=mytrait,DTH=mypheno$DTH)
str(plot_data)

plot_data=plot_data[!(is.na(plot_data$geno.S2D_63397906)|is.na(plot_data$geno.S2D_63397962)|is.na(plot_data$geno.S5B_479718074)),]
rownames(plot_data)=NULL

plot_data$new_geno=paste(plot_data$geno.S2D_63397906,plot_data$geno.S2D_63397962,plot_data$geno.S5B_479718074,sep="+")
plot_data$new_2geno1=paste(plot_data$geno.S2D_63397906,plot_data$geno.S2D_63397962,sep="+")
plot_data$new_2geno2=paste(plot_data$geno.S2D_63397906,plot_data$geno.S5B_479718074,sep="+")
plot_data$new_2geno3=paste(plot_data$geno.S2D_63397962,plot_data$geno.S5B_479718074,sep="+")


k=compare_means(pheno ~ new_geno,  data = plot_data)
mygroup=cbind(c(k$group1[c(3)]),c(k$group2[c(3)]))
my_comparisons <-split(mygroup, seq(nrow(mygroup)))

png(file=paste(mytitle,"_",mymarker1,"_",mymarker2,"_",mymarker3,".png",sep=""),width=8,height=10,units="in",res=300)
ggplot(data=subset(plot_data,!is.na(new_geno)),aes(x=new_geno,y=pheno))+
  geom_boxplot(fill=c("blueviolet","red","chartreuse2"),outlier.size = -1,lwd=3)+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",hide.ns=T,vjust=0.6,size=25,
                    label.y=c(0.255), bracket.size = 2.5)+
  theme(axis.text=element_text(size=40,face="bold"),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 3))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  ggtitle(mytitle)+
  theme(plot.title = element_text(color="blue", size=50, face="bold",hjust=0.5))
dev.off()


k=compare_means(pheno ~ new_2geno2,  data = plot_data)
mygroup=cbind(c(k$group1[c(3)]),c(k$group2[c(3)]))
my_comparisons <-split(mygroup, seq(nrow(mygroup)))

png(file=paste(mytitle,"_",mymarker1,"_",mymarker3,".png",sep=""),width=8,height=10,units="in",res=300)
ggplot(data=subset(plot_data,!is.na(new_2geno2)),aes(x=new_2geno2,y=pheno))+
  geom_boxplot(fill=c("blueviolet","red","chartreuse2"),outlier.size = -1,lwd=3)+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",hide.ns=T,vjust=0.6,size=25,
                     label.y=c(0.255), bracket.size = 2.5)+
  theme(axis.text=element_text(size=40,face="bold"),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 3))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  ggtitle(mytitle)+
  theme(plot.title = element_text(color="blue", size=50, face="bold",hjust=0.5))
dev.off()


