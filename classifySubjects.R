#########################################
# classifySubjects.R - Classify subjects
# Author: Hoai Tuong Nguyen
# Created: 30/05/2013
# Modified: 20/06/2013
# CMD: R --no-save --no-restore --slave -f classifySubjects.R "--args input.dir='?' output.dir='?'"
#########################################

#Install on cluster
#install.packages("",repos="http://bioconductor.org/packages/2.12/bioc",lib="/home/nguyen/R") 

#First read and parse the arguments for on the fly run
args=(commandArgs(TRUE))

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}


library(mc)

#default working path for inner run
#setwd("MetaCARDIS/")
#setwd("C:/Users/TUONG/Documents/metacardis/7.scripts/")

#default input and output paths
if (!exists("input.dir"))
  input.dir="../1.data/metacardis"
if (!exists("output.dir"))
  output.dir="../3.results/CLASSIFICATION/Subjects/2013-06-17"

dir.create(output.dir, showWarnings = FALSE)


#Load imputed data
#seqcount.data<-get(load(sprintf("%s/DonneeSequences_t0_global_counting_table.csv.gz.opt.RData",input.dir)))
clinic.data<-get(load(sprintf("%s/T0ClinicalParametersMicrObesja_All.RData",input.dir)))
ge.data.t0<-get(load(sprintf("%s/GE_t0.RData",input.dir)))
#ge.data.t0<-read.table(sprintf("%s/GE_t0.csv",output.dir),header=T,sep=",")
go.data<-get(load(sprintf("%s/DiogenesFunctionsPerProbeMatrix.rdata",input.dir)))

names(which(go.data[,"GO:0070904"]==1))
names(which(go.data[8134,]==1))

#Pre-processing data
if(F){
  clinic.data.1<-get(load(sprintf("%s/DonneeCliniques0_6_12Sem_3avril2010.xls.imputed.RData",input.dir)))
  tmp<-clinic.data[sapply(1:49, function(x) which(clinic.data$Names == clinic.data.1$NOM[x])),]
  heigh_clin_m<-clinic.data.1$Taille[1:49]
  names(heigh_clin_m)
  clinic.data<-cbind(tmp[,1:5],heigh_clin_m,tmp[,6:ncol(tmp)])
  save(clinic.data,file=sprintf("%s/T0ClinicalParametersMicrObesja_All.RData",input.dir))
  write.table(clinic.data,sprintf("%s/T0ClinicalParametersMicrObesja_All.csv",input.dir),row.names=F,sep=";")
  
  
  #Genes data
  ge.data<-read.table(sprintf("%s/data.d.all.txt",input.dir),sep="\t",header=T)
  ge.annot.data<-read.table(sprintf("%s/sourceBatch236740_Entrez.txt",input.dir),sep="\t",header=T)
  
  
  
  Probe_Id<-ge.data[-1,1]
  Entrez_Gene_ID<-ge.data[-1,3]
  
  Gene_Names<-select(hgu95av2.db, Entrez_Gene_ID, "GENENAME", "ENTREZID")
  Entrez_Gene_ID_2<-select(hgu95av2.db, Probe_Id, "GENENAME", "PROBEID")
  #ReMOI
  #http://remoat.sysbiol.cam.ac.uk/search.php
  
  #Convetion Source:
  #http://www.shodhaka.com/cgi-bin/startbioinfo/simpleresources.pl?tn=Gene%20ID%20conversion
  
  
  ge.data.t0<-as.data.frame(ge.data[-1,which(substr(names(ge.data),9,10)=="0")])
  
  ge.data.t0.names<-gsub(".J0", "", names(ge.data.t0))
  ge.data.t0.names<-gsub(".", "-", ge.data.t0.names,fixed=TRUE)
  colnames(ge.data.t0)<-ge.data.t0.names
  rownames(ge.data.t0)<-Probe_Id
  ge.data.t0<-t(ge.data.t0)
  
  ge.data.t0<-ge.data.t0[-27,]

  
  n.na<-sapply(1:ncol(ge.data.t0),function(x) sum(is.na(ge.data.t0[,x]))) 
  ge.data.t0.nona<-ge.data.t0[,which(n.na==0)]

  save(ge.data.t0,file=sprintf("%s/GE_t0.RData",input.dir))
  
  write.table(ge.data.t0,file=sprintf("%s/GE_t0.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  
  write.table(unique(Entrez_Gene_ID),file=sprintf("%s/Entrez_Gene_ID.csv",input.dir),col.names=F,row.names=F,sep=",",quote=F)
  
  
  lm.classified.wfmh2.t0<-read.csv(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
  lm.classified.adfm.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
  lm.classified.adwh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
  lm.classified.adfmh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))
  lm.classified.adh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir))
  
  lm.classified.adbmifma.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir))
  
  lw.classified.wfmh2.t0<-read.csv(sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
  lw.classified.adfm.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
  lw.classified.adwh2.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
  lw.classified.adfmh2.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))
  lw.classified.adh2.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir))
  
  lw.classified.adbmifma.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-Fatmass-Age_points.csv",output.dir))
  
  
  tmp<-clinic.data[sapply(1:nrow(ge.data.t0), function(x) which(clinic.data$Names == lm.classified.wfmh2.t0$Name)),]
  
  ge.lm.classified.wfmh2.t0<-lm.classified.wfmh2.t0[lm.classified.wfmh2.t0$Name %in% rownames(ge.data.t0),]
  ge.lm.classified.adfm.t0<-lm.classified.adfm.t0[lm.classified.adfm.t0$Name %in% rownames(ge.data.t0),]
  ge.lm.classified.adwh2.t0<-lm.classified.adwh2.t0[lm.classified.adwh2.t0$Name %in% rownames(ge.data.t0),]
  ge.lm.classified.adfmh2.t0<-lm.classified.adfmh2.t0[lm.classified.adfmh2.t0$Name %in% rownames(ge.data.t0),]
  ge.lm.classified.adh2.t0<-lm.classified.adh2.t0[lm.classified.adh2.t0$Name %in% rownames(ge.data.t0),]
  
  ge.lm.classified.adbmifma.t0<-lm.classified.adbmifma.t0[lm.classified.adbmifma.t0$Name %in% rownames(ge.data.t0),]
  
  ge.lw.classified.wfmh2.t0<-lw.classified.wfmh2.t0[lw.classified.wfmh2.t0$Name %in% rownames(ge.data.t0),]
  ge.lw.classified.adfm.t0<-lw.classified.adfm.t0[lw.classified.adfm.t0$Name %in% rownames(ge.data.t0),]
  ge.lw.classified.adwh2.t0<-lw.classified.adwh2.t0[lw.classified.adwh2.t0$Name %in% rownames(ge.data.t0),]
  ge.lw.classified.adfmh2.t0<-lw.classified.adfmh2.t0[lw.classified.adfmh2.t0$Name %in% rownames(ge.data.t0),]
  ge.lw.classified.adh2.t0<-lw.classified.adh2.t0[lw.classified.adh2.t0$Name %in% rownames(ge.data.t0),]
  
  ge.lw.classified.adbmifma.t0<-lw.classified.adbmifma.t0[lw.classified.adbmifma.t0$Name %in% rownames(ge.data.t0),]

  
  write.table(ge.lm.classified.wfmh2.t0,file=sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  write.table(ge.lm.classified.adfm.t0,file=sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  write.table(ge.lm.classified.adwh2.t0,file=sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  write.table(ge.lm.classified.adfmh2.t0,file=sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  write.table(ge.lm.classified.adh2.t0,file=sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  
  write.table(ge.lm.classified.adbmifma.t0,file=sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  
  write.table(ge.lw.classified.wfmh2.t0,file=sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  write.table(ge.lw.classified.adfm.t0,file=sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  write.table(ge.lw.classified.adwh2.t0,file=sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  write.table(ge.lw.classified.adfmh2.t0,file=sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  write.table(ge.lw.classified.adh2.t0,file=sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  
  write.table(ge.lw.classified.adbmifma.t0,file=sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir),col.names=T,row.names=F,sep=",",quote=F)
  
  
  
  
  
}



#Annotation
#"diam_clin" = adipocyte diametre
#"tour.de.taille_clin"= Waist
#"MG.kg_clin" = Fat mas en KG (absolu)
#"MG.percent_clin" = Fat mass/body weight %
#"Taille"= Height


########################################
# PLOT CORRELATION AND ADD LOWESS LINE #
########################################

attach(clinic.data)

#t0 - "Adipocyte diameter" vs "Fat mass" by LOWESS
reg.plot.mc(fatmass_clin_kg,diam_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),         
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="Fat mass",ylab="Adipocyte diameter",
            legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))

#t0 - "Waist" vs "Fat mass/Height*Height" by LOWESS
x=fatmass_clin_kg/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,waist.circumference_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="Fat mass/Height*Height",ylab="Waist",
            legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))



#t0 - "Adipocyte diameter" vs "Waist/Height*Height" by LOWESS
x=waist.circumference_clin_cm/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,diam_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="Waist/Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))


#t0 - "Adipocyte diameter" vs "Fat mass/Height*Height" by LOWESS
x=fatmass_clin_kg/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,diam_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="Fat mass/Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))


#t0 - "Adipocyte diameter" vs "Height*Height" by LOWESS
x=heigh_clin_m*heigh_clin_m
reg.plot.mc(x,diam_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir))


#t0 - "Adipocyte diameter" vs "BMI~fatmass~Age" by LOWESS
x1<-fitted(loess(clinic.data$BMI_clin_kg.per.m~clinic.data$fatmass_clin_kg))
x<-fitted(loess(x1~clinic.data$Age_years))
reg.plot.mc(x,diam_clin_cm,
            type="lowess",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - LOWESS",
            xlab="BMI~fatmass~Age",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-fatmass-Age.pdf",output.dir),
            pointsfile=sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-fatmass-Age_points.csv",output.dir))



########################################
# PLOT CORRELATION AND ADD LINEAR LINE #
########################################


#t0 - "Adipocyte diameter" vs "Fat mass" by Linear Regression
reg.plot.mc(fatmass_clin_kg,diam_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),         
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Fat mass",ylab="Adipocyte diameter",
            legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))

#t0 - "Waist" vs "Fat mass/Height*Height" by Linear Regression
x=fatmass_clin_kg/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,waist.circumference_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Fat mass/Height*Height",ylab="Waist",
            legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))



#t0 - "Adipocyte diameter" vs "Waist/Height*Height" by Linear Regression
x=waist.circumference_clin_cm/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,diam_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Waist/Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))




#t0 - "Adipocyte diameter" vs "Fat mass/Height*Height" by Linear Regression
x=fatmass_clin_kg/(heigh_clin_m*heigh_clin_m)
reg.plot.mc(x,diam_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Fat mass/Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))

#t0 - "Adipocyte diameter" vs "Height*Height" by Linear Regression
x=heigh_clin_m*heigh_clin_m
reg.plot.mc(x,diam_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="Height*Height",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Height2.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir))



#t0 - "Adipocyte diameter" vs "BMI~fatmass" by Linear Regression
x<-fitted(loess(clinic.data$BMI_clin_kg.per.m~clinic.data$fatmass_clin_kg))
reg.plot.mc(x,diam_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="BMI~fatmass",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass_points.csv",output.dir))


#t0 - "Adipocyte diameter" vs "BMI~fatmass" by Linear Regression
x1<-fitted(loess(clinic.data$BMI_clin_kg.per.m~clinic.data$fatmass_clin_kg))
x<-fitted(loess(x1~clinic.data$Age_years))
reg.plot.mc(x,diam_clin_cm,
            type="lm",
            pch=ifelse(Sexe=="M", 0, 1),
            subjects=as.vector(Names),
            title="CORRELATION (T0) - Linear Regression",
            xlab="BMI~fatmass~Age",ylab="Adipocyte diameter",
            legend.topleft=list(title="SHAPE",pch=c(1,0),label=c("Female","Male"),col=c("black","black")),
            imgfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age.pdf",output.dir),
            pointsfile=sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir))



########################################################
# BOXPLOT FOR CLASSIFIED GROUPS WITH LINEAR REGRESSION #
########################################################

lm.classified.wfmh2.t0<-read.csv(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
lm.classified.adfm.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
lm.classified.adwh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
lm.classified.adfmh2.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))

lm.classified.adbmifm.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass_points.csv",output.dir))
lm.classified.adbmifma.t0<-read.csv(sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir))

pdf(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.adh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()




pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lm.classified.adbmifm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()










pdf(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lm.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lm.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lm.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lm.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()












pdf(sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lm_t0_Waist_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()

pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.adh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()


pdf(sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lm.classified.adbmifma.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()
















########################################################
# BOXPLOT FOR CLASSIFIED GROUPS WITH LOWESS REGRESSION #
########################################################


lw.classified.wfmh2.t0<-read.csv(sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
lw.classified.adfm.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
lw.classified.adwh2.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
lw.classified.adfmh2.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))

lw.classified.adbmifm.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-Fatmass_points.csv",output.dir))
lw.classified.adbmifma.t0<-read.csv(sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-Fatmass-Age_points.csv",output.dir))


pdf(sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.adbmifm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()


pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="t",class=lw.classified.adbmifma.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()










pdf(sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lw.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lw.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lw.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="mwu",class=lw.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()












pdf(sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/lw_t0_Waist_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()


pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.adh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(3,5:ncol(clinic.data)),function(x) boxplot.class.mc(data=clinic.data,x,type="auto",class=lw.classified.adbmifma.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()













#######################
# GENE EXPRESSION     #
#######################


#ge.data.t0.nona<-na.omit(ge.data.t0)
n.na<-sapply(1:ncol(ge.data.t0),function(x) sum(is.na(ge.data.t0[,x]))) 
ge.data.t0.nona<-ge.data.t0[,which(n.na==0)]

########################################################
# BOXPLOT FOR CLASSIFIED GROUPS WITH LINEAR REGRESSION #
########################################################

ge.lm.classified.wfmh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lm.classified.adfm.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
ge.lm.classified.adwh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
ge.lm.classified.adfmh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lm.classified.adh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir))

ge.lm.classified.adbmifma.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir))



pdf(sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="t",class=ge.lm.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="t",class=ge.lm.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="t",class=ge.lm.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="t",class=ge.lm.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()











pdf(sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="mwu",class=ge.lm.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="mwu",class=ge.lm.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="mwu",class=ge.lm.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="mwu",class=ge.lm.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()












pdf(sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lm.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lm.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lm.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lm.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lm.classified.adh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()




pdf(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lm.classified.adbmifma.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()















########################################################
# BOXPLOT FOR CLASSIFIED GROUPS WITH LOWESS REGRESSION #
########################################################


ge.lw.classified.wfmh2.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lw.classified.adfm.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
ge.lw.classified.adwh2.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
ge.lw.classified.adfmh2.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lw.classified.adbmifma.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir))


pdf(sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(203:205,function(x) boxplot.class.mc(data=ge.data.t0,x,type="t",class=as.factor(ge.lw.classified.wfmh2.t0[,4]),xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="t",class=ge.lw.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="t",class=ge.lw.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="t",class=ge.lw.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()


#03/07/2013
pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lw.classified.adh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_t-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="t",class=ge.lw.classified.adbmifma.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()









pdf(sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="mwu",class=ge.lw.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="mwu",class=ge.lw.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="mwu",class=ge.lw.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_mwu-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="mwu",class=ge.lw.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()












pdf(sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile=sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lw.classified.wfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lw.classified.adfm.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lw.classified.adwh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()



pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lw.classified.adfmh2.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()


pdf(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test_boxplot.pdf",output.dir))
outfile<-sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(1:ncol(ge.data.t0.nona),function(x) boxplot.class.mc(data=ge.data.t0.nona,x,type="auto",class=ge.lw.classified.adbmifma.t0[,4],xlab="Risk (0=low, 1=high)",outfile=outfile))
dev.off()





data<-clinic.data[,-c(1,2,4)]

plot(data[,1],data[,3])

le<-loess(data[,3]~data[,1])
y.predict<-predict(le,data[,3])

plot(data[,3],fitted(le))


ncol(data)


x<-fitted(loess(clinic.data$BMI_clin_kg.per.m~clinic.data$fatmass_clin_kg))

y<-fitted(loess(clinic.data$diam_clin_cm~clinic.data$fatmass_clin_kg))

plot(clinic.data$diam_clin_cm,x)

plot(clinic.data$BMI_clin_kg.per.m,clinic.data$fatmass_clin_kg)



#07/07

genes.updown<-function(genes,prode.entrez,high,low,stars,fileprefix,output.dir){
  m<-mean(as.numeric(as.matrix(genes)))
  res<-as.data.frame(rbind(c("Probe_ID","High_Risk_Group","Low_Risk_Group"),cbind(colnames(genes),ifelse(colMeans(genes[high,])>m,1,0),ifelse(colMeans(genes[low])>m,1,0))))
  
  #high.up
  ge.high.up<-res[rownames(res)%in%stars & res[,2]==1,]
  high.up<-t(genes[,colnames(genes)%in%as.character(ge.high.up[,1])])
  rownames(high.up)<-prode.entrez[prode.entrez[,1]%in%rownames(high.up),2]
  high.up<-high.up[na.omit(rownames(high.up)),]
  
  #high.up
  ge.high.up<-res[rownames(res)%in%stars & res[,2]==1,]
  high.up<-t(genes[,colnames(genes)%in%as.character(ge.high.up[,1])])
  rownames(high.up)<-prode.entrez[prode.entrez[,1]%in%rownames(high.up),2]
  high.up<-high.up[na.omit(rownames(high.up)),]
  
  #high.up
  ge.high.down<-res[rownames(res)%in%stars & res[,2]==0,]
  high.down<-t(genes[,colnames(genes)%in%as.character(ge.high.down[,1])])
  rownames(high.down)<-prode.entrez[prode.entrez[,1]%in%rownames(high.down),2]
  high.down<-high.down[na.omit(rownames(high.down)),]
  
  #low.up
  ge.low.up<-res[rownames(res)%in%stars & res[,3]==1,]
  low.up<-t(genes[,colnames(genes)%in%as.character(ge.low.up[,1])])
  rownames(low.up)<-prode.entrez[prode.entrez[,1]%in%rownames(low.up),2]
  low.up<-low.up[na.omit(rownames(low.up)),]  
  
  #low.down
  ge.low.down<-res[rownames(res)%in%stars & res[,3]==0,]
  low.down<-t(genes[,colnames(genes)%in%as.character(ge.low.down[,1])])
  rownames(low.down)<-prode.entrez[prode.entrez[,1]%in%rownames(low.down),2]
  low.down<-low.down[na.omit(rownames(low.down)),]
  
  high.refs<-c(rownames(high.up),rownames(high.down))
  low.refs<-c(rownames(low.up),rownames(low.down))

  high.up.frame<-as.data.frame(high.up)
  colnames(high.up.frame)[1]<-"GeneID"
  high.down.frame<-as.data.frame(high.down)
  colnames(high.down.frame)[1]<-"GeneID"
  
  low.up.frame<-as.data.frame(low.up)
  colnames(low.up.frame)[1]<-"GeneID"  
  low.down.frame<-as.data.frame(low.down)
  colnames(low.down.frame)[1]<-"GeneID"
  

  write.table(high.up,file=sprintf("%s/%s.high.up.genes.txt",output.dir,fileprefix),col.names=F,row.names=T,sep="\t",quote=F)
  write.table(high.down,file=sprintf("%s/%s.high.down.genes.txt",output.dir,fileprefix),col.names=F,row.names=T,sep="\t",quote=F)
  write.table(low.up,file=sprintf("%s/%s.low.up.genes.txt",output.dir,fileprefix),col.names=F,row.names=T,sep="\t",quote=F)
  write.table(low.down,file=sprintf("%s/%s.low.down.genes.txt",output.dir,fileprefix),col.names=F,row.names=T,sep="\t",quote=F)
  write.table(high.refs,file=sprintf("%s/%s.high.refs.genes.txt",output.dir,fileprefix),col.names=F,row.names=F,sep="\t",quote=F)
  write.table(low.refs,file=sprintf("%s/%s.low.refs.genes.txt",output.dir,fileprefix),col.names=F,row.names=F,sep="\t",quote=F)
  write.table(res,file=sprintf("%s/updown.genes.%s.txt",output.dir,fileprefix),col.names=F,row.names=F,sep="\t",quote=F)
  
  
  high.up.frame<-as.data.frame(read.table(sprintf("%s/%s.high.up.genes.txt",output.dir,fileprefix),sep="\t",header=F))
  high.down.frame<-as.data.frame(read.table(sprintf("%s/%s.high.down.genes.txt",output.dir,fileprefix),sep="\t",header=F))
  low.up.frame<-as.data.frame(read.table(sprintf("%s/%s.low.up.genes.txt",output.dir,fileprefix),sep="\t",header=F))
  low.down.frame<-as.data.frame(read.table(sprintf("%s/%s.low.down.genes.txt",output.dir,fileprefix),sep="\t",header=F))
  
  
  fun(org="HS", two.lists=TRUE, up.frame=high.up.frame, down.frame=high.down.frame,
      genes.frame=NULL, restrict=TRUE, ref.list=high.ref.list, logged=FALSE,
      discriminant=TRUE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,
      annot.method="specificity", annot.details=TRUE,
      direct=FALSE, enriched=TRUE, fdr=NA, build.annot.net=TRUE,
      coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE,
      hard.th=0.8, soft.th=NA, topological = FALSE, keep.sign=FALSE, level=1,
      annot.clust.method="umilds", annot.prox.measure="dynamical",
      test.recovery=FALSE, test.robust=FALSE, replace.annot=NA,
      build.gene.net=TRUE, gene.clust.method="hclust", gene.net.details=TRUE,
      gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE,fileprefix=sprintf("%s.high.genes",fileprefix))
  
  fun(org="HS", two.lists=TRUE, up.frame=low.up, down.frame=low.down,
      genes.frame=NULL, restrict=TRUE, ref.list=prode.entrez[,1], logged=FALSE,
      discriminant=TRUE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,
      annot.method="specificity", annot.details=TRUE,
      direct=FALSE, enriched=TRUE, fdr=NA, build.annot.net=TRUE,
      coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE,
      hard.th=0.8, soft.th=NA, topological = FALSE, keep.sign=FALSE, level=1,
      annot.clust.method="umilds", annot.prox.measure="dynamical",
      test.recovery=FALSE, test.robust=FALSE, replace.annot=NA,
      build.gene.net=TRUE, gene.clust.method="hclust", gene.net.details=TRUE,
      gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE,fileprefix=sprintf("%s.low.genes",fileprefix))

  
  return(res)
}

ge.data.t0<-get(load(sprintf("%s/GE_t0.RData",input.dir)))

#ge.data.t0.nona<-na.omit(ge.data.t0)
n.na<-sapply(1:ncol(ge.data.t0),function(x) sum(is.na(ge.data.t0[,x]))) 
ge.data.t0.nona<-ge.data.t0[,which(n.na==0)]

ge.lm.classified.adfm.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
ge.lm.classified.adwh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
ge.lm.classified.adfmh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lm.classified.adh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir))
ge.lm.classified.adbmifma.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir))

ge.lm.classified.adfm.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_auto-test.csv",output.dir),sep=";")
ge.lm.classified.adwh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test.csv",output.dir),sep=";")
ge.lm.classified.adfmh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test.csv",output.dir),sep=";")
ge.lm.classified.adh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Height2_auto-test.csv",output.dir),sep=";")
ge.lm.classified.adbmifma.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test.csv",output.dir),sep=";")

ge.lm.classified.adfm.t0.autotest.star<-ge.lm.classified.adfm.t0.autotest[ge.lm.classified.adfm.t0.autotest[,4]!="",]
ge.lm.classified.adwh2.t0.autotest.star<-ge.lm.classified.adwh2.t0.autotest[ge.lm.classified.adwh2.t0.autotest[,4]!="",]
ge.lm.classified.adfmh2.t0.autotest.star<-ge.lm.classified.adfmh2.t0.autotest[ge.lm.classified.adfmh2.t0.autotest[,4]!="",]
ge.lm.classified.adh2.t0.star<-ge.lm.classified.adh2.t0.autotest[ge.lm.classified.adh2.t0.autotest[,4]!="",]
ge.lm.classified.adbmifma.t0.autotest.star<-ge.lm.classified.adbmifma.t0.autotest[ge.lm.classified.adbmifma.t0.autotest[,4]!="",]




ge.data<-read.table(sprintf("%s/data.d.all.txt",input.dir),sep="\t",header=T)
ge.prode.entrez<-ge.data[,2:3]

updown.genes.ge.lm.classified.wfmh2.t0<-genes.updown(genes=ge.data.t0.nona,
                                                     prode.entrez=ge.prode.entrez,
                                                     high=which(ge.lm.classified.adfm.t0[,4]==0),
                                                     low=which(ge.lm.classified.adfm.t0[,4]==1),
                                                     stars=ge.lm.classified.adfm.t0.autotest.star[,1],
                                                     fileprefix="ge.lm.classified.adfm.t0",
                                                     output.dir=output.dir)
updown.genes.ge.lm.classified.adwh2.t0<-genes.updown(genes=ge.data.t0.nona,
                                                     prode.entrez=ge.prode.entrez,
                                                     high=which(ge.lm.classified.adwh2.t0[,4]==0),
                                                     low=which(ge.lm.classified.adwh2.t0[,4]==1),
                                                     stars=ge.lm.classified.adfm.t0.autotest.star[,1],
                                                     fileprefix="ge.lm.classified.adwh2.t0",
                                                     output.dir=output.dir)
updown.genes.ge.lm.classified.adfmh2.t0<-genes.updown(genes=ge.data.t0.nona,
                                                     prode.entrez=ge.prode.entrez,
                                                     high=which(ge.lm.classified.adfmh2.t0[,4]==0),
                                                     low=which(ge.lm.classified.adfmh2.t0[,4]==1),
                                                     stars=ge.lm.classified.adfm.t0.autotest.star[,1],
                                                     fileprefix="ge.lm.classified.adfmh2.t0",
                                                     output.dir=output.dir)
updown.genes.ge.lm.classified.adh2.t0<-genes.updown(genes=ge.data.t0.nona,
                                                    prode.entrez=ge.prode.entrez,
                                                     high=which(ge.lm.classified.adh2.t0[,4]==0),
                                                     low=which(ge.lm.classified.adh2.t0[,4]==1),
                                                     stars=ge.lm.classified.adfm.t0.autotest.star[,1],
                                                     fileprefix="ge.lm.classified.adh2.t0",
                                                     output.dir=output.dir)
updown.genes.ge.lm.classified.adbmifma.t0<-genes.updown(genes=ge.data.t0.nona,
                                                        prode.entrez=ge.prode.entrez,
                                                     high=which(ge.lm.classified.adbmifma.t0[,4]==0),
                                                     low=which(ge.lm.classified.adbmifma.t0[,4]==1),
                                                     stars=ge.lm.classified.adfm.t0.autotest.star[,1],
                                                     fileprefix="ge.lm.classified.adbmifma.t0",
                                                     output.dir=output.dir)




updown.genes.ge.lm.classified.adfm.t0.autotest.star<-updown.genes.ge.lm.classified.adbmifma.t0[updown.genes.ge.lm.classified.adbmifma.t0[,1]%in%ge.lm.classified.adfm.t0.autotest.star[,1],]



updown.genes.ge.lm.classified.adfm.t0.autotest.star<-genes.updown(genes=ge.data.t0.nona[,colnames(ge.data.t0.nona)[which(colnames(ge.data.t0.nona)%in%ge.lm.classified.adfm.t0.autotest.star[,1])]],
                                                        high=which(ge.lm.classified.adfm.t0[,4]==0),
                                                        low=which(ge.lm.classified.adfm.t0[,4]==1),
                                                        filename=sprintf("%s/updown.genes.ge.lm.classified.adfm.t0.autotest.star.csv",output.dir))

updown.genes.ge.lm.classified.adwh2.t0.autotest.star<-genes.updown(genes=ge.data.t0.nona[,colnames(ge.data.t0.nona)[which(colnames(ge.data.t0.nona)%in%ge.lm.classified.adwh2.t0.autotest.star[,1])]],
                                                                  high=which(ge.lm.classified.adwh2.t0[,4]==0),
                                                                  low=which(ge.lm.classified.adwh2.t0[,4]==1),
                                                                  filename=sprintf("%s/updown.genes.ge.lm.classified.adwh2.t0.autotest.star.csv",output.dir))



plot(as.numeric(updown.genes.ge.lm.classified.wfmh2.t0[,2]),
     as.numeric(updown.genes.ge.lm.classified.wfmh2.t0[,3]),
     xlab="Mean expression level of High Risk Group (HRG)",
     ylab="Mean expression level of Low Risk Group (LRG)",
     main="Distribution of gene expression level between groups \n Classification by 'Waist vs Fat mass/Heigth^2'")

plot(as.numeric(updown.genes.ge.lm.classified.adbmifma.t0[,2]),
     as.numeric(updown.genes.ge.lm.classified.adbmifma.t0[,3]),
     xlab="Mean expression level of High Risk Group (HRG)",
     ylab="Mean expression level of Low Risk Group (LRG)",
     main="Distribution of gene expression level between groups \n Classification by 'Waist vs Fat mass/Heigth^2'")



plot(as.numeric(as.character(updown.genes.ge.lm.classified.adbmifma.t0.autotest.star[,2])),
     as.numeric(as.character(updown.genes.ge.lm.classified.adbmifma.t0.autotest.star[,3])),
     xlab="Mean expression level of High Risk Group (HRG)",
     ylab="Mean expression level of Low Risk Group (LRG)",
     xlim=c(0,3000),
     ylim=c(0,3000),
     main="Distribution of gene expression level between groups \n Classification by 'Waist vs Fat mass/Heigth^2'")

plot(as.numeric(as.character(ge.lm.classified.adbmifma.t0.autotest.star.not.star[,2])),
     as.numeric(as.character(ge.lm.classified.adbmifma.t0.autotest.star.not.star[,3])),
     xlab="Mean expression level of High Risk Group (HRG)",
     ylab="Mean expression level of Low Risk Group (LRG)",
     xlim=c(0,3000),
     ylim=c(0,3000),
     main="Distribution of gene expression level between groups \n Classification by 'Waist vs Fat mass/Heigth^2'")




high.ref.list<-read.table(sprintf("%s/Entrez_Gene_ID.csv",input.dir),sep="\t",header=F)
high.ref.list<-data.frame(high.ref.list)
rownames(high.ref.list)<-high.ref.list[,1]

high.up.frame<-as.data.frame(read.table(sprintf("%s/ge.lm.classified.adfm.t0.highup.genes.txt",output.dir),sep="\t",header=F))
high.down.frame<-as.data.frame(read.table(sprintf("%s/ge.lm.classified.adfm.t0.highdown.genes.txt",output.dir),sep="\t",header=F))


fun(org="HS", two.lists=TRUE, up.frame=high.up.frame, down.frame=high.down.frame,
    genes.frame=NULL, restrict=TRUE, ref.list=high.ref.list, logged=FALSE,
    discriminant=TRUE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,
    annot.method="specificity", annot.details=TRUE,
    direct=FALSE, enriched=TRUE, fdr=NA, build.annot.net=TRUE,
    coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE,
    hard.th=0.8, soft.th=NA, topological = FALSE, keep.sign=FALSE, level=1,
    annot.clust.method="umilds", annot.prox.measure="dynamical",
    test.recovery=FALSE, test.robust=FALSE, replace.annot=NA,
    build.gene.net=TRUE, gene.clust.method="hclust", gene.net.details=TRUE,
    gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE,fileprefix=sprintf("%s/%s.high.genes",output.dir,"ge.lm.classified.adfm.t0"))

fun(org="HS", two.lists=TRUE, up.frame=low.up, down.frame=low.down,
    genes.frame=NULL, restrict=TRUE, ref.list=prode.entrez[,1], logged=FALSE,
    discriminant=TRUE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,
    annot.method="specificity", annot.details=TRUE,
    direct=FALSE, enriched=TRUE, fdr=NA, build.annot.net=TRUE,
    coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE,
    hard.th=0.8, soft.th=NA, topological = FALSE, keep.sign=FALSE, level=1,
    annot.clust.method="umilds", annot.prox.measure="dynamical",
    test.recovery=FALSE, test.robust=FALSE, replace.annot=NA,
    build.gene.net=TRUE, gene.clust.method="hclust", gene.net.details=TRUE,
    gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE,fileprefix=sprintf("%s/%s.low.genes",output.dir,"ge.lm.classified.adfm.t0"))

