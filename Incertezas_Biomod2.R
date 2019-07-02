pres<-stack(projections.RF.all,projections.GBM.all,projections.CTA.all,
            projections.GLM.all,projections.GAM.all,projections.ANN.all,
            projections.SRE.all,projections.MARS.all,projections.FDA.all,
            projections.MAXENT.Phillips.all)

fut<-stack(projections.RF.all.2050.rcp45,projections.GBM.all.2050.rcp45,projections.CTA.all.2050.rcp45,
           projections.GLM.all.2050.rcp45,projections.GAM.all.2050.rcp45,projections.ANN.all.2050.rcp45,
           projections.SRE.all.2050.rcp45,projections.MARS.all.2050.rcp45,projections.FDA.all.2050.rcp45,
           projections.MAXENT.Phillips.all.2050.rcp45)

all<-stack(pres,fut)
all.point<-rasterToPoints(all)

RF.p1<-stack(subset(all, grep("RF", names(all))))
GBM.p1<-stack(subset(all, grep("GBM", names(all))))
CTA.p1<-stack(subset(all, grep("CTA", names(all))))
GLM.p1<-stack(subset(all, grep("GLM", names(all))))
GAM.p1<-stack(subset(all, grep("GAM", names(all))))
ANN.p1<-stack(subset(all, grep("ANN", names(all))))
SRE.p1<-stack(subset(all, grep("SRE", names(all))))
MARS.p1<-stack(subset(all, grep("MARS", names(all))))
FDA.p1<-stack(subset(all, grep("FDA", names(all))))
MAXENT.Phillips.p1<-stack(subset(all, grep("MAXENT.Phillips", names(all))))
suit.r<-stack(RF.p1,GBM.p1, CTA.p1, GLM.p1, GAM.p1, ANN.p1, SRE.p1, MARS.p1, FDA.p1, MAXENT.Phillips.p1)
suit.r<-(suit.r/1000)
st<-rasterToPoints(suit.r, na.rm = T)
dim(st)
st<-na.omit(st)
st1<-as.data.frame(suit.r, na.rm =T, xy=F) # suit
head(st1)
coord.1<-st[,1:2]


CA<-mean(stack(projections.2050.rcp45_1_CA,projections.2050.rcp45_2_CA))
CS<-mean(stack(projections.2050.rcp45_1_CS,projections.2050.rcp45_2_CS))
IP<-mean(stack(projections.2050.rcp45_1_IP,projections.2050.rcp45_2_IP))
GCM.stack<-rescale(stack(CA,CS,IP))
GCM.p<-rasterToPoints(GCM.stack)

GCMs<-(na.omit(GCM.p))
dim(GCMs)
##### Avaliando Incertezas ####
# id.1 <- rep(c("pres","fut"), each= nrow(st1))
# suit.1 <- cbind(st1[id.1=='pres',],st1[id.1=='fut',])
metodo.1 <-
  rep(rep(c("RF","GBM","CTA", "GLM","GAM","ANN","SRE","MARS","FDA","MAXENT.Phillips"),each = n.runs * n.conj.pa2), 
      times = 2)
tempo.1 <- rep(c("p","f"), each = n.runs*length(st1))
# idp<-rep(c("PR"), each = 10, times=2)
idf<-rep(c("PR","CA","CS","IP"), each= 2*10)
# GCM1<-rep(c("CA", "CS", "IP"), each = tempo.1)
GCM<-c(idf)
# head(suit)
resu.1<-NULL
for (l in 1:nrow(st1)) {
  suit.l.1 <- data.frame(st1 = as.numeric(st1[l,]), metodo.1 = metodo.1, tempo.1 = tempo.1, GCM = GCM)
  modelo.aov.1 <- summary.aov(lm(st1~metodo.1%in%tempo.1+tempo.1+tempo.1%in%GCM+GCM, data = suit.l.1))
  resu.1<-rbind(resu.1, modelo.aov.1[[1]][,3])
}
head(resu.1)
dim(resu.1)
colnames(resu.1)<-c("MStempo","GCM", "MSmetodo","MSGCM","resíduo")
resu.1.prop <- resu.1/rowSums(resu.1)
head(resu.1.prop)
summary(resu.1.prop)
# hist(resu.1.prop[,2])
resu.1.r <- rasterFromXYZ(cbind(coord.1, resu.1.prop))
plot(resu.1.r)
