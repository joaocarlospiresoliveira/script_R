# install.packages("devtools", dep = T) # instalar pacote necessrio para instalar o 'spfilt'
devtools::install_github("diogosbr/spfilt", dep = T)
library(spfilt)
library(openxlsx)
spp<-read.table(file.choose(), header=T, sep=";")

dim(spp)
x<-filt(spp, inverted = TRUE, shape.municipios = NULL) #por padrao usa Municipios IBGE

