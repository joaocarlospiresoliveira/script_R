library(raster, quietly = T)
library(virtualspecies)

# Chamar variáveis ambientais (aqui é o Worldclim)
worldclim <- getData("worldclim", var = "bio", res = 10)
worldclim
names(worldclim)
worldclim<-list.files('./', pattern = '.tif', full.names = T)
worldclim<-stack(worldclim)
# The first four layers
worldclim[[1:4]]#acessar subconjuto de variaveis
# Layers bio1 and bio12 (annual mean temperature and annual precipitation)
worldclim[[c("bio1", "bio12")]]
# Plotting bio1 and bio12 and bio3
par(oma = c(0.1, 0.1, 0.1, 2.1))
plot(worldclim[[c("bio1", "bio12")]])
# Suitability of the environment for bio1 = 15 °C
dnorm(x = 150, mean = 250, sd = 50)

# my.parameters <- formatFunctions(bio1 = c(fun = 'dnorm', mean = 250, sd = 50),
#                                  bio12 = c(fun = 'dnorm', mean = 4000, sd = 2000))
# 
w1<-na.exclude(values(worldclim[[1]]))
w2<-na.exclude(values(worldclim[[2]]))
w3<-na.exclude(values(worldclim[[3]]))
w4<-na.exclude(values(worldclim[[4]]))
w5<-na.exclude(values(worldclim[[5]]))
w6<-na.exclude(values(worldclim[[6]]))
w7<-na.exclude(values(worldclim[[7]]))
w8<-na.exclude(values(worldclim[[8]]))
w9<-na.exclude(values(worldclim[[9]]))
w10<-na.exclude(values(worldclim[[10]]))
w11<-na.exclude(values(worldclim[[11]]))
w12<-na.exclude(values(worldclim[[12]]))
w13<-na.exclude(values(worldclim[[13]]))
w14<-na.exclude(values(worldclim[[14]]))
w15<-na.exclude(values(worldclim[[15]]))
w16<-na.exclude(values(worldclim[[16]]))
w17<-na.exclude(values(worldclim[[17]]))
w18<-na.exclude(values(worldclim[[18]]))
w19<-na.exclude(values(worldclim[[19]]))
my.stack<-rescale(worldclim)




my.parameters <- formatFunctions(bio1 = c(fun = "logisticFun", 
                                          alpha = -12.7, beta = 68),
                                bio2 = c(fun = "linearFun", 
                                          a = -0.03, b = 191.2),
                                bio3 = c(fun = "dnorm", 
                                          mean = 86.4, sd = 19.1),
                                bio4 = c(fun = "logisticFun", 
                                          alpha = 2198.5, beta = 11381.4))                                
# Gerando uma especie virtual
my.first.species <-
  generateSpFromFun(
    raster.stack = worldclim[[c(
      "bio1",
      "bio2",
      'bio3',
      'bio4'
    )]],
    parameters = my.parameters,
    formula = NULL,
    species.type = "multiplicative",
    rescale = FALSE,
    plot = TRUE
  )

# Convertendo para P/A

my.first.species
pa1$pa.raster <- convertToPA(my.first.species, plot = TRUE, PA.method = "probability", beta = 0.5,
                   alpha = -0.01)


# Fazendo amostragem da VE 
presence.points <- sampleOccurrences(pa1,
                                     n = 500, # The number of points to sample
                                     type = "presence only")
vc.pa<-pa1$pa.raster
point<-data.frame(presence.points$sample.points)
vc<-(my.first.species$suitab.raster)
writeRaster(
  vc.pa,
  filename = "Current Climate_RF.bin.tif",
  formato = "GTiff",
  overwrite = TRUE
)



# Virtual Especies com PCA
my.stack <- worldclim

my.pca.species <- generateSpFromPCA(raster.stack = my.stack,niche.breadth = "any",axes = c(1, 2),
                                    sample.points = F, nb.points = NULL)

par(mfrow = c(1, 3))
plotResponse(my.pca.species, axes = c(1, 2))
plotResponse(my.pca.species, axes = c(1, 3))
plotResponse(my.pca.species, axes = c(2, 3))

############# Você também pode usar as funções básicas fornecidas com o pacote: #####
#### Função linear: formatFunctions(bio1 = c(fun = 'linearFun', a = 1, b = 0)) \ [f (bio1) = a * bio1 + b \]
#### Função quadratica : formatFunctions(bio1 = c(fun = 'quadraticFun', a = 1, b = 2, c = 0)) \ [f (bio1) = a \ times bio1 ^ 2 + b \ times bio1 + c \]
#### Função logística: formatFunctions(bio1 = c(fun = 'logisticFun', beta = 150, alpha = -5)) \ [f (bio1) = \ frac {1} {1 + e ^ {\ frac {bio1 - \ beta} {\ alpha}}} \]

linear.function <- function(x, a, b)
{
  a*x + b
}

linear.function(x = 0.5, a = 2, b = 1)