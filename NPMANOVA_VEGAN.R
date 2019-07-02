data(dune)
data(dune.env)

eixos<-data.frame(data[3:16])
npm.result<-adonis(eixos ~ Grups, data=data, permutations=999, method = "euclidean", by='terms')
summary(npm.result)
aov_tab<-npm.result$aov.tab
adonis(dune ~ Management*A1, data=dune.env, permutations=999, by='terms')

### Example of use with strata, for nested (e.g., block) designs.

dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
dat
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
total <- Agropyron + Schizachyrium
library(lattice)
dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium)
mod <- metaMDS(Y)
plot(mod)
### Hulls show treatment
ordihull(mod, group=dat$NO3, show="0")
ordihull(mod, group=dat$NO3, show="10", col=3)
### Spider shows fields
ordispider(mod, group=dat$field, lty=3, col="red")

### Correct hypothesis test (with strata)
adonis(Y ~ NO3, data=dat, strata=dat$field, perm=1e3)

### Incorrect (no strata)
adonis(Y ~ NO3, data=dat, perm=1e3)
