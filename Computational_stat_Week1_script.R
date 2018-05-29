# Scripts Computational Statistics Week 1

# Problème 2.1

# a) Création de 1000 variables aléatoires et tracé d'un histogramme 
par(mfrow = c(2,1))
X = runif(1000)
plot(X,main="Tracé du nuage de point issu de 1000 variables alétoires uniformes")
hist(X,freq=FALSE)


# b) Génération de n variables aléatoires et tracé des pairs pour autocorrélation
plot(X[1:999],X[2:1000],main="Tracé  X[k+1] et en fonction de X[k]")

Box.test(X, lag = 1)
Box.test (X, lag = 1, type = "Ljung")

# Problème 2.2 

# a) Création de 1000 variables aléatoires binomiales à partir de lois uniformes
p=0.2
n=25
P <-c()
k=0
pk=0
sk=0
while (sk<0.99999) {
  pk=choose(n,k)*p^k*(1-p)^(n-k)
  sk <- sk + pk
  P[k+1]=sk
  k <-k+1
}
P
U = runif(1000)

aleagen <- function(u,P) {
  # Permet d'associer à une probabilité u le rang k dans une liste P croissante tel que u<P[k]
  # La fonction retourne le k associé 
  k=1
  repeat {
    if (u<P[k]) break
    k<-k+1
  }
  return(k-1)
}

Bin <-c()
for (k in seq(1,length(U))) {
  a=U[k]
  Bin[k]=aleagen(a,P)
}

hist(Bin,main="Tracé de l'histogramme de 1000 variables aléatoires binomiales",freq=FALSE)

plot(1:15,dbinom(1:15,n,p),main="Tracé de la densité d'une binomiale pour chaque valeur possible",type="p")

hist(rbinom(1000,n,p),main="Histogramme de 1000 v.a binomiales avec le générateur de R",freq=FALSE)

# b) Création de 5000 variables aléatoires logarithmiques à partir de lois uniformes et sa "densité"

L <-c()
p=0.2
k=1
pk=0
lk=0
while (lk<0.99999) {
  pk=(-(1-p)^k)/(k*log(p))
  lk <- lk + pk
  L[k]=lk
  k <-k+1
}

U=runif((5000))

Xlog <-c()
for (k in seq(1,length(U))) {
  a=U[k]
  Xlog[k]=aleagen(a,L)
}

hist(Xlog,main="Histogramme de 5000 v.a logarithmiques",freq=FALSE)
plot(1:39,(-(1-p)^(1:39))/((1:39)*log(p)),type="p",main="Tracé de la densité de la loi pour chaque valeur possible")
