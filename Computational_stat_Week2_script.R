# Scripts Computational Statistics Week 2

# Problème 3.1 
# a) Traçons les intégrandes  et utilisons une méthode d'intégration de Monte Carlo pour déterminer les intégrales

nt <- 1000
x<-c(0,5,10,20,50)
theta <- seq(0, 60, length = nt)
plot(theta,theta/(1+theta^2)*exp(-(x[1]-(theta))^2/(2)),type="l",main="Tracé de g1(theta)",col="blue",ylab="g1(theta)")
lines(theta,theta/(1+theta^2)*exp(-(x[2]-(theta))^2/(2)),col="red")
lines(theta,theta/(1+theta^2)*exp(-(x[3]-(theta))^2/(2)),col="green")
lines(theta,theta/(1+theta^2)*exp(-(x[4]-(theta))^2/(2)),col="purple")
lines(theta,theta/(1+theta^2)*exp(-(x[5]-(theta))^2/(2)),col="brown")





legend("right", # la position sur le graphique
       c("h1 pour x=0" ,"h1 pour x=5","h1 pour x=10","h1 pour x=20","h1 pour x=50"), # le texte pour chaque courbe
       col=c("blue", "red","green","purple","brown"), # La couleur de chaque courbe
       lwd=c(1,3,1), # L'épaisseur de chaque courbe
       lty=c(1,1,3) # Le type de trait de chaque courbe
)


nt <- 1000
x<-c(0,5,10,20,50)
theta <- seq(0, 60, length = nt)
plot(theta,1/(1+theta^2)*exp(-(x[1]-(theta))^2/(2)),type="l",main="Tracé de g2(theta)",col="blue",ylab="g2(theta)")
lines(theta,1/(1+theta^2)*exp(-(x[2]-(theta))^2/(2)),col="red")
lines(theta,1/(1+theta^2)*exp(-(x[3]-(theta))^2/(2)),col="green")
lines(theta,1/(1+theta^2)*exp(-(x[4]-(theta))^2/(2)),col="purple")
lines(theta,1/(1+theta^2)*exp(-(x[5]-(theta))^2/(2)),col="brown")





legend("right", # la position sur le graphique
       c("h1 pour x=0" ,"h1 pour x=5","h1 pour x=10","h1 pour x=20","h1 pour x=50"), # le texte pour chaque courbe
       col=c("blue", "red","green","purple","brown"), # La couleur de chaque courbe
       lwd=c(1,3,1), # L'épaisseur de chaque courbe
       lty=c(1,1,3) # Le type de trait de chaque courbe
)


delta1=0
delta2=0
n=130000
x=5
k=1
theta <- c(rnorm(n, mean = x, sd = 1))
repeat {
  delta1=delta1 + (theta[k])/(1+theta[k]^2)
  delta2=delta2 + (1)/(1+theta[k]^2)
  k<-k+1
  if (k==n) break
}
delta1=delta1/n
delta2=delta2/n


#b) Nous allons contrôler la précision de notre estimation à 95% à partir de la variance empirique pour avoir trois chiffres significatifs

vh1 = 0
vh2 = 0
k=1
repeat {
  vh1=vh1 + ((theta[k])/(1+theta[k]^2)-delta1)
  vh2=vh2 + ((1)/(1+theta[k]^2)-delta2)
  k<-k+1
  if (k==n) break
}
vh1=vh1/(n-1)
vh2=vh2/(n-1)

# On calcule les intervalles de confiance à n fixé IC1=[a1,b1] , IC2=[a2,b2]

a1=delta1-1.96*(sqrt(vh1)/sqrt(n))
b1=delta1+1.96*(sqrt(vh1)/sqrt(n))
a2=delta2-1.96*(sqrt(vh2)/sqrt(n))
b2=delta2+1.96*(sqrt(vh2)/sqrt(n))
a1
b1
a2
b2
c=a1/b2
d=b1/a2
c
d
d-c
# On crée la fonction associé au numérateur
f1 = function(y){
  return(y/(1+y^2))
}
# On crée la fonction associé au dénominateur
f2 = function(y){
  return(1/(1+y^2))
}

Nsim=10^4
y=rnorm(Nsim,x,1)

plot(cumsum(f1(y))/1:Nsim,type="l")
Nsim=10^4
y=rnorm(Nsim,x,1)
plot(cumsum(f2(y))/1:Nsim,type="l")

# Problème 5.1

#a) Utilisation de la fonction optimize de R

f = function(x){
  return((cos(50*x)+sin(20*x))^2)
}

x51=seq(from=0,to=pi/5,by=0.001)
plot(x51,f(x51),type="l",col="blue",main="Tracé de f(x)",xlab="x",ylab="f(x)")

optimize(f,x51,lower=min(x51),upper=max(x51),maximum=TRUE)

recherche_sto = function(f,U){
  temp<-c()
  k=1
  repeat {
    if (k==length(U)) break
    temp[k]=f(U[k])
    k<-k+1
  }
  result1=max(temp) #On cherche le maximum de vraisemblance de f
  rang=which(temp==result1) #Pour obtenir la valeur de notre échantillon uniforme maximisant la vraisemblance f
  return(U[rang]) # La fonction retourne la valeur qui maximise la vraisemblance  
}

Uf = runif(1000,min=0,max=pi/5)

recherche_sto(f,Uf)

# Problème 5.2

#a) On considère la vraisemblance  de  $Z=min(X,Y)$ avec $X \thicksim \mathcal{N}(\theta,\sigma^2)$ et $Y  \thicksim \mathcal{N}(\mu,\tau^2)$ et on cherche à estimer le maximum par rapport à $\theta$ puis à $\tau$ par deux méthodes différentes

pas = 0.1
theta = seq(from=-10, to=10, by = pas)
likelyztheta = function(theta){
  z=0.05
  sigma=1
  tau=1
  mu=0.5   
  return((1 - pnorm((z - theta) / sigma)) *
           (1 / tau) * dnorm((z - mu) / tau) +
           (1 - pnorm((z - mu) / tau)) *
           (1 / sigma) * dnorm((z - theta) / sigma))
}
plot(theta, likelyztheta(theta),
     type = "l", lwd = 3,
     main = "Vraisemblance de Z par rapport à theta",
     xlab = "theta", ylab = "Lz(theta)", col = "red")

optimize(likelyztheta,theta,lower=min(theta),upper=max(theta),maximum=TRUE)$maximum

tau = seq(from=-15, to=15, by = pas)
likelyztau = function(tau){
  z=0.05
  sigma=1
  mu=0.5
  theta=0.585200266251986
  return((1 - pnorm((z - theta) / sigma)) *
           (1 / tau) * dnorm((z - mu) / tau) +
           (1 - pnorm((z - mu) / tau)) *
           (1 / sigma) * dnorm((z - theta) / sigma))
}
plot(tau, likelyztau(tau),
     type = "l", lwd = 3,
     main = "Vraisemblance de Z par rapport à tau",
     xlab = "tau", ylab = "Lz(tau)", col = "red")

optimize(likelyztau,mu,lower=min(tau),upper=max(tau),maximum=TRUE)$maximum

Utheta = runif(10000,min=-5,max=5)
recherche_sto(likelyztheta,Utheta)
Utau = runif(10000,min=-15,max=15)
recherche_sto(likelyztau,Utau)

#b) On considère ici la vraisemblance de   $W \thicksim p\mathcal{N}(\theta,\sigma^2)+(1-p)\mathcal{N}(\mu,\tau^2)$ et on cherche à estimer le maximum par rapport à $\theta$ puis à $\tau$ par deux méthodes différentes :
theta2 = seq(from=-10, to=10, by = pas)
likelywtheta = function(theta){
  sigma=1
  mu=2
  tau=1
  p=0.5
  w=5
  
  return(p * dnorm(w, theta, sigma) +
           (1 - p) * dnorm(w, mu, tau))
}

plot(theta2, likelywtheta(theta2),
     type = "l", lwd = 1,
     main = "Vraisemblance de W par rapport à theta",
     xlab = "theta", ylab = "LW(theta)", col = "blue")

optimize(likelywtheta,theta2,lower=min(theta2),upper=max(theta2),maximum=TRUE)$maximum

tau2 = seq(from=0, to=20, by = pas)
likelywtau = function(tau){
  sigma=1
  mu=2
  theta=-2
  p=0.5
  w=5
  
  return(p * dnorm(w, theta, sigma) +
           (1 - p) * dnorm(w, mu, tau))
}

plot(tau2, likelywtau(tau2),
     type = "l", lwd = 1,
     main = "Vraisemblance de W par rapport à tau",
     xlab = "tau", ylab = "LW(tau)", col = "blue")

tau,tau2,lower=min(tau2),upper=max(tau2),maximum=TRUE)$maximum

Utheta = runif(10000,min=-5,max=5)

recherche_sto(likelywtheta,Utheta)
Utau = runif(10000,min=0,max=20)
recherche_sto(likelywtau,Utau)






