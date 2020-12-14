#### Simulação - subestimação dados incompletos ####
## Grupo de estudo Modelo não Lineares UFLA
# Autor: Carlos Antônio Zarzar
# Data: 15/12/2020
# E-mail: carloszarzar_@hotmail.com
# Youtube: https://www.youtube.com/channel/UC0aJ_xgty6efYmCy-9PJUYA?view_as=subscriber
#----------------------#----------------------#----------------------#----------------------
# Script: Simulação do peso do camarão em um cultivo
# de baixa densidade utilizando dados limitadas a partir de uma 
# curva Micahelis Menten para crescimento (Artigo: 
# A generalized Michaelis-Menten equation for the analysis of growth
# López et al. 2000). 
# E o segundo objetivo  do script é
# modelar outras curvas a partir desses dados simulados.
# Michaelis-Menten equation
# w(t)= (w0*beta^kappa + w1*x^kappa)/(beta^kappa + x^kappa) ::: MM growth curve
#----------------------#----------------------#----------------------#----------------------
rm(list = ls())
#================
# Simulando dados
#================
set.seed(135)
## Parametros
w0 = 0.2
w1 = 90.58
kappa = 1.2 # 1.277
beta = 47.7 # 41.24
# Funcao Micahelis Menten
mich <- function(x,w0=0.2,w1=90,kappa,beta){
  y <-  (w0*beta^kappa + w1*x^kappa)/(beta^kappa + x^kappa)   
  return(y)
}
x <- seq(0,100,1)
y <- mich(x,w0,w1,kappa,beta)
summary(y)
## Gráfico
curve(mich(x,w0,w1,kappa,beta),0,max(x))
# Selecionando dados limitados ou incompletos
df <- data.frame(time = x, peso = mich(x,w0,w1,kappa,beta)+rnorm(n = length(x),sd=0.08*y,)) # estudar outros sigma 0.25
df[which(df$time<=7),]
points(df$time,df$peso, col=2)
points(df$time[which(df$time<=7)],df$peso[which(df$time<=7)], col=2, pch=19)

#----------------------#----------------------#----------------------#----------------------
## Modelando diferentes curvas de crescimento
library(RColorBrewer)
my_colors = brewer.pal(8, "Dark2") 
# Logistico equation
imc <- 18 # 
df_lim <- df[which(df$time<=imc),]
n1 <- nls(formula= peso ~ A/(1+exp(K*(B-time))), data=df_lim,
          start=list(A=60,B=9,K=0.3), trace=TRUE)
summary(n1)
plot(peso~time, data=df, cex=1,cex.axis=1.3,cex.lab=1.5,
     xlab="Weeks after hatching",
     ylab="Shrimp weight (g)",pch=1)
points(df_lim$time,df_lim$peso, col=2, pch=19)
1-deviance(n1)/deviance(lm(peso~1, df_lim))
A=coef(n1)[1]; B=coef(n1)[2]; K=coef(n1)[3]
curve( A/(1+exp(K*(B-x))), add=TRUE, col=my_colors[1], lty=26, lwd=3.5, type="l")
predict(n1,newdata = data.frame(time=c(36,48)))

# Gompertz
n2 <- nls(formula= peso ~ A*exp(-exp(K*(time-B))), data=df_lim,
          start=list(A=70,B=18,K=0.08), trace=TRUE)
summary(n2)
1-deviance(n2)/deviance(lm(peso~1, df_lim))
A=coef(n2)[1]; B=coef(n2)[2]; K=coef(n2)[3]
curve( A*exp(-exp(K*(x-B))), add=TRUE, col=my_colors[7], lty=2,lwd=3.5)
predict(n2,newdata = data.frame(time=c(36,48)))

# Original
curve(mich(x,w0,w1,kappa,beta),0,max(x), add=TRUE, lwd=2.5)

# von Bertalanffy
df_lim <- df[which(df$time<=imc),]
n3 <- nls(formula= peso ~ A*(1-exp(-K*(time-B)))^3, data=df_lim,
          start=list(A=73,B=5,K=0.05), trace=TRUE)
summary(n3)
1-deviance(n3)/deviance(lm(peso~1, df_lim))
A=coef(n3)[1]; B=coef(n3)[2]; K=coef(n3)[3]
curve( A*(1-exp(-K*(x-K)))^3, add=TRUE, col=my_colors[3], lty=2,lwd=3.5)
predict(n3,newdata = data.frame(time=c(36,48)))

# Richards
df_lim <- df[which(df$time<=imc),]
n6 <- nls(formula= peso ~ A*(1+(D-1)*exp(-K*(time-B)))^(1/(1-D)), data=df_lim,
          start=list(A=75,D=0.4,K=0.04,B=12.5), trace=TRUE)
summary(n6)
1-deviance(n6)/deviance(lm(peso~1, df_lim))
A=coef(n6)[1]; D=coef(n6)[2]; K=coef(n6)[3]; B=coef(n6)[4]
curve( A*(1+(D-1)*exp(-K*(x-B)))^(1/(1-D)), add=TRUE, col=my_colors[4], lty=2,lwd=3.5)
predict(n6,newdata = data.frame(time=c(36,48)))
coef(n6)

# Weibull
df_lim <- df[which(df$time<=imc),]
n4 <- nls(formula= peso ~ A*(1-exp(-B*(time^D))), data=df_lim,
          start=list(A=66,B=.02,D=1.07), trace=TRUE)
summary(n4)
1-deviance(n4)/deviance(lm(peso~1, df_lim))
A=coef(n4)[1]; B=coef(n4)[2]; D=coef(n4)[3]
curve( A*(1-exp(-B*(x^D))), add=TRUE, col=my_colors[5], lty=2,lwd=3.5)
predict(n4,newdata = data.frame(time=c(36,48)))
coef(n4)

# Morgan Mercer Flodin
df_lim <- df[which(df$time<=imc),]
n5 <- nls(formula= peso ~ A-((A-B)/(1+(K*time)^D)), data=df_lim,
          start=list(A=78,B=.21,K=0.03,D=1.45), trace=TRUE)
summary(n5)
1-deviance(n5)/deviance(lm(peso~1, df_lim))
A=coef(n5)[1]; B=coef(n5)[2]; K= coef(n5)[3] ; D=coef(n5)[4]
curve( A-((A-B)/(1+(K*x)^D)), add=TRUE, col=my_colors[6], lty=2,lwd=3.5)
predict(n5,newdata = data.frame(time=c(36,48)))
coef(n5)

# Michaelis-Menten equation growth curve
df_lim <- df[which(df$time<=imc),]
n0 <- nls(formula=peso ~ (0.2*B^K + w*time^K)/(B^K + time^K), data=df_lim,
          start=list(w=80,B=40,K=1.2), trace=TRUE)
summary(n0)
w=coef(n0)[1]; B=coef(n0)[2]; K=coef(n0)[3]
curve( (0.2*x^K + w*x^K)/(B^K + x^K)+0.2, add=TRUE, col=my_colors[2], lty=2, lwd=3.5)
coef(n0)

points(df_lim$time,df_lim$peso, col=2, pch=19)

legend("bottom", c("Original curve","Michaelis-Menten","Morgan.M.Flodin","Richard"),  lty=rep(1,4),
       inset=c(0,1.1), xpd=TRUE, horiz=TRUE, bty="n",col=c("black",my_colors[c(2,6,4)]),lwd=4, cex=1.2
)
legend("bottom", c("Weibull","von Bertalanffy","Gompertz","Logistc"),  lty=rep(1,4),
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n",col=my_colors[c(5,3,7,1)],lwd=4, cex=1.2
)

legend(0.5,70,c("Complete data", "Limited data"), pch = c(1,19),
       col=c(1,2),cex=1.2,box.lty=0)

# Salve Simulation.png - plot size 800 x 490




