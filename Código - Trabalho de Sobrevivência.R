# TRABALHO - ANÁLISE DE SOBREVIVÊNCIA 

library(survival)
library(tidyverse)
library(AdequacyModel)
library(readr)
library(survminer)


# LEITURA DO BANCO DE DADOS ----

dados <- read_table("lung_cancer.txt")
dados

dados$tratamento <- dados$tratamento %>% as.factor()
dados$celula <- dados$celula %>% as.factor()
dados$terapia_anterior <- dados$terapia_anterior %>% as.factor()
censura = dados$censura %>% as.numeric()
tempo = dados$tempo %>% as.numeric()

# Função de Sobrevivência - Kaplain Meier ----

km <- survfit(Surv(tempo,censura) ~ 1, data = dados) # Estimação de KM 
summary(km)

plot(km, conf.int=F, mark.time=T, lwd=2, col="#082296",
     xlab= "Tempo", ylab="Função de Sobrevivência" )

ggsurvplot(km , data = dados, palette = c("#082296"),
           conf.int = T, size = 1.2, linetype = c("solid"), 
           xlim = c(0,1000), censor = T,    
           xlab = "\nTempo", legend = "none",
           ylab = "Probabilidade de Sobrevivência\n",
           break.time.by = 100,
           censor.shape="+", censor.size = 6,
           censor.color = "red", break.y.by = 0.1,
           ggtheme = theme_minimal(base_size = 14))
           



# Função de Risco Acumulada - Nelson Aalen ----
ena <- survfit(coxph(Surv(tempo,censura) ~ 1,method="breslow", data = dados))
summary(ena)

plot(ena,conf.int=F,fun="cumhaz" )

ena$cumhaz
ggsurvplot(ena , data = dados, fun = "cumhaz", palette = c("#960808"),
           conf.int = F, size = 1.2, linetype = c("solid"), 
           xlim = c(0,1000), censor = F,    
           xlab = "\nTempo", legend = "none",
           ylab = "Função de Risco Acumulada\n",
           break.time.by = 100, ylim = c(0,6),
           ggtheme = theme_minimal(base_size = 14))


# Gráfico TTT
TTT(dados$tempo, col="#960808",lwd=2.5,grid=T,lty=2)




# Análise descritiva - TRATAMENTO ----

km2 <- survfit(Surv(tempo,censura) ~ tratamento, data = dados)
summary(km2)

plot(km2, conf.int=F, mark.time=T, lwd=2, col=c("#082296","#960808"),
     xlab= "Tempo", ylab="Função de Sobrevivência")
legend("topright", bty = "n", legend = c("Padrão","Teste"), fill = c("#082296", "#960808"))

ggsurvplot(km2 , data = dados, palette = c("#082296","#960808"),
           conf.int = F, size = 1.2, linetype = c("solid","solid"), 
           xlim = c(0,1000), censor = F,    
           xlab = "\nTempo", legend = "bottom",
           ylab = "Probabilidade de Sobrevivência\n",
           break.time.by = 100, break.y.by = 0.1,
           ggtheme = theme_minimal(base_size = 14), 
           legend.title = "Tratamento",
           legend.labs = c("Padrão", "Teste"))






# Gráfico TTT
dad2<-data.frame(tempo = dados$tempo,tratamento = dados$tratamento)
tempo1<-dad2[dad2$tratamento == "1",]
TTT(tempo1$tempo, col="#960808", lwd=2.5, grid=TRUE, lty=2)
tempo2<-dad2[dad2$tratamento == "2",]
TTT(tempo2$tempo, col="#960808", lwd=2.5, grid=TRUE, lty=2)


# Teste de Comparação entre curvas

survdiff(Surv(tempo,censura) ~ tratamento, data = dados, rho = 0) # Log rank - Riscos proporcionais
survdiff(Surv(tempo,censura) ~ tratamento, data = dados, rho = 1) # Wilcoxon

#Teste de Schoenfeld - Riscos proporcionais
fit1 <- coxph(Surv(tempo, censura) ~ tratamento, data = dados)
cox.zph(fit1)



# Análise descritiva - CELULA ----

km2 <- survfit(Surv(tempo,censura) ~ celula, data = dados)
summary(km2)

plot(km2, conf.int=F, mark.time=T, lwd=2, col=c("#082296", "black","#960808","#960881"),
     xlab= "Tempo", ylab="Função de Sobrevivência")
legend("topright", bty = "n", legend = c("Escamosa","Pequena","Adeno","Grande"), fill = c("#082296", "#089633","#960808","#960881"))


ggsurvplot(km2 , data = dados, palette = c("#082296", "black","#960808","#960881"),
           conf.int = F, size = 1.2, linetype = c("solid","solid","solid","solid"), 
           xlim = c(0,1000), censor = F,    
           xlab = "\nTempo", legend = "bottom",
           ylab = "Probabilidade de Sobrevivência\n",
           break.time.by = 100, break.y.by = 0.1,
           ggtheme = theme_minimal(base_size = 14), 
           legend.title = "Célula",
           legend.labs = c("Escamosa","Pequena","Adeno","Grande"))




# Gráfico TTT
dad2<-data.frame(tempo = dados$tempo, celula = dados$celula)
tempo1<-dad2[dad2$celula == "1",]
TTT(tempo1$tempo, col="#960808", lwd=2.5, grid=TRUE, lty=2)
tempo2<-dad2[dad2$celula == "2",]
TTT(tempo2$tempo, col="#960808", lwd=2.5, grid=TRUE, lty=2)
tempo3<-dad2[dad2$celula == "3",]
TTT(tempo3$tempo, col="#960808", lwd=2.5, grid=TRUE, lty=2)
tempo4<-dad2[dad2$celula == "4",]
TTT(tempo4$tempo, col="#960808", lwd=2.5, grid=TRUE, lty=2)


# Teste de Comparação entre curvas
survdiff(Surv(tempo,censura) ~ celula, data = dados, rho = 0) # Log rank - Riscos proporcionais
survdiff(Surv(tempo,censura) ~ celula, data = dados, rho = 1) # Wilcoxon

#Teste de Schoenfeld - Riscos proporcionais
fit1 <- coxph(Surv(tempo, censura) ~ celula, data = dados)
cox.zph(fit1)




# Análise descritiva - TERAPIA ANTERIOR ----

km2 <- survfit(Surv(tempo,censura) ~ terapia_anterior, data = dados)
summary(km2)

plot(km2, conf.int=F, mark.time=T, lwd=2, col=c("#960808","#082296"),
     xlab= "Tempo", ylab="Função de Sobrevivência")
legend("topright", bty = "n", legend = c("Não","Sim"), fill = c("#960808","#082296"))

ggsurvplot(km2 , data = dados, palette = c("#960808","#082296"),
           conf.int = F, size = 1.2, linetype = c("solid","solid"), 
           xlim = c(0,1000), censor = F,    
           xlab = "\nTempo", legend = "bottom",
           ylab = "Probabilidade de Sobrevivência\n",
           break.time.by = 100, break.y.by = 0.1,
           ggtheme = theme_minimal(base_size = 14), 
           legend.title = "Terapia Anterior",
           legend.labs = c("Não","Sim"))


# Gráfico TTT
dad2<-data.frame(tempo = dados$tempo,terapia_anterior = dados$terapia_anterior)
tempo1<-dad2[dad2$terapia_anterior == "0",]
TTT(tempo1$tempo, col="#960808", lwd=2.5, grid=TRUE, lty=2)
tempo2<-dad2[dad2$terapia_anterior == "10",]
TTT(tempo2$tempo, col="#960808", lwd=2.5, grid=TRUE, lty=2)


# Teste de Comparação entre curva

survdiff(Surv(tempo,censura) ~ terapia_anterior, data = dados, rho = 0) # Log rank - Riscos proporcionais
survdiff(Surv(tempo,censura) ~ terapia_anterior, data = dados, rho = 1) # Wilcoxon

#Teste de Schoenfeld - Riscos proporcionais
fit1 <- coxph(Surv(tempo, censura) ~ terapia_anterior, data = dados)
cox.zph(fit1)


#Teste de Schoenfeld - Riscos proporcionais - Outras variáveis
fit1 <- coxph(Surv(tempo, censura) ~ score, data = dados)
cox.zph(fit1)

fit1 <- coxph(Surv(tempo, censura) ~ meses_diag, data = dados)
cox.zph(fit1)

fit1 <- coxph(Surv(tempo, censura) ~ idade, data = dados)
cox.zph(fit1)





# Modelo de regressão paramétrico ----


# Distribuição weibull, log normal e log logistica

km <- survfit(Surv(tempo,censura) ~ 1, data = dados)
summary(km)

skm = km$surv

# Distribuição weibull
mwweb = survreg(Surv(tempo,censura) ~ 1, dist = "weibull", data = dados)
summary(mwweb)

alfa_web = exp(mwweb$coefficients)
gama_web = 1 / mwweb$scale

st_weibull = exp(-(km$time/alfa_web)^gama_web)
st_weibull

# Distribuição log normal
mwlognorm = survreg(Surv(tempo,censura) ~ 1, dist = "lognormal", data = dados)
summary(mwlognorm)

mu = mwlognorm$coefficients
sdev = mwlognorm$scale

st_lognorm = 1 - pnorm((log(km$time) - mu) / sdev)
st_lognorm

# Distribuição log logistica
mwloglogist = survreg(Surv(tempo,censura) ~ 1, dist = "loglogistic", data = dados)
summary(mwloglogist)

alfa_loglogist = exp(mwloglogist$coefficients)
gama_loglogist = 1 / mwloglogist$scale

st_mwloglogist = 1 / (1 + (km$time/alfa_loglogist)^gama_loglogist)
st_mwloglogist

# Método 1A

par(mfrow = c(2,2))

plot(km, conf.int=F, mark.time=F, lwd=3, col="black",
     xlab= "Tempo", ylab="Função de Sobrevivência", ylim = c(0,1))
lines(c(0,km$time),c(1,st_weibull), lwd=3, col="#082296")
legend("topright", bty = "n", legend = c("Kaplan Meier","Weibull"), fill = c("black", "#082296"))


plot(km, conf.int=F, mark.time=F, lwd=3, col="black",
     xlab= "Tempo", ylab="Função de Sobrevivência", ylim = c(0,1))
lines(c(0,km$time),c(1,st_lognorm), lwd=3, col="#960808")
legend("topright", bty = "n", legend = c("Kaplan Meier","Log-Normal"), fill = c("black", "#960808"))

plot(km, conf.int=F, mark.time=F, lwd=3, col="black",
     xlab= "Tempo", ylab="Função de Sobrevivência", ylim = c(0,1))
lines(c(0,km$time),c(1,st_mwloglogist), lwd=3, col="#960881")
legend("topright", bty = "n", legend = c("Kaplan Meier","Log-Logístico"), fill = c("black", "#960881"))

plot(km, conf.int=F, mark.time=F, lwd=3, col="black",
     xlab= "Tempo", ylab="Função de Sobrevivência", ylim = c(0,1))
lines(c(0,km$time),c(1,st_weibull), lwd=3, col="#082296")
lines(c(0,km$time),c(1,st_lognorm), lwd=3, col="#960808")
lines(c(0,km$time),c(1,st_mwloglogist), lwd=3, col="#960881")
legend("topright", bty = "n", legend = c("Kaplan Meier","Weibull","Log-Normal","Log-Logístico"), fill = c("black", "#082296", "#960808","#960881"))

par(mfrow = c(1,1))



# AIC - Criério de Akaike
extractAIC(mwweb)
extractAIC(mwlognorm)
extractAIC(mwloglogist)

# AIC - Criério de Akaike Corrigido
MuMIn::AICc(mwweb)
MuMIn::AICc(mwlognorm)
MuMIn::AICc(mwloglogist)


# BIC - Critério de Informação Bayesiana
BIC(mwweb)
BIC(mwlognorm)
BIC(mwloglogist)




# SELEÇÃO DO MODELO - Distribuição Weibull ----

modelo1 <- survreg(Surv(tempo,censura) ~ tratamento, dist = "weibull", data = dados)
summary(modelo1)
anova(modelo1)

modelo2 <- survreg(Surv(tempo,censura) ~ celula, dist = "weibull", data = dados)
summary(modelo2)
anova(modelo2)

modelo3 <- survreg(Surv(tempo,censura) ~ score, dist = "weibull", data = dados)
summary(modelo3)
anova(modelo3)

modelo4 <- survreg(Surv(tempo,censura) ~ meses_diag, dist = "weibull", data = dados)
summary(modelo4)
anova(modelo4)

modelo5 <- survreg(Surv(tempo,censura) ~ idade, dist = "weibull", data = dados)
summary(modelo5)
anova(modelo5)

modelo6 <- survreg(Surv(tempo,censura) ~ terapia_anterior, dist = "weibull", data = dados)
summary(modelo6)
anova(modelo6)

modelo7 <- survreg(Surv(tempo,censura) ~ celula + score, dist = "weibull", data = dados)
summary(modelo7)
anova(modelo7)
anova(modelo2, modelo7)
anova(modelo3, modelo7)

modelo8 <- survreg(Surv(tempo,censura) ~ celula + score + tratamento, dist = "weibull", data = dados)
summary(modelo8)
anova(modelo8)
anova(modelo7, modelo8)

modelo9 <- survreg(Surv(tempo,censura) ~ celula + score + meses_diag, dist = "weibull", data = dados)
summary(modelo9)
anova(modelo9)
anova(modelo7, modelo9)


modelo10 <- survreg(Surv(tempo,censura) ~ celula + score + idade, dist = "weibull", data = dados)
summary(modelo10)
anova(modelo10)
anova(modelo7, modelo10)


modelo11 <- survreg(Surv(tempo,censura) ~ celula + score + terapia_anterior, dist = "weibull", data = dados)
summary(modelo11)
anova(modelo11)
anova(modelo7, modelo11)

modelo12 <- survreg(Surv(tempo,censura) ~ ., dist = "weibull", data = dados)
summary(modelo12)
anova(modelo12)
anova(modelo7, modelo12)

modelo13 <- survreg(Surv(tempo,censura) ~ celula + score + terapia_anterior + idade + tratamento, dist = "weibull", data = dados)
summary(modelo13)
anova(modelo13)
anova(modelo7, modelo13)

modelo14 <- survreg(Surv(tempo,censura) ~ celula + score + idade + tratamento, dist = "weibull", data = dados)
summary(modelo14)
anova(modelo14)
anova(modelo7, modelo14)

modelo15 <- survreg(Surv(tempo,censura) ~ celula + score + tratamento + tratamento:celula, dist = "weibull", data = dados)
summary(modelo15)
anova(modelo15)
anova(modelo7, modelo15)

modelo16 <- survreg(Surv(tempo,censura) ~ celula + score + celula:score, dist = "weibull", data = dados)
summary(modelo16)
anova(modelo16)
anova(modelo7, modelo16)

modelo17 <- survreg(Surv(tempo,censura) ~ celula + score + celula:score + tratamento + tratamento:celula, dist = "weibull", data = dados)
summary(modelo17)
anova(modelo17)
anova(modelo7, modelo17)

# AIC - Criério de Akaike
extractAIC(modelo7)
extractAIC(modelo15)


# AIC - Criério de Akaike Corrigido
MuMIn::AICc(modelo7)
MuMIn::AICc(modelo15)

# BIC - Critério de Informação Bayesiana
BIC(modelo7)
BIC(modelo15)



# ANÁLISE DE RESÍDUOS ----

#Residuo Martingal
martingal = resid(modelo7,type="martingale")

#Resíduo Cox-Snell
coxsnell <- censura - resid(modelo7,type="martingale")
res <- abs(coxsnell)

#Residuo Deviance
deviance = resid(modelo7,type="deviance")


#Resíduo Cox-Snell
summary(modelo7)

coxsnell = (tempo * exp(- modelo15$linear.predictors))^(1/modelo15$scale)
coxsnell

#Residuo Martingal
martingal = censura - coxsnell
martingal

#Residuo Deviance
deviance = resid(modelo7,type="deviance")
deviance


ekm = survfit(Surv(sort(coxsnell),censura) ~ 1, type = c("kaplan-meier"))
summary(ekm) 

st = ekm$surv
t = ekm$time
sexp1 = exp(-t)

par(mfrow = c(1,2))
plot(st,sexp1, xlab = "S(e): Kaplan-Meier", pch = 16,
     ylab = "S(e): Exponecial Padrão", col="#082296")
lines(st, st, lty = 1, col = "#960808",lwd=3)

plot(ekm, conf.int=F, lwd=3, col="#082296",
     xlab= "Resíduos de Cox-Snell", ylab="S(e) estimada", lty = c(1,1))
lines(t, sexp1, lty = 1, col = "#960808",lwd=3)
legend(3,1,lty = c(1,1), bty = "n", cex = 0.9, lwd=2,
       c("Kaplan-Meier","Exponencial Padrão"), col = c("#082296","#960808"))

sort(coxsnell)
ekm$cumhaz

par(mfrow = c(1,1))
plot(sort(coxsnell), c(ekm$cumhaz), xlim = c(0,8), 
     xlab = "Resíduos de Cox-Snell", ylab = "Função Acumulada do RCS")
lines(sort(coxsnell), sort(coxsnell), lty = 1, col = "#960808",lwd=2)

plot(abs(martingal))
plot(deviance)

plot(martingal, tempo)
plot(abs(martingal), log(tempo))

plot(deviance, tempo)
plot(deviance, log(tempo))





