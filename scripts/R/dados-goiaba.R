# Definições gerais -------------------------------------------------------
rm(list = ls(all.names = TRUE))
# setwd('/home/andrefelipe/Dropbox/4° Série/Modelos Não Lineares/trabalhos/trabalho2')
setwd('C:/Users/User/Dropbox/4° Série/Modelos Não Lineares/trabalhos/trabalho2')
# options(OutDec = ',')

bib <- c('nlstools', 'nlme', 'xtable', 'car', 'fitdistrplus')
sapply(bib, require, character.only = T)
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}

# Leitura dos dados -------------------------------------------------------
dados <- read.table("http://www.leg.ufpr.br/~walmes/data/goiaba.txt", header=TRUE, sep="\t")
# write.table(x = dados, file = 'goiaba.txt', sep = '\t', quote = F, row.names = F)

# Análise descritiva ------------------------------------------------------
pdf('hist.pdf', width = 9, height = 6)
par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.6)
hist(dados$peso, col = 'gray', main = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = ''); box()
axis(1, seq(0, 350, l = 6)); axis(2, seq(0, 80, l = 5))
mtext('Massa fresca do fruto (g)', 1, line =2, cex = 1.8)
mtext('Frequência', 2, line =2, cex = 1.8)
graphics.off()
x <- seq(min(dados$peso), max(dados$peso))
lines(x, dnorm(x, mean(dados$peso), sd(dados$peso)))
fit <- fitdist(dados$peso, 'gumbel', start = c(1, 1))
lines(x, dinvgauss(x, fit$estimate[1], fit$estimate[2]))

pdf(file = "dispersao.pdf", width = 8, height = 6)
par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.6)
plot(peso ~ daa, data = dados, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', cex = 0.8,
     xlim = c(10, 140), ylim = c(6.2, 304.2))
lines(smooth.spline(dados$peso ~ dados$daa), lwd = 2, col = 'blue')
mtext("Dias após a antese", side = 1, line = 2.0, cex = 1.8)
mtext("Massa fresca do fruto (g)", side = 2, line = 2, cex = 1.8)
abline(h=seq(6.2, 304.2, l = 5), v=seq(12, 140, l = 5), col = "gray", lty = "dotted")
axis(1, seq(12, 140, l = 5))
axis(2, seq(6.2, 304.2, l = 5), labels = FF(seq(6.2, 304.2, l = 5), 1))
graphics.off()

pdf(file = "boxplot.pdf", width = 10.5, height = 6.5)
par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.6)
boxplot(peso ~ daa, data = dados,  xlab = '', ylab = '', cex = 0.5, col = 'gray')
mtext("Dias após a antese", side = 1, line = 2.0, cex = 1.8)
mtext("Massa fresca do fruto (g)", side = 2, line = 2, cex = 1.8)
graphics.off()

media   <- with(dados, tapply(peso, daa, mean, na.rm = T))
mediana <- with(dados, tapply(peso, daa, median, na.rm = T))
dp      <- with(dados, tapply(peso, daa, sd, na.rm = T))
n       <- with(dados, tapply(peso, daa, length))
xtable(FF((cbind(media, mediana, dp)), 2))


pdf(file = "desvio-dias.pdf", width = 8, height = 6)
par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.6)
Rx <- range(dp)
plot(dp ~ unique(dados$daa), xlim = c(10, 140), xlab = '', ylab = '', xaxt = 'n', 
     yaxt = 'n')
mtext("Dias após a antese", side = 1, line = 2.0, cex = 1.8)
mtext("Desvio padrão da massa fresca", side = 2, line = 2, cex = 1.8)
abline(h=seq(Rx[1], Rx[2], l = 5), v=seq(12, 140, l = 5), col = "gray", lty = "dotted")
axis(1, seq(12, 140, l = 5))
axis(2, seq(Rx[1], Rx[2], l = 5), labels = FF(seq(Rx[1], Rx[2], l = 5), 1))
graphics.off()



bartlett.test(peso ~ factor(daa), data = dados)
leveneTest(peso ~ factor(daa), data = dados)

# Modelos -----------------------------------------------------------------
gompertz <- peso ~ alpha * exp(-exp(beta - gamma * daa))
logistic <- peso ~ alpha * (1 + exp(beta - gamma * daa))^(-1)

# Chutes iniciais ---------------------------------------------------------
alpha0     <- 180 # assintóta quando x-> infinity
c.gompertz <- abs(coef(lm(log(-log(peso / alpha0)) ~ daa, data = dados)))
c.logistic <- abs(coef(lm(logit(peso / alpha0)  ~ daa, data = dados)))

# Ajuste dos modelos homocedásticos ---------------------------------------
chute1   <- list(alpha = 138, beta = 4, gamma = 0.09)

homo.gompertz <- nls(formula = gompertz, data = dados, start = chute1)
homo.logistic <- nls(formula = logistic, data = dados, start = chute1)


# Ajuste dos modelos heterocedásticos -------------------------------------
hete.gompertz <- gnls(gompertz, data = dados, start = chute1, weights = varExp(0.1, form = ~ daa))
hete.logistic <- gnls(logistic, data = dados, start = chute1, weights = varPower(0.1, form = ~ daa))
chute1   <- list(alpha = 138, beta = 4, gamma = 0.09)
hete.logistic <- gnls(logistic, data = dados, start = c(alpha = 232, beta = 7.49, gamma = 0.07378),
                      weights = varExp(0.02, form = ~ daa))
summary(hete.logistic)

?varPower()



m0 <- nls(peso~A-(A-D)*exp(-exp(C*(daa-B))), data=dados,
          start=c(A=200, C=0.09, B=105, D=7))
summary(m0)


A-(A-D)*exp(-exp(C*(daa-B)))




Chick.6 <- subset(ChickWeight, (Chick == 6) & (Time > 0))
SSweibull(Chick.6$Time, 160, 115, -5.5, 2.5)   # response only
Asym <- 160; Drop <- 115; lrc <- -5.5; pwr <- 2.5
SSweibull(Chick.6$Time, Asym, Drop, lrc, pwr)  # response and gradient
getInitial(weight ~ SSweibull(Time, Asym, Drop, lrc, pwr), data = Chick.6)
getInitial(peso ~ SSweibull(daa, Asym, Drop, lrc, pwr), data = dados)


## Initial values are in fact the converged values
fm1 <- nls(weight ~ SSweibull(Time, Asym, Drop, lrc, pwr), data = Chick.6)
summary(fm1)




rmse <- function(yhat, n, p) sum((preditos$peso - yhat)^2) / (n - p)

preditos <- read.delim('pred.txt')
head(preditos)
n <- nrow(preditos); p <- 6
(rmse1 <- rmse(yhat = preditos$predhm1, n = n, p = 5))
(rmse2 <- rmse(yhat = preditos$predhm2, n = n, p = 5))

(rmse3 <- rmse(yhat = preditos$predht1, n = n, p = 6))
(rmse4 <- rmse(yhat = preditos$predht2, n = n, p = 6))

(rmse5 <- rmse(yhat = preditos$predht3, n = n, p = 6))
(rmse6 <- rmse(yhat = preditos$predht4, n = n, p = 6))










require(gnlm)

f <- c(215, 1485, 5331, 10649, 14959, 11929, 6678, 2092, 342)
y <- seq(0,8)
fit.dist(y, f, "binomial", plot=TRUE, xlab="Number of males",
         main="Distribution of males in families of 8 children")
#
f <- c(1,1,6,3,4,3,9,6,5,16,4,11,6,11,3,4,5,6,4,4,5,1,1,4,1,2,
       0,2,0,0,1)
y <- seq(1100,4100,by=100)
fit.dist(y, f, "normal", delta=100, plot=TRUE,
         xlab="Monthly salary (dollars)",
         main="Distribution of women mathematicians
         '
         salaries")
fit.dist(y, f, "log normal", delta=100, plot=TRUE, add=TRUE, lty=3)
fit.dist(y, f, "logistic", delta=100, exact=FALSE, plot=TRUE, add=TRUE, lty=2)


sex <- c(rep(0,10),rep(1,10))
sexf <- gl(2,10)
age <- c(8,10,12,12,8,7,16,7,9,11,8,9,14,12,12,11,7,7,7,12)
y <- cbind(c(9.2, 7.3,13.0, 6.9, 3.9,14.9,17.8, 4.8, 6.4, 3.3,17.2,
             14.4,17.0, 5.0,17.3, 3.8,19.4, 5.0, 2.0,19.0),
           c(0,1,0,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1))


# or equivalently
gnlr(y, dist="inverse Gauss", mu=~exp(b0+b1*sex+b2*age),
     pmu=list(b0=3,b1=0,b2=0), pshape=-1)











head(dados)
sum(is.na(dados))
dados
library(gamlss.nl)
data(la)
mod1 <- nlgamlss(y = PET60, mu.fo = ~bflow * (1 - (1 - exp(p1)) * exp(-p2/bflow)), 
                 sigma.formula = ~1, mu.start = c(-0.9, 90), sigma.start = log(0.1),
                 nu.start = 0, tau.start = log(2.5), family = BCPE, data = la)
mod1 <- nlgamlss(y = peso, mu.formula = ~ theta1 - (theta1 - theta2) * exp(-exp(theta3 * (daa - theta4))),
                 sigma.formula = ~1, mu.start = c(190, 18, 0.09, 107), sigma.start = 1, family = RG, 
                 data = na.omit(dados))

?dGU()
dBCPE()


fxtheta <- 5*x
lambda  <- 1.5
sigma   <- 2.0
sigmai  <- sigma * x^lambda
ei <- rgumbel(174, 0, sigmai)
yobs  <- fxtheta + ei
yhat  <- fxtheta + (-digamma(1)) * sigmai

ei <- (yobs - yhat) / sigmai
mean(ei)


par <- coef(fitdist(y, 'gumbel', start = list(location = 1, scale = 1)))
Ey  <- par[1] + (-digamma(1)) * par[2]
z1  <- (y - Ey) / par[2] 
coef(fitdist(z1, 'gumbel', start = list(location = 1, scale = 1)))
z2  <- (y - par[1]) / par[2] 
coef(fitdist(z2, 'gumbel', start = list(location = 1, scale = 1)))























