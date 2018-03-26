##------------------------------------------------------------------
rm(list = ls())
library(ggplot2)
library(nlstools)
library(broom)
library(xtable)
library(car)

options(OutDec = ',')
setwd('C:/Users/User/Dropbox/4° Série/Modelos Não Lineares/trabalhos')
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
# Conjunto de dados -------------------------------------------------------

dilution <-rep(c(rep(30,2), rep(90,2), rep(270,2), rep(810,2), rep(2430,2), rep(7290,2), rep(21869,2), rep(65609,2)), 2)
elisa <- data.frame(logd=log10(dilution),
                    OD=c(1.909, 1.956, 1.856, 1.876, 1.838, 1.841, 1.579,
                         1.584, 1.057, 1.072, 0.566, 0.561, 0.225, 0.229,
                         0.072, 0.114, 1.886, 1.880, 1.853, 1.870, 1.747,
                         1.772, 1.424, 1.406, 0.781, 0.759, 0.377, 0.376,
                         0.153, 0.138, 0.053, 0.058),
                    month=c(rep("may", 16), rep("june", 16)))
elisa$month <- factor(elisa$month)
elisa$month <- relevel(elisa$month, 'may')
# Análise descritiva ------------------------------------------------------

ggplot(elisa, aes(x = logd, y = OD, col = month)) +
  geom_point(size = 2) +
  scale_x_continuous(limits = c(1.4, 5), breaks = seq(1.4, 5, l = 5)) +
  scale_color_manual(labels = c("Junho", "Maio"), values = c('royalblue', 'orangered')) +
  labs(x = 'log-diluição', y = 'Densidade óptica', col = '') +
  theme_bw() +
  theme(text = element_text(size = 20), panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(linetype = "dotted")) 
ggsave(filename = 'dispersao.pdf', width = 8, height = 6)

pdf(file = "dispersao2.pdf", width = 8, height = 6)
par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.6)
plot(OD ~ logd, col = c('royalblue', 'orangered')[month], data = elisa, pch = 19, xlim = c(1.4, 5), ylim = c(0, 2),
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex = 0.6)
mtext("log-diluição", side = 1, line = 2.0, cex = 1.8)
mtext("Densidade óptica", side = 2, line = 2, cex = 1.8)
abline(h=seq(0, 2, l = 5), v=seq(1.4, 5, l = 5), col = "gray", lty = "dotted")
axis(1, seq(1.4, 5, l = 5), labels = FF(seq(1.4, 5, l = 5), 1))
axis(2, seq(0, 2, l = 5), labels = FF(seq(0, 2, l = 5), 1))
legend('topright', legend = c('Maio', 'Junho'),  
       col = c('orangered', 'royalblue'), pch = 19, bty = 'n', inset = 0.02, cex = 0.6)
graphics.off()

# Ajuste do modelo --------------------------------------------------------

mod.c <- OD ~ theta1[month] + (theta2[month] - theta1[month]) / (1 + exp(theta3[month] * (logd - theta4[month])))
mod.r <- OD ~ theta1 + (theta2 - theta1) / (1 + exp(theta3 * (logd - theta4)))
ini.c <- list(theta1 = c(0.05, 0.05), theta2 = c(1.9, 1.9), theta3 = c(2.5, 2.5), theta4 = c(3.2, 3.2))
ini.r <- list(theta1 = 0.05, theta2 = 1.9, theta3 = 2.5, theta4 = 3.2)

elisa.c <- nls(formula = mod.c, data = elisa, start = ini.c)
summary(elisa.c)

elisa.r <- nls(formula = mod.r, data = elisa, start = ini.r)
summary(elisa.r)



# Ajuste separado e IC bootstrap ------------------------------------------

elisa.may  <- nls(formula = mod.r, data = subset(elisa, month == 'may'), start = ini.r)
summary(elisa.may)
elisa.june <- nls(formula = mod.r, data = subset(elisa, month == 'june'), start = ini.r)
summary(elisa.june)

boot.may  <- nlsBoot(elisa.may)
boot.june <- nlsBoot(elisa.june)

est <- rbind(tidy(elisa.may), tidy(elisa.june))[, -c(4, 5)]
icb <- rbind(boot.may$bootCI, boot.june$bootCI)[, -1]
fim <- cbind(est, icb)
print(xtable(fim, digits = 4))


# Plots com as curvas ajustadas -------------------------------------------
fx <- function(x, theta)
{
  theta1 <- theta[1]; theta2 <- theta[2]; theta3 <- theta[3]; theta4 <- theta[4]
  theta1 + (theta2 - theta1) / (1 + exp(theta3 * (x - theta4)))
}
x          <- seq(1.4, 5, l = 1000)
theta.may  <- coef(elisa.may)
theta.june <- coef(elisa.june)
fx.may     <- fx(x, theta.may)
fx.june    <- fx(x, theta.june)

pdf(file = "may-june.pdf", width = 8, height = 6)
par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.6)
plot(OD ~ logd, col = c('royalblue', 'orangered')[month], data = elisa, pch = 19, xlim = c(1.4, 5), ylim = c(0, 2),
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex = 0.6)
lines(fx.may ~ x, lwd = 2, col = 'orangered')
lines(fx.june ~ x, lwd = 2, col = 'royalblue')
mtext("log-diluição", side = 1, line = 2.0, cex = 1.8)
mtext("Densidade óptica", side = 2, line = 2, cex = 1.8)
abline(h=seq(0, 2, l = 5), v=seq(1.4, 5, l = 5), col = "gray", lty = "dotted")
axis(1, seq(1.4, 5, l = 5), labels = FF(seq(1.4, 5, l = 5), 1))
axis(2, seq(0, 2, l = 5), labels = FF(seq(0, 2, l = 5), 1))
legend('topright', legend = c('Maio', 'Junho'),  
       col = c('orangered', 'royalblue'), lwd = 1.4, pch = 19, bty = 'n', inset = 0.02, cex = 0.6)
graphics.off()


# Plots das curvas ajustadas com as bandas de confiança -------------------
mod.r1 <- deriv3( ~ theta1 + (theta2 - theta1) / (1 + exp(theta3 * (x - theta4))), 
                  c('theta1', 'theta2', 'theta3', 'theta4'), function(x, theta1, theta2, theta3, theta4){NULL})
df.pred      <- data.frame(logd = seq(1.4, 5, l = 500))
fxnew.may    <- mod.r1(x = df.pred$logd, theta1 = theta.may[1], theta2 = theta.may[2], theta3 = theta.may[3], theta4 = theta.may[4] )
df.pred$ynew <- c(fxnew.may)
F0           <- attr(fxnew.may, 'gradient')
U            <- chol(vcov(elisa.may))  
df.pred$se   <- sqrt(apply(F0 %*% t(U), 1, function(x) sum(x^2)))
df.pred      <- transform(df.pred, lwr=ynew - qnorm(0.975) * se, upr=ynew + qnorm(0.975) * se)

pdf(file = "may.pdf", width = 11, height = 7)
par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.8)
with(elisa, plot(logd, OD, col = 'black', pch = 19, cex = 0.8, xlim = c(1.4, 5), ylim = c(0, 2), xaxt = 'n', 
                   yaxt = 'n', xlab = '', ylab = '', lwd = 2))
with(df.pred, lines(logd, ynew, col = 'royalblue', pch = 19, cex = 0.6))
with(df.pred, lines(logd, lwr, lty = 2, col = 'royalblue'))
with(df.pred, lines(logd, upr, lty = 2, col = 'royalblue'))
mtext("log-diluição", side = 1, line = 2.0, cex = 1.8)
mtext("Densidade óptica", side = 2, line = 2, cex = 1.8)
abline(h=seq(0, 2, l = 5), v=seq(1.4, 5, l = 5), col = "gray", lty = "dotted")
axis(1, seq(1.4, 5, l = 5), labels = FF(seq(1.4, 5, l = 5), 1))
axis(2, seq(0, 2, l = 5), labels = FF(seq(0, 2, l = 5), 1))
title('Maio', cex = 1.5)
graphics.off()


df.pred      <- data.frame(logd = seq(1.4, 5, l = 500))
fxnew.june   <- mod.r1(x = df.pred$logd, theta1 = theta.june[1], theta2 = theta.june[2], theta3 = theta.june[3], theta4 = theta.june[4] )
df.pred$ynew <- c(fxnew.june)
F0           <- attr(fxnew.may, 'gradient')
U            <- chol(vcov(elisa.may))  
df.pred$se   <- sqrt(apply(F0 %*% t(U), 1, function(x) sum(x^2)))
df.pred      <- transform(df.pred, lwr=ynew - qnorm(0.975) * se, upr=ynew + qnorm(0.975) * se)

pdf(file = "june.pdf", width = 11, height = 7)
par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.8)
with(elisa, plot(logd, OD, col = 'black', pch = 19, cex = 0.8, xlim = c(1.4, 5), ylim = c(0, 2), xaxt = 'n', 
                 yaxt = 'n', xlab = '', ylab = '', lwd = 2))
with(df.pred, lines(logd, ynew, col = 'royalblue', pch = 19, cex = 0.6))
with(df.pred, lines(logd, lwr, lty = 2, col = 'royalblue'))
with(df.pred, lines(logd, upr, lty = 2, col = 'royalblue'))
# plotfit(elisa.june, smooth = T, col.fit = 'royalblue', pch.obs = 19, cex = 0.6, xlim = c(1.4, 5), ylim = c(0, 2),
#         xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', lwd = 2)
mtext("log-diluição", side = 1, line = 2.0, cex = 1.8)
mtext("Densidade óptica", side = 2, line = 2, cex = 1.8)
abline(h=seq(0, 2, l = 5), v=seq(1.4, 5, l = 5), col = "gray", lty = "dotted")
axis(1, seq(1.4, 5, l = 5), labels = FF(seq(1.4, 5, l = 5), 1))
axis(2, seq(0, 2, l = 5), labels = FF(seq(0, 2, l = 5), 1))
title('Junho', cex = 1.5)
graphics.off()


# Testando a hipótese (1) -------------------------------------------------
modH0 <- OD ~ theta1 + (theta2  - theta1)/(1 + exp(theta3 * (logd - theta4[month])))
modH1 <- OD ~ theta1[month] + (theta2[month] - theta1[month])/(1 + exp(theta3[month] * (logd - theta4[month])))

fitH0 <- nls(modH0, data = elisa, start = list(theta1 = 0.05, theta2 = 1.9, theta3 = 2.5, theta4 = c(3.2, 3.2)))
fitH1 <- nls(modH1, data = elisa, start = list(theta1 = c(0.05, 0.05), theta2 = c(1.9, 1.9), theta3 = c(2.5, 2.5), theta4 = c(3.2, 3.2)))

anova(fitH0, fitH1)
ll_H0 <- as.numeric(logLik(fitH0))
ll_H1 <- as.numeric(logLik(fitH1))
S_LR  <- 2 * (ll_H1 - ll_H0)
c(S_LR, pchisq(S_LR, df = 3, lower.tail = F))

theta_hat <- diff(coef(fitH1)[-c(7, 8)])[-c(2, 4)]
I_may     <- vcov(fitH1)[c('theta11', 'theta21', 'theta31'), c('theta11', 'theta21', 'theta31')]
I_june    <- vcov(fitH1)[c('theta12', 'theta22', 'theta32'), c('theta12', 'theta22', 'theta32')]
I_hat     <- I_may + I_june
S_W       <- theta_hat %*% solve(I_hat) %*% theta_hat
c(S_W, pchisq(S_W, df = 3, lower.tail = F))

?deltaMethod


# Intervalo de confiança Bootstrap para rho -------------------------------
rho <- function(theta4_may, theta4_june)
{
  10^(-(theta4_may - theta4_june))
}
rho.hat   <- rho(coef(elisa.c)['theta41'], coef(elisa.c)['theta42'])
set.seed(1212)
B                <- 1000
boot.may         <- nlsBoot(nls = elisa.may,  niter = B)
boot.june        <- nlsBoot(nls = elisa.june, niter = B)
theta4.may.boot  <- boot.may$coefboot[, 4]
theta4.june.boot <- boot.june$coefboot[, 4]
rho.boot         <- ldply(lapply(1:B, function(i) rho(theta4_may = theta4.may.boot[i], theta4_june = theta4.june.boot[i])))$V1
qf               <- quantile(rho.boot, c(0.025, 0.975), na.rm = T)


pdf(file = "boot-hist.pdf", width = 8, height = 6)
par(mar = c(3.2, 3.2, 1.5, 1.5), cex = 1.4)
hist(rho.boot, probability = T, xlim = c(0.52, 0.7), xaxt = 'n', yaxt = 'n',
      xlab = '', ylab = '', main = '');box()
axis(side = 1, at = seq(0.52, .7, l = 5), FF(seq(0.52, 0.7, l = 5), 1))
axis(side = 2, at = seq(1, 14, l = 5), FF(seq(1, 14, l = 5), 1))
abline(v = rho.hat, lty = 2, col = 'red')
abline(v = qf[1], lty = 3)
abline(v = qf[2], lty = 3)
graphics.off()







mod.c      <- OD ~ theta1[month] + (theta2[month] - theta1[month]) / (1 + exp(theta3[month] * (logd - theta4[month])))
ini.c      <- list(theta1 = c(0.05, 0.05), theta2 = c(1.9, 1.9), theta3 = c(2.5, 2.5), theta4 = c(3.2, 3.2))
elisa.c    <- nls(formula = mod.c, data = elisa, start = ini.c)
rho.hat    <- rho(coef(elisa.c)['theta41'], coef(elisa.c)['theta42'])
Y.hat      <- fitted(elisa.c)
ei.hat     <- scale(elisa$OD - Y.hat)
elisa$Yhat <- Y.hat
elisa.boot <- elisa
B          <- 1000
rho.boot   <- c()
j          <- 1
while(j < B)
{
  ei.star       <- sample(ei.hat, replace = T)
  elisa.boot$OD <- elisa$Y.hat + ei.star
  fit           <- coef(nls(formula = mod.c, data = elisa.boot, start = ini.c))[c('theta41', 'theta42')]
  rho.boot[j]   <- rho(fit['theta41'], fit['theta42'])
  j <- j + 1
  cat('\n', j)
}


