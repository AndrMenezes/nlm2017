library(VGAM); library(fitdistrplus); library(goftest)
rm(list = ls(all.names = TRUE))

setwd('C:/Users/User/Dropbox/4° Série/Modelos Não Lineares/trabalhos/trabalho2')

preditos <- read.delim('ypred.txt')
preditos2 <- read.delim('ypred2.txt')
tail(preditos)
tail(preditos2)


# Funções -----------------------------------------------------------------
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}

gof <- function(x, p, ...)
{
  xks       <- ks.test(x  = x, p, ...)
  xcvm      <- cvm.test(x = x, null = p, ...)
  xad       <- ad.test(x  = x, null = p, ...)
  statistic <- unname(c(xks$statistic, xcvm$statistic, xad$statistic))
  pvalue    <- c(xks$p.value, xcvm$p.value, xad$p.value)
  paste0(FF(statistic, 3), " (", FF(pvalue, 3), ")")
}

residuos <- function(yhat, sigma, var = 'cte', lambda = 0, dist = 'normal')
{
  x <- preditos$daa;  yobs <- preditos$peso
  if(var == 'pot') sigma_i <- sigma * x^lambda
  if(var == 'exp') sigma_i <- sigma * exp(lambda * x)
  if(var == 'cte') sigma_i <- sigma
  
  res  <- (yobs - yhat) / sigma_i
  qemp <- sort(res)
  if(dist == 'normal') 
  {
    qteo <- qnorm(ppoints(length(yhat))) 
    tt   <- gof(qemp, p = 'pnorm')
  }
  if(dist == 'gumbel')
  {
    qteo <- qgumbel(ppoints(length(yhat)), location = digamma(1))
    tt   <- gof(qemp, p = 'pgumbel', location = digamma(1))
  }
  
  Rx <- range(qteo); seqRx <- seq(Rx[1], Rx[2], l = 5); seqRx <- ifelse(seqRx >= -0.009 & seqRx < 0, 0, seqRx)
  Ry <- range(qemp); seqRy <- seq(Ry[1], Ry[2], l = 5); seqRy <- ifelse(seqRy >= -0.009 & seqRy < 0, 0, seqRy)
  pdf(paste0('qqplot-', dist, '-', var, '.pdf'), width = 9, height = 6)
  par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.6)
  plot(qemp ~ qteo, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', xlim = Rx, ylim = Ry, cex = 0.6)
  abline(h=seqRy, v=seqRx, col = "gray", lty = "dotted"); abline(a = 0, b = 1)
  axis(1, seqRx, FF(seqRx, 2)); axis(2, seqRy, FF(seqRy, 2))
  mtext("Quantis teóricos", side = 1, line = 2.0, cex = 1.8)
  mtext("Resíduos", side = 2, line = 2, cex = 1.8)  
  if(dist == 'normal') title('Normal', cex = 1.5)
  if(dist == 'gumbel') title('Gumbel', cex = 1.5)
  graphics.off()  
  
  Rx <- c(10, 140); seqRx <- seq(Rx[1], Rx[2], l = 6)
  Ry <- range(res); seqRy <- seq(Ry[1], Ry[2], l = 5); seqRy <- ifelse(seqRy >= -0.009 & seqRy < 0, 0, seqRy)
  pdf(paste0('res-', dist, '-', var, '.pdf'), width = 9, height = 6)
  par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.6)
  plot(res ~ x, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', xlim = Rx, ylim = Ry, cex = 0.8)
  abline(h=seqRy, v=seqRx, col = "gray", lty = "dotted")
  axis(1, seqRx, FF(seqRx, 0)); axis(2, seqRy, FF(seqRy, 2))
  mtext("Dias após a antese", side = 1, line = 2.0, cex = 1.8)
  mtext("Resíduos", side = 2, line = 2, cex = 1.8)  
  if(dist == 'normal') title('Normal', cex = 1.5)
  if(dist == 'gumbel') title('Gumbel', cex = 1.5)
  graphics.off()  

  return(tt)
}
mNHo  <- residuos(yhat = preditos$predhm1, sigma = 24.6528, var = 'cte', dist = 'normal')
mGHo  <- residuos(yhat = preditos$predhm2, sigma = 20.9402, var = 'cte', dist = 'gumbel')

mNHtp <- residuos(yhat = preditos$predht1, sigma = 0.0132, var = 'pot', lambda = 1.6371, dist = 'normal')
mGHtp <- residuos(yhat = preditos$predht1, sigma = 0.0146, var = 'pot', lambda = 1.5878, dist = 'gumbel')

mNHte <- residuos(yhat = preditos$predht3, sigma = 1.0583, var = 'exp', lambda = 0.03011, dist = 'normal')
mGHte <- residuos(yhat = preditos$predht3, sigma = 0.9000, var = 'exp', lambda = 0.0305, dist = 'gumbel')

TT <- cbind(mNHo, mGHo, mNHtp, mGHtp, mNHte, mGHte)
xtable::xtable(TT)










mNHo  <- residuos(yobs = preditos$peso, yhat = preditos$predhm1, sigma = 24.6528, var = 'cte', dist = 'normal')
ks.test(mNHo, 'pnorm')
fitdist(mNHo, 'norm')

mGHo  <- residuos(yobs = preditos$peso, yhat = preditos$predhm2, sigma = 20.9402, var = 'cte', dist = 'gumbel')
ks.test(mGHo, 'pgumbel', location = digamma(1))
fitdist(mGHo, 'gumbel', start = list(location = 1, scale = 1))

mNHtp <- residuos(yobs = preditos$peso, yhat = preditos$predht1, sigma = 0.0132, var = 'pot',
                  x = preditos$daa, lambda = 1.6371, dist = 'normal')
ks.test(mNHtp, 'pnorm')

mGHtp <- residuos(yobs = preditos$peso, yhat = preditos$predht1, sigma = 0.0146, var = 'pot',
                  x = preditos$daa, lambda = 1.5878, dist = 'gumbel')
ks.test(mGHtp, 'pgumbel', location = digamma(1))
fit = fitdist(mGHtp, 'gumbel', start = list(location = 1, scale = 1))

mNHte <- residuos(yobs = preditos$peso, yhat = preditos$predht3, sigma = 1.0583, var = 'exp',
                  x = preditos$daa, lambda = 0.03011, dist = 'normal')
ks.test(mNHte, 'pnorm')

mGHte <- residuos(yobs = preditos$peso, yhat = preditos$predht3, sigma = 0.9000, var = 'exp',
                  x = preditos$daa, lambda = 0.0305, dist = 'gumbel')
ks.test(mGHte, 'pgumbel', location = digamma(1))
fitdist(mGHte, 'gumbel', start = list(location = 1, scale = 1))



