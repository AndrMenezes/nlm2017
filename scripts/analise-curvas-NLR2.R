rm(list = ls(all.names = TRUE))

setwd('C:/Users/User/Dropbox/4° Série/Modelos Não Lineares/trabalhos/trabalho2')

dados <- read.delim('predicao.txt')
pred  <- dados[dados$dev == 0, ][, -c(2, 3)]
x     <- dados[dados$dev == 1, ]$daa
y     <- dados[dados$dev == 1, ]$peso

FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}

curvas <- function(pNormal, pGumbel, var = 'Constante', j = 1)
{
  pdf(paste0('curvas', j, '.pdf'), width = 9, height = 6)
  par(mar = c(3.2, 3.2, 0.6, 0.6), cex = 1.6)
  plot(y ~ x, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', xlim = c(10, 140), ylim = c(6.2, 304.2), cex = 0.8)
  mtext("Dias após a antese", side = 1, line = 2.0, cex = 1.8)
  mtext("Massa fresca do fruto (g)", side = 2, line = 2, cex = 1.8)
  abline(h=seq(6.2, 304.2, l = 5), v=seq(12, 140, l = 5), col = "gray", lty = "dotted")
  axis(1, seq(12, 140, l = 5))
  axis(2, seq(6.2, 304.2, l = 5), labels = FF(seq(6.2, 304.2, l = 5), 0))
  lines(x = pred$daa, pNormal,  lwd = 2, col = 'royalblue')
  lines(x = pred$daa, pGumbel,  lwd = 2, col = 'red')
  legend(8.8, 310, legend = c('Normal', 'Gumbel'), col = c('royalblue', 'red'), lwd = 2, inset = 0.04, bty = 'n')
  #title(var, cex = 1.5)
  graphics.off()
}

curvas(pNormal = pred$predhm1, pGumbel = pred$predhm2)
curvas(pNormal = pred$predht1, pGumbel = pred$predht2, var = 'Potência', j = 2)
curvas(pNormal = pred$predht3, pGumbel = pred$predht4, var = 'Exponencial', j = 3)



with(pred, lines(x = daa, lowerhm1, lty = 2, col = 'royalblue'))
with(pred, lines(x = daa, upperhm1, lty = 2, col = 'royalblue'))

with(pred, lines(x = daa, lowerhm2, lty = 2, col = 'red'))
with(pred, lines(x = daa, upperhm2, lty = 2, col = 'red'))


plot(y ~ x)
with(pred, lines(x = daa, predht1))
with(pred, lines(x = daa, lowerht1, lty = 'dotted'))
with(pred, lines(x = daa, upperht1, lty = 'dotted'))





dados2     <- dados[, c(1, 2)]
dados2$dev <- 1
nrow(dados2)
ndados     <- cbind(daa = seq(13, 138, l = 1000), peso = 0, dev = 0)
dados3     <- rbind(dados2, ndados)
write.table(x = dados3, file = 'ndados.txt', quote = F, sep = '\t', row.names = F)


