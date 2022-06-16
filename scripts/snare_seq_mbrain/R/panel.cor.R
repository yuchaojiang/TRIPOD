# For pairs() function to plot pairwise plot
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor=NULL, ...)
{ usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
filter=!is.na(y) & !is.na(x) & !is.nan(x) & !is.nan(y) & !is.infinite(x) & !is.infinite(y)
r <- cor(x[filter], y[filter])
txt <- format(c(r, 0.123456789), digits = digits)[1]
txt <- paste0(prefix, 'r = ', txt)
cex.cor <- 0.5/strwidth(txt)
text(0.5, 0.5, txt, cex = cex.cor * r*2)}
