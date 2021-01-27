#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(error=expression(NULL))
##Scale control urine
urine_cont <- read.csv(file = args[1], header=TRUE, sep="\t")
x <- urine_cont$Count
y <- urine_cont$control_urine
df <- data.frame(x=x, y=y)
model <- lm(y ~ x, data = df)
expec <- coef(model)[[2]]*x + coef(model)[[1]]
new_tab <- cbind(x,y,expec)
ste <- summary(model)[["coefficients"]][1,2]
dev <- (y-expec)/ste
taxon <- urine_cont$Taxa
results <- cbind.data.frame(taxon, new_tab, dev)
write.table(results, args[3], col.names=NA)
urine_usd <- read.csv(file = args[2], header=TRUE, sep="\t")
x <- urine_usd$Count
y <- urine_usd$USD_urine
df <- data.frame(x=x, y=y)
model <- lm(y ~ x, data = df)
expec <- coef(model)[[2]]*x + coef(model)[[1]]
new_tab <- cbind(x,y,expec)
ste <- summary(model)[["coefficients"]][1,2]
dev <- (y-expec)/ste
taxon <- urine_usd$Taxa
results <- cbind.data.frame(taxon, new_tab, dev)
write.table(results, args[4], col.names=NA)
stool_cont <- read.csv(file = args[5], header=TRUE, sep="\t")
x <- stool_cont$Count
y <- stool_cont$control_stool
df <- data.frame(x=x, y=y)
model <- lm(y ~ x, data = df)
expec <- coef(model)[[2]]*x + coef(model)[[1]]
new_tab <- cbind(x,y,expec)
ste <- summary(model)[["coefficients"]][1,2]
dev <- (y-expec)/ste
taxon <- stool_cont$Taxa
results <- cbind.data.frame(taxon, new_tab, dev)
write.table(results, args[6], col.names=NA)
stool_usd <- read.csv(file = args[7], header=TRUE, sep="\t")
x <- stool_usd$Count
y <- stool_usd$USD_stool
df <- data.frame(x=x, y=y)
model <- lm(y ~ x, data = df)
expec <- coef(model)[[2]]*x + coef(model)[[1]]
new_tab <- cbind(x,y,expec)
ste <- summary(model)[["coefficients"]][1,2]
dev <- (y-expec)/ste
taxon <- stool_usd$Taxa
results <- cbind.data.frame(taxon, new_tab, dev)
write.table(results, args[8], col.names=NA)
