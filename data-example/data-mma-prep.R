rm(list=ls())
libraries_check <- c("data.table","mma")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# load dataset
data(weight_behavior)
summary(weight_behavior)
data_raw <- data.table(weight_behavior)
setkey(data_raw)
summary(data_raw)
dim(data_raw)

data_raw <- data_raw[age > 0]
data_raw[, sex := ifelse(sex=="F",yes=1L,no=0L)]
data_raw[race=="", race := NA] # missing values
boxplot(data_raw$exercises) # data-entry error?
data_raw[exercises>100, exercises := NA]
data_raw[, exercises := exercises*1.0]
data_raw[, snack := ifelse(snack==1, 1L, 0L)]
data_raw[, sports := ifelse(sports==1, 1L, 0L)]
data_raw[, sweat := sweat*1.0]
data_raw[, id := 1:nrow(data_raw)]

setnames(data_raw,
         old=c("sex",
               "tvhours","cmpthours","cellhours","exercises","sweat",
               "sports","snack",
               "bmi"),
         new=c("A",paste0("M",1:7),"Y"))
setkey(data_raw)

# remove missing and erroneous observations
data <- data.table(na.omit(data_raw))
dim(data)
nrow(data_raw)-nrow(data)
setcolorder(data,order(names(data)))
setkey(data)
summary(data)
lapply(data,class)
hist(data$Y)

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- r
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex=1.5)
}
# Customize upper panel
upper.panel<-function(x, y){
  x_jitter <- x + runif(n=length(x),min=-.2,max=.2)
  y_jitter <- y + runif(n=length(x),min=-.2,max=.2)
  points(x_jitter,y_jitter, pch = 19, cex=.75)
}
# Create the plots
png("data-mma-pairs-A1.png",width=800,height=600)
pairs(data[A==1,c(paste0("M",1:5),"Y"), with=FALSE], 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="Female")
dev.off()
png("data-mma-pairs-A0.png",width=800,height=600)
pairs(data[A==0,c(paste0("M",1:5),"Y"), with=FALSE], 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="Male")
dev.off()

save(data,file="data-mma.Rdata")