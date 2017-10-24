#读取数据
genotypeNumber <- read.table('genotypeNumber.dat',header = T)
kind <- rep(0:1,each = 500)
genotypeNumber <- data.frame(kind,genotypeNumber)
genotypeNumber <- lapply(genotypeNumber, as.factor)
genotypeNumber <- as.data.frame(genotypeNumber)

head(genotypeNumber,n = 1)
str(genotypeNumber)
sum(is.na(genotypeNumber))

library(randomForest)
library(e1071)
library(gmodels)

# question2 ---------------------------------------------------------------


#变量选择
genotypeNumber <- read.table('genotypeNumber.dat',header = T)
ill <- rep(0:1,each = 500)
genotypeNumber <- data.frame(ill,genotypeNumber)
genotypeNumber <- lapply(genotypeNumber, as.factor)
genotypeNumber <- as.data.frame(genotypeNumber)

pstd <- 0.05
pvalue <- 0
varsel <- NA
for (j in 2:9445) {
  rs.crosstable <- CrossTable(genotypeNumber$ill,genotypeNumber[,j])
  #rs1.crosstable$t
  #chisq.res <- list()
  for (i in 1:3) {
    tempframe <- data.frame(v1 = c(rs.crosstable$t[1,i],rs.crosstable$t[2,i]), 
                            V2 = 500 - c(rs.crosstable$t[1,i],rs.crosstable$t[2,i]))
    chisq.res <- chisq.test(tempframe)
    pvalue[i] <- chisq.res$p.value
    # print(chisq.test(tempframe))
  }
  if(pvalue[1] <= pstd | pvalue[2] <= pstd | pvalue[3] <= pstd)
    varsel[j] <- 1
  else
    varsel[j] <- 0
}

#选择显著性位点
write.table(colnames(genotypeNumber)[which(varsel == 1)],file = './doc/firstselection.txt') 

#随机森林
set.seed(1)
genotypeNumber.fsel <- genotypeNumber[,c(1,which(varsel == 1))]
model <- randomForest(ill~.,data = genotypeNumber.fsel,importance = T,ntree = 600)
impor <- model$importance
MDG <- data.frame(rs = colnames(genotypeNumber.fsel)[2:871],MeanDecreaseGini = impor[,4])
MDG$rs <- as.character(MDG$rs)
length(which(MDG[,2] >= 0.9))
order(MDG[,2],decreasing = T)[1:13]
write.csv(impor,file = './doc/MeanDecreaseGine3.csv')

#位点MDG
dev.new()
plot(x = 1:nrow(MDG),y = MDG$MeanDecreaseGini,pch = 1,xaxt = 'n',xlab = '位点',ylab = 'MDG',main = '位点MeanDecreaseGini')
axis(side = 1,at = 1:nrow(MDG))
abline(h = 0.9,col = 'red',lwd = 2)
#text(x = which.max(MDG[,2]) + 1,y = MDG['rs2273298',2]-0.1,labels = 'rs2273298' )
text(x = c(order(MDG[,2],decreasing = T)[1:13]),y = c(MDG[order(MDG[,2],decreasing = T)[1:13],2]-0.05),
     labels = c(MDG[order(MDG[,2],decreasing = T)[1:13],1]))
points(x = c(order(MDG[,2],decreasing = T)[1:13]), y = c(MDG[order(MDG[,2],decreasing = T)[1:13],2]),col = 'red')

#选取前13个位点
genotypeNumber.ssel <- genotypeNumber.fsel[,c(1,order(MDG[,2],decreasing = T)[1:13]+1)]

#SVM验证
set.seed(1)
temp1 <- sample(1:500,350)
temp2 <- sample(501:1000,350)
genotypeNumber.ssel.train <- genotypeNumber.ssel[c(temp1,temp2),]
genotypeNumber.ssel.test <- genotypeNumber.ssel[-c(temp1,temp2),]

model <- svm(ill~.,data = genotypeNumber.ssel.train)
#str(model)
modelpred <- predict(model,genotypeNumber.ssel.test)
model.err <- sum(modelpred != genotypeNumber.ssel.test$ill)/nrow(genotypeNumber.ssel.test)

#SVM参数优化
set.seed(1)
best.perf <- 0
for (i in 2:14) {
  if(i != 10 & 1 != 13 & i != 3){
    tune.res <- tune(svm,ill~.,data = genotypeNumber.ssel[,-c(i,10,13,3)],
                 ranges = list(kernel = c('radial','linear','sigmoid','polynomial'),
                               cost = c(0.1,1,5,10)))
tune.res.summ <- summary(tune.res)
best.perf[i] <- tune.res.summ$best.performance
  }
}
order(best.perf)
best.perf

set.seed(1)
tune.res <- tune(svm,ill~.,data = genotypeNumber.ssel[,1:5],
                 ranges = list(kernel = c('radial','linear','sigmoid','polynomial'),
                               cost = c(0.1,1,5,10)))
tune.res.summ <- summary(tune.res)

genotypeNumber.tsel <- genotypeNumber.ssel[,-c(3,10,13)]
write.table(colnames(genotypeNumber.tsel)[2:11],file = './doc/rstop10.txt')

#卡方检验
chisq.res <- 0
for (j in 2:11) {
  rs.crosstable <- CrossTable(genotypeNumber.tsel$ill,genotypeNumber.tsel[,j])
  chisq.res[j-1] <- chisq.test(rs.crosstable$t)$p.value
}
which(chisq.res > 0.01)
write.table(colnames(genotypeNumber.tsel)[-c(1,8)],file = './doc/rsname_2.txt')

#最终SVM模型
set.seed(1)
tune.res <- tune(svm,ill~.,data = genotypeNumber.tsel[,-8],
                 ranges = list(kernel = c('radial','linear','sigmoid','polynomial'),
                               cost = c(0.1,1,5,10)))
tune.res.summ <- summary(tune.res)

set.seed(1)
temp1 <- sample(1:500,350)
temp2 <- sample(501:1000,350)
train <- genotypeNumber.tsel[c(temp1,temp2),-8]
test <- genotypeNumber.tsel[-c(temp1,temp2),-8]
model <- svm(ill~.,data = train,kernel = 'polynomial', cost = 10, probability = T)
pred <- predict(model,test,probability = T)
err <- sum(pred != test$ill)/nrow(test)
#ROC曲线
library(ROCR)
fitted <- attributes(pred)$probabilities[1:300,2]
prediction <- prediction(predictions = fitted,test$ill)
perf <- performance(prediction,measure = 'tpr',x.measure = 'fpr')
dev.new()
plot(perf,main = 'SVM模型ROC曲线', col ='red',lwd = 3)
auc <- performance(prediction,measure = 'auc')@y.values


# question4 ---------------------------------------------------------------


#聚类
multi_phenos <- read.table('multi_phenos_new.txt')
multi_phenos <- as.data.frame(lapply(multi_phenos, as.factor))
str(multi_phenos)
?hclust
?dist
#层次聚类
set.seed(1)
phenos.hclust <- hclust(d = dist(multi_phenos,method = 'euclidean'))
dev.new()
plot(phenos.hclust, hang = 0, labels = F,main = '相关性状层次聚类', xlab = '样本', ylab = '高度')
abline( h = 3.05,col = 'red',lwd = 2)
?hclust
str(phenos.hclust)
#动态聚类
?kmeans
set.seed(1)
phenos.kmeans <- kmeans(multi_phenos,centers = 6)
str(phenos.kmeans)

phenos.kmeans$centers
#保存聚类结果
kind <- phenos.kmeans$cluster
table(kind)

multiphenos.rs <- genotypeNumber.fsel[,-1]
multiphenos.rs <- data.frame(kind,multiphenos.rs)
multiphenos.rs$kind <- as.factor(multiphenos.rs$kind)

#随机森林
set.seed(1)
model <- randomForest(kind~.,data = multiphenos.rs,importance = T,ntree = 600)
err <- mean(model$err.rate)
impor <- model$importance
MDG <- data.frame(rs = colnames(multiphenos.rs)[2:871],MeanDecreaseGini = impor[,8])
MDG$rs <- as.character(MDG$rs)
length(which(MDG[,2] >= 1))
order(MDG[,2],decreasing = T)[1:28]
write.csv(impor,file = './doc/MeanDecreaseGine5.csv')

#位点MDG
dev.new()
plot(x = 1:nrow(MDG),y = MDG$MeanDecreaseGini,pch = 1,xaxt = 'n',xlab = '位点',ylab = 'MDG',main = '位点MeanDecreaseGini')
axis(side = 1,at = 1:nrow(MDG))
abline(h = 1.2,col = 'red',lwd = 2)
#text(x = which.max(MDG[,2]) + 1,y = MDG['rs2273298',2]-0.1,labels = 'rs2273298' )
text(x = c(order(MDG[,2],decreasing = T)[1:28]),y = c(MDG[order(MDG[,2],decreasing = T)[1:28],2]-0.05),
     labels = c(MDG[order(MDG[,2],decreasing = T)[1:28],1]))
points(x = c(order(MDG[,2],decreasing = T)[1:28]), y = c(MDG[order(MDG[,2],decreasing = T)[1:28],2]),col = 'red')

#选取前28个位点
multiphenos.rs.2 <- multiphenos.rs[,c(1,order(MDG[,2],decreasing = T)[1:28]+1)]


#SVM变量选择
set.seed(1)
tune.res <- tune(svm,kind~.,data = multiphenos.rs.2,
                 ranges = list(kernel = c('radial','linear','sigmoid','polynomial'),
                               cost = c(0.1,1,5,10)))
tune.res.summ <- summary(tune.res)


set.seed(1)
best.perf <- 0
fbest.perf <- 0
namesdel <- 0
temp.data <- multiphenos.rs.2
for (j in 1:28) {
  for (i in 1:(29-j)) {
    tune.res <- tune(svm,kind~.,data = temp.data[,-(i+1)],
                 ranges = list(kernel = c('sigmoid'),
                               cost = c(5)))
    #multiphenos.rs.2.temp <- multiphenos.rs.2
    tune.res.summ <- summary(tune.res)
    best.perf[i] <- tune.res.summ$best.performance
}
fbest.perf[j] <- min(best.perf)
numdel <- order(best.perf)[1]+1
namesdel[j] <- colnames(temp.data)[numdel]
temp.data <- temp.data[,-numdel]
}
fbest.perf
namesdel

set.seed(1)
model <- tune(svm,kind~.,data = multiphenos.rs.2[,1:25],
                 ranges = list(kernel = c('radial','linear','sigmoid','polynomial'),
                               cost = c(0.1,1,5,10)))
tune.res.summ <- summary(model)


#卡方检验
chisq.res <- 0
for (j in 2:29) {
  rs.crosstable <- CrossTable(multiphenos.rs.2$kind,multiphenos.rs.2[,j])
  chisq.res[j-1] <- chisq.test(rs.crosstable$t)$p.value
}
which(chisq.res > 0.03)

set.seed(1)
model <- tune(svm,kind~.,data = multiphenos.rs.2[,-c(which(chisq.res > 0.03)+1)],
              ranges = list(kernel = c('radial','linear','sigmoid','polynomial'),
                            cost = c(0.1,1,5,10)))
tune.res.summ <- summary(model)
multiphenos.rs.3 <- multiphenos.rs.2[,-c(which(chisq.res > 0.03)+1)]
colnames(multiphenos.rs.3)[2:12]

chisq.res <- 0
for (j in 2:12) {
  rs.crosstable <- CrossTable(multiphenos.rs.3$kind,multiphenos.rs.3[,j])
  chisq.res[j-1] <- chisq.test(rs.crosstable$t)$p.value
}
chisq.res
which(chisq.res > 0.01)
write.table(colnames(multiphenos.rs.3)[-c(1,which(chisq.res > 0.01)+1)],file = './doc/rsname_4.txt')

