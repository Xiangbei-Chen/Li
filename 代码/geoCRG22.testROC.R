
#install.packages("caret")
#install.packages("DALEX")
#install.packages("ggplot2")
#install.packages("randomForest")
#install.packages("kernlab")
#install.packages("pROC")
#install.packages("xgboost")


#引用包
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

set.seed(123)      #设置种子
inputFile="GSE44770_normalize.txt"         #表达数据文件
geneFile="importanceGene.SVM.txt"      #基因列表文件
setwd("F://WHW//Test_AZ_2022.12.31.6th//25.testROC_GSE44770 #可用")      #设置工作目录

#读取表达数据文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

#读取基因列表文件, 提取疾病特征基因的表达量
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#获取样品分组信息
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

#对数据进行分组
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

#选择模型
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
if(geneFile=="importanceGene.RF.txt"){
	#随机森林树模型
	model=train(Type ~ ., data = train, method='rf', trControl = control)
}else if(geneFile=="importanceGene.SVM.txt"){
	#SVM机器学习模型
	model=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)
}else if(geneFile=="importanceGene.XGB.txt"){
	#XGB模型
	model=train(Type ~., data = train, method = "xgbDART", trControl=control)
}else if(geneFile=="importanceGene.GLM.txt"){
	#GLM模型
	model=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)
}

#绘制ROC曲线
yTest=ifelse(test$Type=="Control", 0, 1)
pred1=predict(model, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file="ROC_SVM.pdf", width=5, height=5)
plot(roc1, print.auc=T, legacy.axes=T, main="", col="red")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()


