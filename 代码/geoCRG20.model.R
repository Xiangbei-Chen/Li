
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
inputFile="normalize.txt"      #表达数据文件
geneFile="interGenes.txt"      #基因列表文件
setwd("C://Users//41939//Desktop//大创项目文件夹//子宫内膜异位症//Reanalysis_GSE51981//22.model")      #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#读取基因列表文件,提取交集核心基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
row.names(data)=gsub("-", "_", row.names(data))

#获取样品分组信息
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

#对数据进行分组
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

#RF随机森林树模型
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
mod_rf = train(Type ~ ., data = train, method='rf', trControl = control)

#SVM机器学习模型
mod_svm=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)

#XGB模型
mod_xgb=train(Type ~., data = train, method = "xgbDART", trControl=control)

#GLM模型
mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)


#定义预测函数
p_fun=function(object, newdata){
	predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(test$Type=="Control", 0, 1)

#RF随机森林树模型预测结果
explainer_rf=explain(mod_rf, label = "RF",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_rf=model_performance(explainer_rf)
#SVM机器学习模型预测结果
explainer_svm=explain(mod_svm, label = "SVM",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_svm=model_performance(explainer_svm)
#XGB模型预测结果
explainer_xgb=explain(mod_xgb, label = "XGB",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_xgb=model_performance(explainer_xgb)
#GLM模型预测结果
explainer_glm=explain(mod_glm, label = "GLM",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_glm=model_performance(explainer_glm)

#绘制四种方法的残差反向累计分布图
pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm)
print(p1)
dev.off()

#绘制四种方法的残差箱线图
pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()


#绘制ROC曲线
pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
pred3=predict(mod_xgb, newdata=test, type="prob")
pred4=predict(mod_glm, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="red")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="blue", add=T)
plot(roc3, print.auc=F, legacy.axes=T, main="", col="green", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="yellow", add=T)
legend('bottomright',
	   c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
	     paste0('SVM: ',sprintf("%.03f",roc2$auc)),
	     paste0('XGB: ',sprintf("%.03f",roc3$auc)),
	     paste0('GLM: ',sprintf("%.03f",roc4$auc))),
	   col=c("red","blue","green","yellow"), lwd=2, bty = 'n')
dev.off()

#对四种方法进行基因的重要性分析,得到四种方法基因重要性评分
importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_xgb<-variable_importance(
  explainer_xgb,
  loss_function = loss_root_mean_square
)
#绘制基因重要性图形
pdf(file="importance.pdf", width=7, height=10)
plot(importance_rf[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_svm[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_xgb[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_glm[c(1,(ncol(data)-8):(ncol(data)+1)),])
dev.off()
#输出重要性评分最高的基因
geneNum=5     #设置基因的数目
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_xgb[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.XGB.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)


