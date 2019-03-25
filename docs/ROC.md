## ROC 曲线
在信号检测理论中，接收者操作特征曲线（receiver operating characteristic curve）是一种坐标图式的分析工具，用于（1）选择最佳的信号侦测模型，舍弃次佳的模型
（2）在同一模型中设定最佳阈值。

ROC分析的是二元分类模型，也就是输出结果只有两种类别的模型，例如：（阳性/阴性）（有病/没病）等。

**术语**

| * | * | 真实值 | 真实值 | 总数 |
| :--- | :--- | :--- | :--- | :--- | 
| * | * | p(positive) | n(negative) | * |
| 预测值 | p' | 真阳性（TP） | 伪阳性（FP） | P' |
| 预测值 | n' | 伪阴性（FN） | 真阴性（TN） | N' |
| * | 总数 | P | N | 

真阳性（TP）:正确的肯定，诊断为有，实际上也有高血压。

伪阳性（FP）:错误的肯定，诊断为有，实际上却没有高血压。第一型错误。

真阴性（TN）:正确的否定，诊断为没有，实际上也没有高血压。

伪阴性（FN）:错误的否定，诊断为没有，实际上却有高血压。第二型错误。

真阳性率（TPR, true positive rate）,又称为 命中率（hit rate），敏感度（sensitivity）

TPR = TP / P = TP / (TP+FN)

伪阳性率（FPR, false positive rate）,又称为 命中率，假警报率 (false alarm rate)

FPR = FP / N = FP / (FP + TN)

真阴性率（TNR）, 又称为 特异度（SPC, specificity）

SPC = TN / N = TN / (FP + TN) = 1 - FPR

准确度 (ACC, accuracy)

ACC = (TP + TN) / (P + N), 即：(真阳性+真阴性) / 总样本数

阳性预测值 (PPV)：PPV = TP / (TP + FP)

阴性预测值 (NPV)：NPV = TN / (TN + FN)

假发现率 (FDR)：FDR = FP / (FP + TP)






