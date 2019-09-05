setwd("C:\\workplace\\文章\\iRPBGO\\整理的数据程序\\")
##################################
#计算表达差异EcoMAC_matrixg由170302EcoIN基因二级B号分离.py处理后手工加上标题获得
#此矩阵中有基因ID冗余，需要进一步处理
##################################
data_ori<-read.table('EcoMAC_matrixg.txt',header = T,sep = '\t')
#处理ID冗余问题
table_data<-table(data_ori[,1])
redu_gene_list<-names(table_data[table_data>1])
redu_gene_location<-c()
redu_gene_matrix<-as.data.frame(matrix(nrow =length(redu_gene_list),ncol=dim(data_ori)[2] ))
colnames(redu_gene_matrix)<-colnames(data_ori)

for (redui in 1:length(redu_gene_list))
{
    gene_location<-which(data_ori[,1]==redu_gene_list[redui])
    redu_gene_location<-c(redu_gene_location,gene_location)
    genedata<-data_ori[gene_location,]
    genedata_mean<-apply(genedata[,2:dim(data_ori)[2]],2,function(x){
        mean(x) 
    })  
    redu_gene_matrix[redui,1]=redu_gene_list[redui]
    redu_gene_matrix[redui,2:dim(data_ori)[2]]=genedata_mean
    
}
unidata<-data_ori[-redu_gene_location,]

data<-rbind(unidata,redu_gene_matrix)

#设置对照试验
control_name<-read.table('control_data_M9_e.txt',sep = '\t')
control_location<-c()
for (i in 1:length(control_name[,1]))
{
    control_location<-c(control_location,which(colnames(data)==control_name[i,1]) ) 
}
case_data<-data[,-c(1,control_location)]
control_data<-data[,control_location]
row.names(case_data)<-data[,1]
use_control_data<-matrix(apply(control_data,1,mean),nrow = dim(data)[1],dimnames = list(data[,1],'WT_M9'))
##################################
#计算表达差异
##################################
diff_expr_data<-apply(case_data,2,function(x){
    x-use_control_data[,1] 
})

#表达量变化超过2倍
fold_diff_express<-apply(diff_expr_data,1,function(x){
    sum(abs(x)>=log2(2))
})
##################################
#差异基因坐标
##################################
diff_express_data<-fold_diff_express[fold_diff_express>(dim(case_data)[2]/2)]
diff_expre_gene<-names(diff_express_data)
diff_gene_location<-c()
for (i in 1:length(diff_expre_gene))
{
    diff_gene_location<-c(diff_gene_location,which(row.names(case_data)==diff_expre_gene[i]) ) 
}

final_diff_data<-case_data[diff_gene_location,]
not_diff_data<-case_data[-diff_gene_location,]

control_final_diff_data<-control_data[diff_gene_location,]
control_not_diff_data<-control_data[-diff_gene_location,]
##################################
#求离散系数
##################################
Coefficient_of_Variation<-apply(not_diff_data,1,function(x){
    sd(x)/mean(x)
})

Coefficient_of_Variation2<-apply(control_not_diff_data,1,function(x){
    sd(x)/mean(x)
})
Coefficient_of_Variation3<-apply(control_final_diff_data,1,function(x){
    sd(x)/mean(x)
})
Coefficient_of_Variation4<-apply(final_diff_data,1,function(x){
    sd(x)/mean(x)
})

not_diff<-as.matrix(Coefficient_of_Variation,col=1)
diff<-as.matrix(Coefficient_of_Variation4,col=1)
cv_matrix<-rbind(not_diff,diff)
colnames(cv_matrix)<-'CV_values'
notdiff_tag<-as.matrix(rep(0,dim(not_diff)[1]),col=1)
diff_tag<-as.matrix(rep(1,dim(diff)[1]),col=1)
cv_tag<-rbind(notdiff_tag,diff_tag)
colnames(cv_tag)<-'CV_tag'
cv_draw<-cbind(cv_matrix,cv_tag)
cv_draw<-as.data.frame(cv_draw,stringsAsFactors=F)
cv_draw$CV_tag<-factor(cv_draw$CV_tag)
library(plyr)
cv_draw$CV_tag<-revalue(cv_draw$CV_tag,c('0'='Not diff expression','1'='Diff expression'))
#density=每个组里观测值的密度（占整体的百分数/组宽）
##################################
#确定MF,CC,BP对应的CV，DAVID_MF_divid_GO_term.txt这些数据通过EcoIN往DAVID进行映射获得
#DAVIDI注释分类.py 进行BP,CC,MF提取
##################################
MF_CV_TERM_NUM<-c()
BP_CV_TERM_NUM<-c()
CC_CV_TERM_NUM<-c()
MF_CV_GENE_NUM<-c()
BP_CV_GENE_NUM<-c()
CC_CV_GENE_NUM<-c()
total_CV_TERM_NUM<-c()
total_CV_GENE_NUM<-c()
setwd("C:\\workplace\\文章\\iRPBGO\\新流程\\EcoMAC")
MF_data<-read.table('DAVID_MF_divid_GO_term.txt',sep = '\t',stringsAsFactors = F,quote = '')
MF_go_gene_location<-c()
for (i in seq(0.01,0.20,0.005)){
    print(i)
    cut_CV<-Coefficient_of_Variation[Coefficient_of_Variation<=i]
    print (length(names(cut_CV)))
    for (j in 1:length(names(cut_CV)))
    {
        MF_go_gene_location<-c(MF_go_gene_location,which(MF_data[,1]==names(cut_CV)[j]) ) 
        
    }
    MF_go_gene_location<-unique(MF_go_gene_location)
    new_MF_data<-MF_data[MF_go_gene_location,]
    new_MF_table_term<-table(new_MF_data[,2])
    MF_CV_TERM_NUM<-c(MF_CV_TERM_NUM,length(new_MF_table_term))
    new_MF_table_gene<-table(new_MF_data[,1])
    MF_CV_GENE_NUM<-c(MF_CV_GENE_NUM,length(new_MF_table_gene))

}

MF_table_gene<-table(MF_data[,1])
MF_table_term<-table(MF_data[,2])
MF_table_termname<-table(MF_data[,3])

print(length(MF_table_term))

##################################
#做曲线图,选取合适的cv点，图中可以选择0.075
##################################
library(ggplot2)
library(gridExtra)
MF_dataframe<-data.frame(cv=seq(0.01,0.20,0.005),category=rep('MF',length(MF_CV_TERM_NUM)),termnum=MF_CV_TERM_NUM,genenum=MF_CV_GENE_NUM,gene2term_ratio=MF_CV_GENE_NUM/MF_CV_TERM_NUM,total_gene_ratio=MF_CV_GENE_NUM/length(MF_table_gene),total_term_ratio=MF_CV_TERM_NUM/length(MF_table_term))
p3<-ggplot(MF_dataframe,aes(cv,gene2term_ratio))+geom_path(aes(colour=factor(category)))+facet_grid(category~.,scales = 'free_y')#+geom_vline(xintercept = 0.1)
p4<-ggplot(MF_dataframe,aes(cv,total_gene_ratio))+geom_path(aes(colour=factor(category)))+facet_grid(category~.,scales = 'free_y')
p5<-ggplot(MF_dataframe,aes(cv,total_term_ratio))+geom_path(aes(colour=factor(category)))+facet_grid(category~.,scales = 'free_y')

grid.arrange(p4,p5,nrow=2)
