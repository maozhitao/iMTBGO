setwd("..\\iMTBGO\\")

data_ori<-read.table('EcoMAC_matrixg.txt',header = T,sep = '\t')
#ID redundant
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

#control
control_name<-read.table('control_data_M9_e.txt',sep = '\t')
control_location<-c()
for (i in 1:length(control_name[,1]))
{
    control_location<-c(control_location,which(colnames(data)==control_name[i,1]) ) 
}
case_data<-data[,-c(1,control_location)]
control_data<-data[,control_location]
row.names(control_data)<-data[,1]
row.names(case_data)<-data[,1]
use_control_data<-matrix(apply(control_data,1,mean),nrow = dim(data)[1],dimnames = list(data[,1],'WT_M9'))

write.table(case_data,"case_matrix.txt",quote = F,sep='\t',col.names = T,row.names = T)
write.table(control_data,"control_matrix.txt",quote = F,sep='\t',col.names = T,row.names = T)

##################################
#diff expresiion
##################################
diff_expr_data<-apply(case_data,2,function(x){
    x-use_control_data[,1] 
})

#2 fold-change
fold_diff_express<-apply(diff_expr_data,1,function(x){
    sum(abs(x)>=log2(2))
})

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
#CV
##################################
Coefficient_of_Variation<-apply(not_diff_data,1,function(x){
    sd(x)/mean(x)
})
Coefficient_sd<-apply(not_diff_data,1,function(x){
    sd(x)
})
Coefficient_mean<-apply(not_diff_data,1,function(x){
    mean(x)
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

Coefficient_sd4<-apply(final_diff_data,1,function(x){
    sd(x)
})
Coefficient_mean4<-apply(final_diff_data,1,function(x){
    mean(x)
})


not_diff<-as.matrix(Coefficient_sd,col=1)
diff<-as.matrix(Coefficient_sd4,col=1)
sd_matrix<-rbind(not_diff,diff)
colnames(sd_matrix)<-'sd_values'

not_diff<-as.matrix(Coefficient_mean,col=1)
diff<-as.matrix(Coefficient_mean4,col=1)
mean_matrix<-rbind(not_diff,diff)
colnames(mean_matrix)<-'mean_values'

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
cv_draw$CV_tag<-revalue(cv_draw$CV_tag,c('0'='Stable expression','1'='Not stable expression'))
library(ggplot2)
ggplot(cv_draw,aes(x=CV_values))+geom_freqpoly(aes(y=..density..,colour=CV_tag),binwidth=0.01)
ggplot(cv_draw,aes(x=CV_values))+geom_freqpoly(aes(y=..count..,colour=CV_tag),binwidth=0.01)
ggplot(cv_draw,aes(x=CV_values))+geom_histogram(aes(fill=CV_tag),binwidth=0.01,position = 'fill')
ggplot(cv_draw,aes(x=CV_values))+geom_histogram(aes(y=..density..),binwidth=0.01)+facet_grid(CV_tag~.)
ggplot(cv_draw,aes(x=CV_values))+geom_histogram(aes(y=..count..),binwidth=0.01)+facet_grid(CV_tag~.)
ggplot(cv_draw,aes(x=CV_values),fill=CV_tag)+geom_density(alpha=0.3)+facet_grid(CV_tag~.)


GO_type='MF'#'MF''total'
MF_CV_TERM_NUM<-c()
MF_CV_GENE_NUM<-c()
filename<-paste('DAVID',GO_type,'divid_GO_term.txt',sep="_");
MF_data<-read.table(filename,sep = '\t',stringsAsFactors = F,quote = '')
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
    print(MF_CV_TERM_NUM)
    new_MF_table_gene<-table(new_MF_data[,1])
    MF_CV_GENE_NUM<-c(MF_CV_GENE_NUM,length(new_MF_table_gene))
}

MF_table_gene<-table(MF_data[,1])
MF_table_term<-table(MF_data[,2])
MF_table_termname<-table(MF_data[,3])
print(length(MF_table_term))


library(ggplot2)
library(gridExtra)
MF_dataframe<-data.frame(cv=seq(0.01,0.20,0.005),term=rep('Go_term',length(MF_CV_TERM_NUM)),termnum=MF_CV_TERM_NUM,genenum=MF_CV_GENE_NUM,gene2term_ratio=MF_CV_GENE_NUM/MF_CV_TERM_NUM,BP_gene_ratio=MF_CV_GENE_NUM/length(MF_table_gene),BP_term_ratio=MF_CV_TERM_NUM/length(MF_table_term))
p3<-ggplot(MF_dataframe,aes(cv,gene2term_ratio))+geom_path()+facet_grid(term~.,scales = 'free_y')
p4<-ggplot(MF_dataframe,aes(cv,BP_gene_ratio))+geom_path(aes(colour=factor(term)))+facet_grid(term~.,scales = 'free_y')
p5<-ggplot(MF_dataframe,aes(cv,BP_term_ratio))+geom_path(aes(colour=factor(term)))+facet_grid(term~.,scales = 'free_y')
p6<-ggplot(cv_draw,aes(x=CV_values))+geom_freqpoly(aes(y=..density..,colour=CV_tag),binwidth=0.01)
grid.arrange(p4,p5,p3,nrow=3)
grid.arrange(p4,p5,p3,p6,nrow=2)
grid.arrange(p4,p5,p3,p6,nrow=2)

#cv_thrd=0.1
cv=seq(0.04,0.20,0.005)
for (cv_thrd in cv){
MF_term_location<-c()
miss_MF_table_term<-c()
MF_go_gene_location<-c()
cut_CV<-Coefficient_of_Variation[Coefficient_of_Variation<=cv_thrd]
cut_CV_no_use<-Coefficient_of_Variation[Coefficient_of_Variation>cv_thrd]
for (i in 1:length(names(cut_CV)))
{
    MF_go_gene_location<-c(MF_go_gene_location,which(MF_data[,1]==names(cut_CV)[i]) ) 
}
#MF_go_gene_location<-unique(MF_go_gene_location)
new_MF_data<-MF_data[MF_go_gene_location,]
MF_table_term<-table(MF_data[,2])
new_MF_table_term<-table(new_MF_data[,2])
not_MF_data<-MF_data[-MF_go_gene_location,]
not_MF_table_term<-table(not_MF_data[,2])

setwd("..\\iMTBGO\\sensitivity_analysis\\")
MF_data_out<-data.frame()
MF_marker_gene_data_out<-data.frame()
#
for (i in 1:length(names(new_MF_table_term)))
{
    
    MF_location<-c()
    MF_location<-c(MF_location,which(new_MF_data[,2]==names(new_MF_table_term)[i])) 
    GO_location<-c(MF_location,which(MF_data[,2]==names(new_MF_table_term)[i])) 
    
    cut_cv_loc<-c()
    for (mf_gene in new_MF_data[MF_location,1]){
        #print(mf_gene)
        cut_cv_loc<-c(cut_cv_loc,which(names(cut_CV)==mf_gene))         
    }

    MF_data_out[i,1]<-names(new_MF_table_term)[i]
    MF_data_out[i,2]<-names(which.min(cut_CV[cut_cv_loc]))
    MF_data_out[i,3]<-min(cut_CV[cut_cv_loc])   
    for (eachgo in GO_location){
        MF_marker_gene_data_out[i,1]<-names(new_MF_table_term)[i]
        MF_marker_gene_data_out[i,2]<-MF_data[eachgo,1]
        MF_marker_gene_data_out[i,3]<-names(which.min(cut_CV[cut_cv_loc]))
        MF_marker_gene_data_out[i,4]<-min(cut_CV[cut_cv_loc])           
    }

}

#term no marker 
MF_data_out[i,1]<-'other'
MF_data_out[i,2]<-names(which.min(Coefficient_of_Variation))
MF_data_out[i,3]<-min(Coefficient_of_Variation)

colnames(MF_data_out)<-c("term_id","marker_gene",'cv_values')
filename<-paste(cv_thrd,GO_type,'marked_gene2go.txt',sep="_");
filename2<-paste(cv_thrd,GO_type,'gene_marker.csv',sep="_");
colnames(MF_marker_gene_data_out)<-c("term_id","gene","marker_gene",'cv_values')
write.table(MF_data_out,filename,quote = F,sep='\t',col.names = T,row.names = F)
write.table(MF_marker_gene_data_out,filename2,quote = F,sep=',',col.names = T,row.names = F)


gene2marker_uni=data.frame()

for (eachgene in unique(MF_marker_gene_data_out[,'gene'])){
    gene2marker_loc<-c()
    gene_loc=which(MF_marker_gene_data_out[,'gene']==eachgene)
    gene2marker_uni[eachgene,1]=eachgene
    gene2marker_uni[eachgene,2]<-MF_marker_gene_data_out[gene_loc,'marker_gene'][which.min(MF_marker_gene_data_out[gene_loc,'cv_values'])]
    gene2marker_uni[eachgene,3]<-min(MF_marker_gene_data_out[gene_loc,'cv_values'])  
}
filename3<-paste(cv_thrd,GO_type,'gene_marker_uni.csv',sep="_");
colnames(gene2marker_uni)<-c("gene","uni_marker_gene",'cv_values')
write.table(gene2marker_uni,filename3,quote = F,sep=',',col.names = T,row.names = F)


}


##################################
#specific profils
##################################
gseid="Ishii2007"
#gseid="Holm2010"

oldwd<-"..\\iMTBGO\\expression_file\\Ishii2007"
#oldwd<-"..\\iMTBGO\\expression_file\\Holm2010"

setwd(oldwd)

filename<-paste(gseid,"expr_id_trans.txt",sep="_");
data_ori<-read.table(filename,header = T,sep = '\t')

table_data<-table(data_ori[,1])
redu_gene_list<-names(table_data[table_data>1])
if (length(redu_gene_list)>0){
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
    tmdata<-rbind(unidata,redu_gene_matrix)
    useeset<-tmdata[,-1]
    row.names(useeset)<-tmdata[,1]
    #useeset<-log2(useeset);
    for (usei in 1:dim(useeset)[1])
    {
        useeset[usei,which(is.na(useeset[usei,]))]=apply(useeset[usei,],1,function(x){
            mean(x[which(!is.na(x))])
        }
        )  
    }  
}else{
    useeset=data_ori[,2:dim(data_ori)[2]]
    row.names(useeset)<-data_ori[,1]
}

##################################
#marker gene
##################################
cv=seq(0.04,0.20,0.005)
for (cv_thrd in cv){
setwd("..\\iMTBGO\\sensitivity_analysis\\")
infilename<-paste(cv_thrd,GO_type,'gene_marker_uni.csv',sep="_");
gene2marker=read.csv(infilename)
rownames(gene2marker)=gene2marker[,'gene']

outfile_list<-c()
outnamefile_list<-c()
GSE_detail<-colnames(useeset)
for (GSEj in 1:length(GSE_detail)){
    avg_use_expr<-as.matrix(useeset[,GSEj],ncol=1)
    row.names(avg_use_expr)<-row.names(useeset)
    colnames(avg_use_expr)<-'avg_expr'
    
    avg_use_expr_plus<-as.data.frame(avg_use_expr)

    norm_value<-c()
    for (geneid in rownames(avg_use_expr)){
        if (is.na(gene2marker[geneid,'uni_marker_gene'])){
            norm_value<-c(norm_value,avg_use_expr['b3033',])
        }
        else{
            norm_value<-c(norm_value,avg_use_expr[gene2marker[geneid,'uni_marker_gene'],])
        }
    }
    avg_use_expr_plus$norm_value<-as.numeric(norm_value) 
    avg_use_expr_plus$expr_ratio<-avg_use_expr_plus$avg_expr/avg_use_expr_plus$norm_value
    max_avg_expr<-max(avg_use_expr_plus$avg_expr)
    avg_use_expr_plus$eflux_expr_ratio<-avg_use_expr_plus$avg_expr/max_avg_expr
    avg_use_expr_plus[is.na(avg_use_expr_plus)]<-0#
    gene_name_file<-as.matrix(row.names(avg_use_expr_plus),ncol=1)
    outgenename<-paste(cv_thrd,GO_type,GSE_detail[GSEj],'gene_name.txt',sep = '_')
    outfilename<-paste(cv_thrd,GO_type,GSE_detail[GSEj],'expr_data.txt',sep = '_')
    setwd("..\\iMTBGO\\sensitivity_analysis\\matlab_infile\\")
    write.table(avg_use_expr_plus,outfilename,quote = FALSE,sep='\t',row.names = F,col.names = F)
    write.table(gene_name_file,outgenename,quote = FALSE,sep='\t',row.names = F,col.names = F)  
    outfile_list<-c(outfile_list,outfilename)
    outnamefile_list<-c(outnamefile_list,outgenename)
}

setwd("..\\iMTBGO\\sensitivity_analysis\\matlab_infile\\")
filename1<-paste(gseid,cv_thrd,GO_type,'outgenenamefile_list.txt',sep="_");
filename2<-paste(gseid,cv_thrd,GO_type,'outfile_list.txt',sep="_");
outnamefile_list<-as.data.frame(outnamefile_list)
write.table(outnamefile_list,filename1,quote = FALSE,sep='\t',row.names = F,col.names = F)
outfile_list<-as.data.frame( outfile_list)
write.table(outfile_list,filename2,quote = FALSE,sep='\t',row.names = F,col.names = F)
}