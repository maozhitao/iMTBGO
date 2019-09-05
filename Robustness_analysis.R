
gseid="Ishii2007"
#gseid="Holm2010"
GO_type='MF'#'total'#
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
#cv_thrd=0.13
cv_thrd=0.175
coef_seq=seq(0,1,0.25)
for (coef in coef_seq){
    setwd("..\\iMTBGO\\sensitivity_analysis\\")
    infilename<-paste(cv_thrd,GO_type,'gene_marker_uni.csv',sep="_");
    gene2marker=read.csv(infilename)
    rownames(gene2marker)=gene2marker[,'gene']
    outfile_list<-c()
    outnamefile_list<-c()
    GSE_detail<-colnames(useeset)
    #random sample
    nrow=dim(useeset)[1]
    ncol=dim(useeset)[2]
    eset_sample=c()
    cc=c()
    for (ncoli in 1:ncol) {
        eset_sample=c(eset_sample,sample(useeset[,ncoli],size=nrow))
        cc=c(cc,ncoli)
    }
    newuseeset=matrix(eset_sample,ncol=ncol,nrow=nrow)
    row.names(newuseeset)<-row.names(useeset)
    colnames(newuseeset)<-colnames(useeset)

    useuseeset=useeset+coef*(newuseeset-useeset)
    for (GSEj in 1:length(GSE_detail)){
        avg_use_expr<-as.matrix(useuseeset[,GSEj],ncol=1)
        row.names(avg_use_expr)<-row.names(useuseeset)
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
        avg_use_expr_plus[is.na(avg_use_expr_plus)]<-0
        gene_name_file<-as.matrix(row.names(avg_use_expr_plus),ncol=1)
        outgenename<-paste(cv_thrd,coef,GO_type,GSE_detail[GSEj],'gene_name.txt',sep = '_')
        outfilename<-paste(cv_thrd,coef,GO_type,GSE_detail[GSEj],'expr_data.txt',sep = '_')
        setwd("..\\iMTBGO\\Robustness_analysis\\matlab_infile\\")
        write.table(avg_use_expr_plus,outfilename,quote = FALSE,sep='\t',row.names = F,col.names = F)
        write.table(gene_name_file,outgenename,quote = FALSE,sep='\t',row.names = F,col.names = F)  
        outfile_list<-c(outfile_list,outfilename)
        outnamefile_list<-c(outnamefile_list,outgenename)
    }
    
    setwd("..\\iMTBGO\\Robustness_analysis\\matlab_infile\\")
    filename1<-paste(gseid,cv_thrd,coef,GO_type,'outgenenamefile_list.txt',sep="_");
    filename2<-paste(gseid,cv_thrd,coef,GO_type,'outfile_list.txt',sep="_");
    outnamefile_list<-as.data.frame(outnamefile_list)
    write.table(outnamefile_list,filename1,quote = FALSE,sep='\t',row.names = F,col.names = F)
    outfile_list<-as.data.frame( outfile_list)
    write.table(outfile_list,filename2,quote = FALSE,sep='\t',row.names = F,col.names = F)
}