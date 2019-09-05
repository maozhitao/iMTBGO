function [newmodel]= model_change(model,pheno)

    if ~isnan(pheno(1))
        model = changeRxnBounds(model,'EX_glc(e)',-pheno(1),'l');
    end
    if ~isnan(pheno(2))  
        if pheno(2)==0
            model = changeRxnBounds(model,'EX_o2(e)',0,'l');%
        else
            model = changeRxnBounds(model,'EX_o2(e)',-pheno(2),'l');
        end
    end
    newmodel=model;