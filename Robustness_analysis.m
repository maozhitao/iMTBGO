clear
cd '..\iMTBGO\'
%exp_name='Holm2010'
%biomass_num=7;
GO_type='MF';
%GO_type='total';
exp_name='Ishii2007'
biomass_num=9;
num_seq=0:0.25:1;
savei=1;
fluxes_error_matrix=[];
biomass_error_matrix=[];
changeCobraSolver('glpk');
numi=0.175
for coef= num_seq
    filename=strcat('..\iMTBGO\Robustness_analysis\matlab_infile\',exp_name,'_',num2str(numi),'_',num2str(coef),'_',GO_type,'_outgenenamefile_list.txt');
    filename2=strcat('..\iMTBGO\Robustness_analysis\matlab_infile\',exp_name,'_',num2str(numi),'_',num2str(coef),'_',GO_type,'_outfile_list.txt');
    filename3=strcat('..\iMTBGO\expression_file\',exp_name,'\phenotype.txt');
    filename4=strcat('..\iMTBGO\',exp_name,'_centre_reactiong.mat');
    path='..\iMTBGO\Robustness_analysis\matlab_infile\';
    outgenenamefile_list=import_myfile(filename);
    outfile_list=import_myfile(filename2);
    phenotype = import_phenotype(filename3);
    load('ijo1366.mat')
    load(filename4);
    
    total_centre_flux={};
    totla_flux={};
    j=0;
    exp_ratio_matrix={};
    for i=1:length(outfile_list)
        expr_file=strcat(path,outfile_list{i});
        genename_file=strcat(path,outgenenamefile_list{i});
        my_gene_ratio=import_exp_data(expr_file);
        genename = import_gene_name(genename_file);
        model1=model;
        [fluxes, exp_ratio,opt_grow,newmodel]= iMTBGO(model1, genename,my_gene_ratio,phenotype(i,:));
        exp_ratio_matrix{i}=exp_ratio;
        for ii=1:length(centre_reaction_list)
            [x , y] = find(strcmp(model.rxns,centre_reaction_list{ii}));
            total_centre_flux{ii,j+1}=fluxes(x,y);
        end
        total_flux{1,i}=fluxes;
        j=j+1;
    end
    filename5=strcat('..\iMTBGO\',exp_name,'_exp_data.mat');
    load(filename5);
    for i=1:length(outfile_list)
        %fluxes error
        cal_error=norm(cell2mat(total_centre_flux(:,i))-cell2mat(exp_data(:,i)))./norm(cell2mat(exp_data(:,i)));
        %biomass error
        cal_error2=abs(cell2mat(total_centre_flux(biomass_num,i))-cell2mat(exp_data(biomass_num,i)))/abs(cell2mat(exp_data(biomass_num,i)));
        fluxes_error_matrix(i,savei)=cal_error;
        biomass_error_matrix(i,savei)=cal_error2;
    end
    savei=savei+1;
    save_file=strcat('..\iMTBGO\Robustness_analysis\iMBTGO\',exp_name,'_',GO_type,'_',num2str(numi),'_robust_results.mat');
    save(save_file)      
end
filename6=strcat('..\iMTBGO\Robustness_analysis\iMBTGO\',exp_name,'_iMBTGO_robust_results');

set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.55]);
Average=mean(fluxes_error_matrix);  
Variance=var(fluxes_error_matrix,1);   
errorbar(num_seq,Average,Variance)    
xlabel('Noise');ylabel('Error');
%saveas(gcf,filename6,'jpg')
save_file=strcat(filename6,'_',GO_type,'_fluxes.jpg');
print(gcf,'-djpeg',save_file) 

set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.55]);
Average=mean(biomass_error_matrix);  %
Variance=var(biomass_error_matrix,1);   %
errorbar(num_seq,Average,Variance)    %
xlabel('Noise');ylabel('Error');
%saveas(gcf,filename6,'jpg')
save_file=strcat(filename6,'_',GO_type,'_biomass.jpg');
print(gcf,'-djpeg',save_file) %

save_file=strcat(filename6,'_',GO_type,'.mat');
save(save_file)  