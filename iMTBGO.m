function [fluxes, exp_ratio,opt_grow,model2]= iMTBGO(model2, gene_names, gene_exp,pheno)
% INPUTS
%       model - cobra model
%       gene_names - gene ids
%       gene_exp - gene expression
%
% OUTPUTS
%       fluxes - flux distribution
% Author: Mr.Mao, 2017 
    eps=1e-9;
    exp_ratio = gene_to_reaction_levels(model2, gene_names, gene_exp);
    exp_ratio(isnan(exp_ratio)) = 1;

    blocked_lb = model2.lb >= 0;
    blocked_ub = model2.ub <= 0;
    levels=exp_ratio;
    model2.lb = -levels;
    model2.lb(blocked_lb) = 0;
    model2.ub = levels;
    model2.ub(blocked_ub) = 0;
    
    model2=model_change(model2,pheno);
    sol = optimizeCbModel(model2);  
    scale_idx = strcmp('EX_glc(e)', model2.rxns);
    if isempty(sol.x)
        fluxes = [];
    else
        fluxes = sol.x* abs(pheno(1) / sol.x(scale_idx));
        fluxes(abs(fluxes)<eps)=0;
    end
    if isempty(sol.f)
        opt_grow = [];
    else
        opt_grow = sol.f* abs(pheno(1) / sol.x(scale_idx));
    end   
end