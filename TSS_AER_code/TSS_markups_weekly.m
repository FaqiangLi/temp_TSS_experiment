function TSS_markups_weekly(theta,inp,steplength,t)

% TSS_markups_weekly calculates elasticities, markups etc., price
% counterfactuals for 2000 households in week t.

% __________ Input arguments ______________________________________________
% theta:            estimated parameters (Dx1)
% inp:              input structure
% S:                9 or 16. The number of firms, depending on level of
%                   aggregation
% steplength:       step length for finite difference derivatives of demand
%                   with respect to price
% t:                week (1 to 156) for which markups etc. are computed
disp('Generating weekly data.')
inp1 = TSS_weekly_input_struct(inp,t);

R = 1;
disp('Starting calculation of markups and elasticities.')
[mcost,pr_w] = TSS_print_elasticities_markups(theta,inp1,R,steplength);
save mcost mcost pr_w
% disp('Warning: loading marginal cost')
% load mcost
disp('Starting calculation of counterfactuals [by exogenous consumer type].')
inp1.print = true;
TSS_counterfactual_exog_dist(theta,inp1,mcost,pr_w,steplength);


