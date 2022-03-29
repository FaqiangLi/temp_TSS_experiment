%clear % clear workspace


clear

seed=12345
rng(seed)

%________________Settings__________________________________________________
N = 2000; % total number of households
outofsample = false; % = true for out-of-sample predictions
%   Moment conditions:       (281) , 
%                            of which
%       Continuous quantity: (146)      of which
%                            household (hh) size (8)
%                            time dummies (5*8 = 40)
%                            price IV (8)
%                            cross cat price IV (10)
%                            log storesize (8)
%                            chaindummies (8*8 = 64)
%                            constant (8)
%                            one-store store pair (8)
%       Discrete cat purchase: (106)    of which
%                            household (hh) size (8)
%                            price IV (8)
%                            cross cat price IV (10)
%                            log storesize (8)
%                            chaindummies (8*8 = 64)
%                            constant (8)
%                            one-store store pair (8)
%       Store pair choice:   (6)        of which
%                            distance (1)
%                            single-store store pair dummy (1)
%                            distance squared (1)
%                            distance*(single-store store pair dummy) (1)
%                            mean price IV across 2*8 categories (1)
%                            (mean price IV across 2*8 categories)/(per capita income) (1)
%      Cross-period moments: (5)
%    Cross-category moments: (2)
cat_comb = cell(4,1);
cat_comb{1} = [3 4];
cat_comb{2} = [8 2];
cat_comb{3} = [1 5 7];
cat_comb{4} = 6;

% cross-category instruments. 
IV_comb = cell(8,1);
IV_comb{1} = [5 7];     %     'Bakery' : F&V, Meat
IV_comb{2} = 8;         %     'Dairy': Milk
IV_comb{3} = 4;         %     'Drink': Dry
IV_comb{4} = 3;         %     'Dry grocery': Drink
IV_comb{5} = [1 7];     %     'Fruit & vegetables': Bakery, Meat
IV_comb{6} = [];        %     'Household goods': none
IV_comb{7} = [1 5];     %     'Meat': Bakery, F&V
IV_comb{8} = 2;         %     'Milk': Dairy

%_______________Get data ready_____________________________________________
% Create structures containing all the data etc. needed for estimation and
% post-estimation calculations
% inp: input arguments used in estimation
% inp2: input arguments used only after estimation, for calculating markups etc.
[inp,inp2] = TSS_input(N,cat_comb,IV_comb,outofsample);
disp('Data are ready.')
% number of parameters to be estimated
N = inp.N; J = inp.J; on = inp.on; S = inp.S; K = inp.K;
stchlabels = inp2.stchcl; L1 = size(stchlabels,1);
hhchlabels = inp2.hhchcl; L2 = size(hhchlabels,1);
rnlabels = inp2.rncl; Lr = size(rnlabels,1);
splabels = inp2.spcl; Lsp = size(splabels,1);
Li = inp.Li; 
% number of estimated parameters in the order in which they enter
D = K ...       % cat. scaling params. and constant
    + L1 ...    % store chars.
    + L2 ...    % household chars. (incl. time dummies)
    + Lr ...    % random terms in continuous utility
    + K ...     % second-order 'diagonal' terms
    + Li ...    % second-order interaction terms
    + 2 ...     % price terms
    + Lsp ...   % store pair terms
    + K*(S-1);  % firm-category effects (for all cats., all firms but ASDA)
disp(['There are ' num2str(inp.LL) ' moment conditions to estimate '...
    num2str(D) ' parameters.'])
inp.D = D;

%_______________Estimation with preliminary weighting matrix_______________

% specify new moments to be off
inp.new_moments = false % do not include cross-period moment conditions

% Simulated method of moment fit metrics
objective = @(theta)TSS_gmm(theta,inp);

% evaluate objective function at the true paramater
load('theta_health_predictable.mat');
theta_DGP = theta_health_predictable;
check_obj_val = objective(theta_DGP);

disp(['At true DGP with W=inv(ZZ), the objective function is ' num2str(check_obj_val) ] )


% Load preliminary estimates to generate initial bounds for optimization
% algorithm. Any Dx1 vector can be used. Alternatively bounds can be specified
% directly. See TSS_optimization.
% initial bound close to true DGP (a cheat)
theta=theta_DGP; % check TSS_optimization_1 for further details

% display initial theta
disp('Initial bound is (true theta +- 0.5 |true theta|) : ')
theta

% preliminary estimation
TSS_optimization_1 % this is a separate script file instead of a function
% file in order to permit flexibility in running the estimation, since
% computation time is very long.

disp('Prelim stage estimates is : ')
theta_prelim_stg = theta;

% show estimates in diary
[theta_header theta_DGP theta_prelim_stg]


%_______________Estimation with 1st-stage weighting matrix_________________

% Add in new moments and evaluate new weighting matrix (note that W2
% incorporates inp.new_moments = true
inp.new_moments = true
W2 = TSS_W_second_stage_modified(theta,inp,inp2);

inp.W = inv(W2);
objective = @(theta)TSS_gmm(theta,inp);

% Re-evaluate objective function at the updated weighting matrix
load('theta_health_predictable.mat');
theta_DGP = theta_health_predictable;
check_obj_val = objective(theta_DGP);

disp(['At true DGP with W=prelim-optimal, the objective function is ' num2str(check_obj_val) ] )

% initial bound for first stage
disp('Initial bound for the first stage is (prelim stage estimates):')
theta


% For first and second stage estimates, I use both TSS_optimization_2
TSS_optimization_2
disp('First stage estimates is ')
theta_first_stg = theta;
[theta_header theta_DGP theta_first_stg]


% Calculate asymtotic var-cov matrix
% caveat to prof. Tybout: in the March 17 email report, this following part
% is skipped. 
acov_1stg = TSS_acov(theta,inp,inp2);
se_1stg = sqrt(diag(2*acov_1stg)); %  2* to correct for simulation noise
theta = reshape(theta,D,1);
se_1stg = reshape(se_1stg,D,1);

[theta_header theta_DGP theta_first_stg se_1stg]

% print parameter estimates and standard errors to excel file named
% 'ParameterEstimates_ddmm_yyyy_HHMM.xlsx'
% faqiang note: this seems to use outdated Matlab command.
% TSS_print_estimates(theta,se_1stg,inp,inp2)




%_______________Estimation with 2nd-stage weighting matrix_________________
inp.new_moments = true
% load first_stage_estimates
W2 = TSS_W_second_stage(theta,inp,inp2);
inp.W = inv(W2);
objective = @(theta)TSS_gmm(theta,inp);

load('theta_health_predictable.mat');
theta_DGP = theta_health_predictable;
check_obj_val = objective(theta_DGP);

disp(['At true DGP with W=firststage-optimal, the objective function is ' num2str(check_obj_val) ] )
disp('Initial bound for the second stage is (first stage estimates):')

TSS_optimization_2
% save second_stage_estimates

%________________Standard errors___________________________________________
acov = TSS_acov(theta,inp,inp2);
se = sqrt(diag(2*acov)); %  2* to correct for simulation noise
theta = reshape(theta,D,1);
se = reshape(se,D,1);
% print parameter estimates and standard errors to excel file named
% 'ParameterEstimates_ddmm_yyyy_HHMM.xlsx'
TSS_print_estimates(theta,se,inp,inp2)

%_________________Predictions______________________________________________
% print predictions to excel file named 'Predictions_ddmm_yyyy_HHMM.xlsx'
% sheet 1: firm/category predictions
% sheet 2: firm-level avg. hh.chars. & market shares in #vis. & rev., by
%           1SS and 2SS
% sheet 3: firm combinations: visitors (obs. & pred., separate tables)
% sheet 4: firm combinations: revenue (obs. & pred., separate tables)
% sheet 5: firm/category predictions - validation sample
TSS_print_predictions(theta,inp,inp2)

%___________week by week calculations (markups etc., counterfactuals)______ 
inp2.A = [];
steplength = 5e-5;
inp.frmcl = inp2.frmcl; 
inp.catcl = inp2.catcl;
inp.chain2 = inp2.chain2;
inp.chain = inp2.chain;

t = 78; % week (out of 156) from sample
disp(['Week number ' num2str(t)])
TSS_markups_weekly(theta,inp,steplength,t)

%_____________ Generate distribution of median markups ____________________
S = 16;     % number of firms
NS = 2000;  % number of draws for confidence intervals
[markups,markups_deleg,markups_etc_full] = TSS_markups_counterf_distribution(theta,2*acov,inp,S,NS,steplength);