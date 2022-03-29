%% Simulate the endogenous variables from the exogenous variables

% One needs to custom one things before running: N_select: This is the size
% for the sample containing all the candidate households for structural
% estimation.

% To see the parameter value, see parparam_predictable_nosample_selection.m

clear



% seed is 123457 throughout in this script
seed = 123457 % this is the seed for generating data
rng(seed);


% The following process are stochastic in generating the data. Therefore
% they are all affected by the seed.
% 1.Generate random household heterogeneity (\nu) in taste used to simulate
% endogenous variables
% 2.Randomly select N_select out of 7500 households to form the estimation
% candidate sample
% 3.Draw choices for all the 7500 households
%


% Diary
DiaryName = ['exprmnt','_',datestr(now,'yyyy_mmdd_HHMM'),'_','content'];
DiaryName = [pwd ,'/diary/', DiaryName]
diary(DiaryName)
generate_data_again=1

% Workspace that stores all the output of the simulation, namely the
% simulated (endogenous) varialbe and the exogenous data for structural
% estimation. The output will be formated as in TSS's paper.
output_dir = [pwd,'/neoexperiment_6000_real_dgpfixedfromraw']

if generate_data_again==1

    % Directory storing input, i.e., raw exogeneous variable used to generate simulated data.
    input_dir = [pwd,'/fq_matrix'];

    % Load the parameters used to simulate data (stored in the main
    % workspace ) 
    cd(pwd);
    load('theta_health_predictable.mat')
    params = theta_health_predictable;   

    % Select the number of the candidate households we expect to feed into
    % the TSS estimation routine.
    N_select = 6000;  % N_select needs to be <7500, the largest number of candidate households for structural estimation

    % Generate data. The function has no output. The formated data will be
    % saved in output_dir.  
    simulate_data_mainscript(params,N_select,input_dir,output_dir,seed)
    
end
