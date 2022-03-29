%%% Generate true DGP parameters, from TSS 

% Load TSS's estiamtes. 
% The reason is that there are some parameter values not provided in the
% paper so we just set those to be specified by their first stage
% estimates.
clear


load('first_stage_estimates_TSS.mat')

% cat taste parameters 
beta0=[2.085 1.418 1.096 1.901 2.665 1.115 2.309];

% baseline category (milk)
% Sept 8: change this to one for this is the normalized category.
contant0=1;

% store parameters (betax)
beta1=0.456; % log floor space

% household parameters (betaz)
%  beta3 and quarter fixed effect are included here,from 'first_stage_estimates_forGenData.mat'
beta2=0.477; % household size
beta3=theta(10); %  I do not know 
beta_quarter=theta(11:14)'; % quarter fixed effect, just use the initial estimate

% unobserved individual heterogeneity
sigma=[0.208 0.920 0.677 1.129]; % cat/store, time varying, cat specific, store&cat specific

% qudratic parameters (cross cat complementarity)
Lambda_diag=[19.852 11.239 3.802 10.536 15.952 4.360 8.901 14.062]*0.7;
Lambda_interact=[1.742 1.368 0.269 0.572 0.076]*0.7; % model 2 complementarity between categoreis 34,28,15,17,57  

% price parameter
alpha1=1.839; % constant
alpha1=alpha1*0.7;
alpha2=32.881; % 1/[weekl income per head]
alpha2=alpha2*0.7;

% shopping cost
gamma1=[7.528 0.440]; % two-store dummy, sd
gamma2=[10.269 0.394]; % distance, sd

% category-firm fixed effect
firmcat = theta(39:end)';
temp_id = find(firmcat<0);
firmcat(temp_id)= 0 - firmcat(temp_id);
% increase Soriana and Aurerea
mulitplier = 7;
firmcat(1:8) = mulitplier*firmcat(1:8) ;
firmcat(9:16) = 1*firmcat(9:16) ;


% delete the confounding 
clear theta fval_ga_mat

% put all parameter into a D*1 vector, D=length(theta)
theta=[beta0 contant0 beta1 beta2 beta3 beta_quarter sigma Lambda_diag Lambda_interact alpha1 alpha2 gamma1 gamma2 firmcat]';
theta_header = ["Bakery(all are TSS's cat @@)";"Diary";"Drink";"Dry Groc";"Fruit&Veg";"HH goods";"Meat"; ...
   "Milk(baseline cat)" ; ...
   "ln(floor size) @@ ";  ...
   "HH size"; ...
   "beta3";  ...
   "Time FE1";  "Time FE2"; "Time FE3"; "Time FE4"; ...
   "sigma_i";   "sigma_it"; "sigma_ik"; "sigma_ijk"; ...
   "Bakery-Bakery @@"; "Dairy-Dairy"; "Drink-Drink"; "Dry Groc-Dry Groc"; "Fruit&Veg-Fruit&Veg"; "HH good-HH good"; "Meat-Meat"; "Milk-Milk"; ...
   "Drink-Dry Groc";  "Diary-Milk";  "Bakery-Fruit&Veg";  "Bakery-Meat";  "Fruit&Veg-Meat"; ... 
   "Price-sensitvity, the constant part";  ... 
   "Price-sensitvity, the variable part";  ... 
   "Two store dummy"; "sd "; ... 
   "Distance" ;  "sd " ;... 
   "firmcat @@" ;   strings(length(firmcat)-1,1); ... 
   ];


% give the theta a different name
theta_health_morepredictable=theta;

save theta_health_morepredictable theta_health_morepredictable theta_header




