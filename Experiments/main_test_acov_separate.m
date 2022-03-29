%% Extra script to get variance covariance matrix separately from the main script


clear

seed=12345
rng(seed)

DiaryName = ['exprmnt_acov','_',datestr(now,'yyyy_mmdd_HHMM'),'_','first_stage_acov'];
DiaryName = [pwd ,'/diary/', DiaryName]
diary(DiaryName)

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


% Get the correct weighting matrix


% Add in new moments and evaluate new weighting matrix (note that W2
% incorporates inp.new_moments = true
outputdir = "/storage/work/f/fxl146/WalmexNutrition/Experiments/Output"
filename = outputdir+'/1stg_2022_0314_1018_Run20of20.mat'
load(filename)
theta_for_weighting = x';
inp.new_moments = true
W2 = TSS_W_second_stage_modified(theta_for_weighting,inp,inp2);
inp.W = inv(W2);
objective = @(theta)TSS_gmm(theta,inp);

% Re-evaluate objective function and compare this to the log file to make
% sure we are loading the correct weighting matrix used during the
% estimation process we expect to evaluate.
load('theta_health_predictable.mat');
theta_DGP = theta_health_predictable;
check_obj_val = objective(theta_DGP);

disp(['At true DGP,the objective function is ' num2str(check_obj_val) ] )
outputdir = "/storage/work/f/fxl146/WalmexNutrition/Experiments/Output"
filename = outputdir+'/2stg_2022_0315_0703_Run20of20.mat'
load(filename)
theta_for_acov = x';
acov_temp = TSS_acov(theta_for_acov,inp,inp2);
se_temp = sqrt(diag(2*acov_temp)); %  2* to correct for simulation noise
theta_for_acov = reshape(theta_for_acov,D,1);
se_temp = reshape(se_temp,D,1);
% print parameter estimates and standard errors to excel file named
% 'ParameterEstimates_ddmm_yyyy_HHMM.xlsx'
[theta_header theta_DGP theta_for_acov se_temp]
 TSS_print_estimates(theta_for_acov,theta_for_acov,inp,inp2)

% 
%  102×4 string array
% 
%     "Bakery(all are TSS's cat @@)"       "2.085"          "1.8182"          "0.052158" 
%     "Diary"                              "1.418"          "1.8686"          "0.034739" 
%     "Drink"                              "1.096"          "1.22721"         "0.024898" 
%     "Dry Groc"                           "1.901"          "2.53504"         "0.055814" 
%     "Fruit&Veg"                          "2.665"          "3.0515"          "0.069554" 
%     "HH goods"                           "1.115"          "1.20648"         "0.028129" 
%     "Meat"                               "2.309"          "2.80069"         "0.053879" 
%     "Milk(baseline cat)"                 "1"              "1.53424"         "0.033641" 
%     "ln(floor size) @@ "                 "0.456"          "0.503321"        "0.0072771"
%     "HH size"                            "0.477"          "0.111132"        "0.0034365"
%     "beta3"                              "0.497866"       "0.381434"        "0.037347" 
%     "Time FE1"                           "0.169787"       "0.176429"        "0.028673" 
%     "Time FE2"                           "0.27256"        "0.183037"        "0.034118" 
%     "Time FE3"                           "0.281607"       "0.302216"        "0.013959" 
%     "Time FE4"                           "0.162685"       "0.152345"        "0.027696" 
%     "sigma_i"                            "0.208"          "0.0698468"       "0.017356" 
%     "sigma_it"                           "0.92"           "0.0177384"       "0.0063583"
%     "sigma_ik"                           "0.677"          "-0.0319006"      "0.023762" 
%     "sigma_ijk"                          "1.129"          "-0.000347353"    "0.21971"  
%     "Bakery-Bakery @@"                   "13.8964"        "66.072"          "3.21"     
%     "Dairy-Dairy"                        "7.8673"         "30.7819"         "1.3084"   
%     "Drink-Drink"                        "2.6614"         "14.7888"         "0.6747"   
%     "Dry Groc-Dry Groc"                  "7.3752"         "37.4628"         "1.6387"   
%     "Fruit&Veg-Fruit&Veg"                "11.1664"        "41.5436"         "2.2887"   
%     "HH good-HH good"                    "3.052"          "12.4992"         "0.59193"  
%     "Meat-Meat"                          "6.2307"         "28.753"          "1.2204"   
%     "Milk-Milk"                          "9.8434"         "16.7472"         "0.84693"  
%     "Drink-Dry Groc"                     "1.2194"         "0.585139"        "0.42236"  
%     "Diary-Milk"                         "0.9576"         "1.70068"         "0.50315"  
%     "Bakery-Fruit&Veg"                   "0.1883"         "0.0526135"       "0.33038"  
%     "Bakery-Meat"                        "0.4004"         "0.117721"        "0.25251"  
%     "Fruit&Veg-Meat"                     "0.0532"         "0.0287567"       "0.86894"  
%     "Price-sensitvity, the constan…"    "1.2873"         "1.90129"         "0.046594" 
%     "Price-sensitvity, the variabl…"    "23.0167"        "0.883177"        "0.32787"  
%     "Two store dummy"                    "7.528"          "4.50604"         "0.53826"  
%     "sd "                                "0.44"           "0.053377"        "0.0071475"
%     "Distance"                           "10.269"         "18.0652"         "5.2662"   
%     "sd "                                "0.394"          "0.0103759"       "0.031039" 
%     "firmcat @@"                         "7.08961"        "12.5168"         "0.37746"  
%     ""                                   "0.143079"       "0.115782"        "0.068215" 
%     ""                                   "4.06846"        "3.63326"         "0.24845"  
%     ""                                   "17.4459"        "26.5577"         "0.66772"  
%     ""                                   "0.0393509"      "0.171296"        "0.068127" 
%     ""                                   "0.0393701"      "0.00266454"      "0.18222"  
%     ""                                   "4.62551"        "6.79638"         "0.25674"  
%     ""                                   "1.10369"        "0.413871"        "0.17221"  
%     ""                                   "0.251555"       "0.420682"        "0.057518" 
%     ""                                   "0.18654"        "-0.0764702"      "0.13111"  
%     ""                                   "0.000553382"    "0.0011825"       "0.081779" 
%     ""                                   "0.0179749"      "0.0170135"       "0.17289"  
%     ""                                   "0.0185715"      "0.0121061"       "0.18272"  
%     ""                                   "0.0687266"      "0.0393204"       "0.064335" 
%     ""                                   "0.0875933"      "0.0726606"       "0.066789" 
%     ""                                   "0.0556822"      "0.137872"        "0.064004" 
%     ""                                   "0.00682609"     "0.00217701"      "0.061265" 
%     ""                                   "0.531805"       "0.709657"        "0.10183"  
%     ""                                   "0.0323393"      "0.0186795"       "0.087689" 
%     ""                                   "0.028235"       "0.0340118"       "0.075378" 
%     ""                                   "0.0772689"      "0.14137"         "0.25027"  
%     ""                                   "0.0765672"      "0.0661956"       "0.082955" 
%     ""                                   "0.154785"       "0.241823"        "0.10098"  
%     ""                                   "0.0424225"      "0.0421235"       "0.048766" 
%     ""                                   "0.0174002"      "0.0539344"       "0.080137" 
%     ""                                   "0.992778"       "1.647"           "0.13356"  
%     ""                                   "0.322439"       "0.0192129"       "0.14882"  
%     ""                                   "0.161729"       "0.024321"        "0.1671"   
%     ""                                   "0.000666965"    "0.00081132"      "0.20506"  
%     ""                                   "0.00335557"     "0.00119346"      "0.10299"  
%     ""                                   "0.0596057"      "0.0367815"       "0.1258"   
%     ""                                   "0.0868413"      "0.0587702"       "0.11181"  
%     ""                                   "0.137376"       "0.0645667"       "0.087444" 
%     ""                                   "0.24776"        "0.17852"         "0.13996"  
%     ""                                   "0.110702"       "0.105937"        "0.05523"  
%     ""                                   "0.0368951"      "0.0125853"       "0.15751"  
%     ""                                   "1.16411"        "2.38842"         "0.19046"  
%     ""                                   "0.0365157"      "0.0146401"       "0.10602"  
%     ""                                   "0.0322072"      "0.0496171"       "0.11352"  
%     ""                                   "0.0792962"      "0.0825911"       "0.093476" 
%     ""                                   "0.104091"       "-0.0798459"      "0.02826"  
%     ""                                   "3.16122"        "5.90486"         "0.15458"  
%     ""                                   "0.00622926"     "0.00628244"      "0.094175" 
%     ""                                   "0.0324675"      "-0.020394"       "0.058234" 
%     ""                                   "0.00753743"     "0.00612686"      "0.25716"  
%     ""                                   "0.0969285"      "0.355654"        "0.072763" 
%     ""                                   "0.467903"       "0.710097"        "0.075723" 
%     ""                                   "0.168985"       "0.505565"        "0.059482" 
%     ""                                   "1.201"          "-0.197046"       "0.078195" 
%     ""                                   "0.600793"       "-0.0689472"      "0.13908"  
%     ""                                   "0.344591"       "0.434219"        "0.097154" 
%     ""                                   "0.716374"       "-0.281739"       "0.15599"  
%     ""                                   "1.17047"        "0.0156384"       "0.16323"  
%     ""                                   "0.0439412"      "0.0310794"       "0.090127" 
%     ""                                   "0.0519911"      "0.00878296"      "0.09907"  
%     ""                                   "0.0926941"      "0.121849"        "0.055612" 
%     ""                                   "0.0926755"      "0.182479"        "0.038982" 
%     ""                                   "0.243258"       "0.55855"         "0.069416" 
%     ""                                   "0.0583391"      "0.137995"        "0.061849" 
%     ""                                   "0.150113"       "-0.101984"       "0.075207" 
%     ""                                   "0.420333"       "0.430673"        "0.20178"  
%     ""                                   "0.0635436"      "0.0513765"       "0.044315" 
%     ""                                   "0.285797"       "-0.121934"       "0.055856" 
%     ""                                   "0.0446408"      "-0.0243226"      "0.045466" 
