
%% Script clip to reformat exogeneous variables.

% input variable for this clip:
% parameter_customized,total_households,input_dir,output_dir,seed

seed_in_use = seed;
rng(seed_in_use);

% Directories
ITAM_input_dir = input_dir;
ITAM_output_dir = output_dir;
ITAM_mfile_dir = pwd;


theta = parameter_customized;

% Switch on or off for distance construction, true for triangular,false for
% two trip method
onetrip_dist = false;

%% import characteristics from the stata

% read tables
cd(ITAM_input_dir)
NTjchars2 = table2array(readtable('NTjchars2_stata.csv'));
NTJfirmnum2 = table2array(readtable('NTJfirmnum2_stata.csv'));
Hhchars2 = table2array(readtable('Hhchars2_stata.csv'));
cat_prices = table2array(readtable('cat_prices_stata.csv'));
cat_iv = cat_prices;

cd(ITAM_mfile_dir)

% test whether the identifiers of different datasets are aligned
test = zeros(4,1);
temp_ntindex = unique(NTjchars2(:,[1 2 3]),'stable','rows');
temp_nindex  = unique(NTJfirmnum2(:,1),'stable');
test(1) = sum(sum(temp_ntindex~=cat_prices(:,1:3)));
test(2) = sum(sum(cat_prices(:,1:3)~=cat_iv(:,1:3)));
test(3) = sum(sum(cat_iv(:,1:3)~=NTJfirmnum2(:,1:3)));
test(4) = sum(sum(temp_nindex~=Hhchars2(:,1)));      % only hhnumber in hhchars
if sum(test)~=0
    disp('Warning: not all input data arrays have the same ordering of observations')
    return
end

% Transform all data to be sort time id store (the original stata files are
% sorted in id time store (March 20))
NTjchars2 = sortrows(NTjchars2,[3 1 5]);
NTJfirmnum2 = sortrows(NTJfirmnum2,[3 1]);
Hhchars2 = sortrows(Hhchars2,1);
cat_prices = sortrows(cat_prices,[3 1]);
cat_iv = sortrows(cat_iv,[3 1]);

%% Set fundamentals

% Size of the total households allowed from Stata. In all the recent
% experiments this number is a constant, 7500.
N=size(Hhchars2,1);

T=size(NTJfirmnum2,1)/N; % Period sampled.
K=8; % Category.
J=30; % Choice set size,

C=J*(J+1)/2; % number of store-store combination
S=9; % number of firm. Currently, let it be 9, identical to their's code.
NT=N*T; % Total number of purchase in the sample
JK=J*K;
% reformat hhchars matrix
hhchars=repmat(Hhchars2,T,1);

%% Reformat matrices



% cat_comb, see table 3 of the paper, this is just to assume
% complementarity between categroies in the Lambda matrix
cat_comb = cell(4,1);
cat_comb{1} = [3 4];
cat_comb{2} = [8 2];
cat_comb{3} = [1 5 7]; % notice this means, 1,5,7 can mutally complementary, i.e. 15, 17 or 57
cat_comb{4} = 6;
% Number of interactive term, in table 3, this is part of Lambda estimate in model 2
Linteract = 0;
L = size(cat_comb,1);
for l=1:L
    c_cat_l = cat_comb{l};
    no_cat = size(c_cat_l,2);
    if no_cat>1
        combs = combnk(c_cat_l,2);
        Linteract = Linteract + size(combs,1);
    end
end



% hh code
hh_codeTN0 = cat_prices(:,1);   
hh_codeN0 = unique(hh_codeTN0,'stable');    % list of unique household codes
N0 = length(hh_codeN0);         % total number of households
hh_ix_1=1:N;
hh_ix_2=1:NT;



% Time indicators
weeks = cat_prices(:,2);
weeks = unique(weeks);
week78 = 78; % middle week of sample period (78/156)



week0 = cat_prices(:,2);
% observations that are in 'week 78' (middle week of sample period)
ix_w78 = find(week0 == week78); % pinpoint the all positions of week 78
period = cat_prices(:,3);
year = floor(week0/100);
week = week0-year*100;
qu1 = week<14;               % quarter dummy
qu2 = week>=14 & week<=26;
qu3 = week>=27 & week<=39;
qu4 = week>=40;
y2 = year==2;             % year dummy
y3 = year==3;
y4 = year==4;
quarter = [qu2 qu3 qu4 (qu4 & y2 | (qu1|qu2|qu3)& y3)  (qu4 & y3 | (qu1|qu2|qu3)& y4)];


% vectors predefined for expanding the matrix
on = ones(NT,1);                    % vectors predefined for expansion
oj = ones(J,1);
ok = ones(K,1);
o2 = ones(2,1);
os = ones(S,1);
oc = ones(C,1);

% Iseful indexing matrices for store pair
temp = true(J,J);                   % store / store pair indices
temp = triu(temp);
ix_JJ2C = false(J*J,1);
ix_JJ2C(temp) = true;
J1 = (1:J)';
storeindex = zeros(J*J,2);          % all combinations (j,j'), 1<=j,j'<=J
storeindex(:,1) = kron(ones(J,1),J1);
storeindex(:,2) = kron(J1,ones(J,1));
storeindex = storeindex(ix_JJ2C,:); % combinations (j,j') s.th 1<=j<=j'<=J
onestop = storeindex(:,1)==storeindex(:,2);

% household size and income in pounds
z = [sum(hhchars(:,[10 11]),2) hhchars(:,end)/10 quarter];
% number of demographic characteristics (including time dummies)
L2 = size(z,2);
z(:,2) = z(:,2)./z(:,1);            % change income to per capita income

% Household-time-store characteristics
ntjchars=NTjchars2;
x0 = ntjchars;                      % [Hhnumber week t Ncommon n salesarea easting northing dist]
ncommon = x0(:,4);                  % the first ncommon stores in i's choice set are there in t=1 and t=2
ncommon = ncommon(x0(:,5)==1,:);    % pick one row for each household (first store for each hh.)
x0(:,6) = log(x0(:,6));             % log of store size (salesarea)
x0(:,7) = x0(:,7)/1000;             % convert to metres. (dist = x0(:,9) is already in metres).
x0(:,8) = x0(:,8)/1000;
storechars = zeros(NT,J,4);
for i=1:N
    for t=1:T
        ix_temp = x0(:,1)==hh_codeN0(hh_ix_1(i)) & x0(:,3)==t;
        storechars(i+N*(t-1),:,:) = x0(ix_temp,6:end);      % logsalesarea easting northing dist
    end
end


x = storechars(:,:,1);
% number of store characteristics enterring into consumer's taste
L1 = size(x,3);
% Distance-related matrix construction
log_stsize = storechars(:,:,1);
east = storechars(:,:,2);
north = storechars(:,:,3);
dist = storechars(:,:,4);
L3 = 2;
dist = reshape(dist,NT,J,1);
dist = dist(:,:,oj);
dist_t = permute(dist,[1 3 2]);         % transpose
north = reshape(north,NT,J,1);
north = north(:,:,oj);
north_t = permute(north,[1 3 2]);       % transpose
east = reshape(east,NT,J,1);
east = east(:,:,oj);
east_t = permute(east,[1 3 2]);         % transpose
% Set two ways of distance construction

xp_triangular = zeros(NT,J,J,L3);       % distance home -> A -> B -> home
for i=1:NT
    % dummy for two stores
    xp_triangular(i,:,:,1) = ~eye(J);
    % triangular distance when store pair c = (A,B)
    distA = squeeze(dist(i,:,:));       % distance home to A
    distB = squeeze(dist_t(i,:,:));     % distance home to B
    % distance A to B.
    % (distAB)^2 = (eastingA-eastingB)^2 + (northingA-northingB)^2
    distAB = squeeze(sqrt( (east(i,:,:) - east_t(i,:,:)).^2  + (north(i,:,:) - north_t(i,:,:)).^2 ));
    % combined distance from home to A to B to home.
    dist_comb = distA + distAB + distB;
    xp_triangular(i,:,:,2) = ~eye(J).*dist_comb + eye(J).*(distA*2);
    % eye(J) gives the diagonal the (j,0) "pairs"
end

xp_twotrips = zeros(NT,J,J,L3);         % distance home -> A -> home -> B -> home
for i=1:NT
    % dummy for two stores
    xp_twotrips(i,:,:,1) = ~eye(J);
    % sum of distance when store pair c = (A,B)
    distA = squeeze(dist(i,:,:));       % distance home to A
    distB = squeeze(dist_t(i,:,:));     % distance home to B
    dist_sum = 2*(distA + distB);
    xp_twotrips(i,:,:,2) = ~eye(J).*dist_sum + eye(J).*(distA*2);
    % eye(J) gives the diagonal the (j,0) "pairs"
end

if onetrip_dist
    xp = xp_triangular;
else
    xp = xp_twotrips;
end
clear xp_triangular xp_twotrips


% Firm number mapping
firmnum = NTJfirmnum2(:,4:end); % NT*J, storing firm number
firmnum9 = zeros(NT,J);
% ASDA                      1
frmcl9{1} = 'Asda';
ix_temp = (firmnum==24);
firmnum9(ix_temp) = 1;
% DISC (ALDI, LIDL, NETTO)  2
frmcl9{2} = 'Discounter';
ix_temp = (firmnum==103);
firmnum9(ix_temp) = 2;
% ICELAND                   3
frmcl9{3} = 'Iceland';
ix_temp = (firmnum==12);
firmnum9(ix_temp) = 3;
% MORRISONS                 4
frmcl9{4} = 'Morrisons';
ix_temp = (firmnum==131);
firmnum9(ix_temp) = 4;
% MS                        5
frmcl9{5} = 'MS';
ix_temp =( firmnum==27);    % trad
firmnum9(ix_temp) = 5;
% OTHER (OTHER, BUDGEN, KWIKSAVE, SAFEWAY, COOP, SOMERFIELD) 6
frmcl9{6} = 'Other';
ix_temp = (firmnum==54);
firmnum9(ix_temp) = 6;
% SAINSBURY                 7
frmcl9{7} = 'Sainsbury';
ix_temp = firmnum==70;
firmnum9(ix_temp) = 7;
% TESCO                     8
frmcl9{8} = 'Tesco';
ix_temp = (firmnum==71);
firmnum9(ix_temp) = 8;
% WAITROSE                  9
frmcl9{9} = 'Waitrose';
ix_temp = firmnum==82;
firmnum9(ix_temp) = 9;

firmnum9(firmnum9==0) = 1;


% used later to add in firm category fixed effect in the utility
chain9 = zeros(NT,J,9);
for s=1:9
    chain9(:,:,s) = firmnum9==s;
end
% NxTxJxKxS. Expand for categories
chain9 = reshape(chain9,NT,J,1,9);
chain9 = chain9(:,:,ok,:);

% generate and save individual heterogeneity
nu = randn(N,J+1+1,K);
nu1 = nu(:,1:J,:);
nu2 = nu(:,J+1,:);
nu3 = randn(NT,1);
nupr = rand(N,1,1);
nupr = sqrt(-2*log(nupr)); % rayleigh distribution
nugamma = randn(N,1,1,2);
U = rand(N,5);
nu_const = randn(N,1);
% There is a whole lot process of generating a nu3_78 in the original
% script and it is a trivial term. See what I commented away below.
% It gives a random shock for those households that have week id
% equal to 78 the same shock received by them in nu3
nu3_78 = randn(N,1);
% hhcode6000 = hh_codeTN0(hh_ix_2);
% hh_78 = hhcode6000(ix_w78);
% for i=1:length(hh_78)
%     nu3_78(hh_78(i) == hh_code) = nu3(ix_w78(i));
% end

cd(ITAM_output_dir)
save nu_genrating_choice.mat nu1 nu2 nu3 nu3_78 nupr nugamma U nu_const


cd(ITAM_mfile_dir)

% Price
p0 = cat_prices(:,4:end);
P = zeros(NT,J,K);
for k=1:K
    temp2 = p0(:,k:K:J*K);
    P(:,:,k) = temp2;
end


%%  Store all data to inp, a struct

% This part is just reversing the data annoucement part of
% simulate_choice.m

% Read data from inp.. (NoComment means scalar)
inp.ix2=ix_JJ2C; % ix*J*J*2*C, upper tri of store pair matrix
inp.price=P; %  NT*J*K, price
inp.C=C; % store combination, C=30*29/2=435
inp.NT=N*T; % total sampled purchase, N=2000,T=3
inp.T=T; % periodsa
inp.J=J; % size of choice set
inp.K=K; % size of category set
inp.S=S; % # firms
inp.on=ones(NT,1); % N*1, on, oj, ok are all matrix repetitors
inp.oj=ones(J,1); % J*1
inp.ok=ones(K,1); % K*1
inp.x=x; % NT*J, storechars(:,:,1), i.e. store size,
inp.xp=xp; % NT*J*J*2, distance matrix, NT*J*J*2 = each i in t * store A * store B * onestore or two store trip
% the first dimension of xp is "whether visit two stores for that trip"
inp.z=z; % N*4, 4 chars of matrix hhchars
inp.L1=L1; % number of store chars enterring into taste mu
inp.L2=L2; % number of hh chars enterring into taste mu
inp.L3=L3; % number of maximum store allowed to visit, L3=2
inp.Li=Linteract; % number of interactive elasticity
inp.chain0 = squeeze(chain9(:,:,1,:));  
inp.cat_comb = cat_comb;
% storeindex is a C*2 matrix, where C=(J*(J-1))/2. It just lists every
% posible combination with the order as
% introduced in its first and second col. Its first col is (1-J)*J times
% and the second is 11..1,22..2, ... JJ...J, with the selection with the
% unique combination of stores.
inp.storeindex=storeindex; % C*2

% All the individual unobserved characteristic
inp.nupr=nupr; % N*1*1
inp.nu_const=nu_const ; % NT*1
inp.nugamma=reshape(nugamma,N,1,1,2); % NT*1*1*2, nugamma = randn(N,1,1,2)
inp.nu1=nu1; % see below
inp.nu2=nu2;
inp.nu3=nu3;


inp.discount = 0; % for use in counterfactuals
inp.fix = 'empty'; % for use when calculating conditional derivatives

