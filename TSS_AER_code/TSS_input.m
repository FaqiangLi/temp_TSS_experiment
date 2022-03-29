function [inp,inp2] = TSS_input(N,cat_comb,IV_comb,outofsample)

if outofsample              % always use N=2000 households for out-of-sample prediction
    N = 2000;
end
if N>2000
    disp('Monitor memory use - a large number of observations is used.')
end
%___________________Settings_______________________________________________
T = 3;                      % number of time periods
prelim_est = false;         % estimation performed to assign stores
first_time = false;         % first time using this sample
onetrip_dist = false;       % use distance of visiting both stores on one trip  
new_moments = true;         % include moment conditions involving correlations across time periods
cat_1st_instr = true;       % use 'store pair contains only 1 store' as an 
                            % instrument in the category moments.
price_instr = false;        % use price IVs (rather than price itself)

if first_time
    prompt = 'Are you sure you want to overwrite the saved estimation and validation samples? Y/N: ';
    str = input(prompt,'s');
    if str ~= 'Y'
        return
    end
end

%___________________Import data____________________________________________
%dir = 'OneSPC/';            % import data from data directory
%dir = 'REVISEDPRICE/';
dir = 'REVISED_DATA_JULY18/';
cat_USIs = table2array(readtable([dir 'cat_USIs.csv']));
cat_spends = table2array(readtable([dir 'cat_spends.csv']));
cat_prices = table2array(readtable([dir 'cat_prices.csv']));
if price_instr
    cat_iv = table2array(readtable([dir 'cat_iv.csv']));
else
    cat_iv = cat_prices;
end
hhchars = table2array(readtable([dir 'HHchars2.csv']));
ntjchars = table2array(readtable([dir 'NTJchars2.csv']));
ntjfirmnum = table2array(readtable([dir 'NTJfirmnum2.csv']));
store_visit = table2array(readtable([dir 'store_visit.csv']));
sameday = table2array(readtable([dir 'sameday.csv']));

hhchars = repmat(hhchars,T,1);      % hhchars are constant across time

%____________test__________________________________________________________
% check that observations are ordered in the same way in all data arrays:
test = zeros(7,1);
test(1) = sum(sum(cat_iv(:,1:3)~=cat_prices(:,1:3)));
test(2) = sum(sum(cat_prices(:,1:3)~=cat_spends(:,1:3)));
test(3) = sum(sum(cat_spends(:,1:3)~=cat_USIs(:,1:3)));
test(4) = sum(sum(cat_USIs(:,1:3)~=ntjfirmnum(:,1:3)));
test(5) = sum(sum(ntjfirmnum(:,1:3)~=store_visit(:,1:3)));
test(6) = sum(sum(sameday(:,1:3)~=store_visit(:,1:3)));
test(7) = sum(sum(sameday(:,1)~=hhchars(:,1)));      % only hhnumber in hhchars
if sum(test)~=0
    disp('Warning: not all data arrays have the same ordering of observations')
    %return
end
% check that there are not observations that are NaN
test = zeros(9,1);
test(1) = sum(isnan(cat_iv(:)));
test(2) = sum(isnan(cat_prices(:)));
test(3) = sum(isnan(cat_spends(:)));
test(4) = sum(isnan(cat_USIs(:)));
test(5) = sum(isnan(ntjfirmnum(:)));
test(6) = sum(isnan(ntjchars(:)));
test(7) = sum(isnan(hhchars(:)));
test(8) = sum(isnan(store_visit(:)));
test(9) = sum(isnan(sameday(:)));
if sum(test)~=0
    disp('Warning: some data array contains NaN entries.')
    %return
end
%____________end test______________________________________________________

% count number of interaction parameters
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

% for use in post-estimation calculations
weeks = cat_spends(:,2);
weeks = unique(weeks);
week78 = weeks(78); % middle week of sample period (78/156)

hh_codeTN0 = cat_spends(:,1);   % first column in cat_spends is househould number: household codes
hh_codeN0 = unique(hh_codeTN0,'stable');    % list of unique household codes
N0 = length(hh_codeN0);         % total number of households
J = 30;                         % stores in each person's choice set
C = J*(J+1)/2;                  % number of store pairs
K = 8;                          % categories of groceries
S = 9;                          % firms (with some aggregation)
NT = N*T;                       % no. of households * no. of time periods
% draw N random households
if first_time               % if new consumers need to be drawn
    % estimation sample
    hh_code = TSS_notin(hh_codeN0,[],N);
    hh_code = sort(hh_code,'ascend');
    save hh_code_estimation_sample hh_code
    % validation sample
    hh_code_insample = hh_code;
    hh_code = TSS_notin(hh_codeN0,hh_code_insample,N);
    hh_code = sort(hh_code,'ascend');
    save hh_code_validation_sample hh_code
    %____________test__________________________________________________________
    if length(unique(hh_code))~=N
        disp('Warning: the wrong number of consumers has been drawn')
    end
    count1 = 0;
    count2 = 0;
    for i=1:N
        if sum(hh_code(i)==hh_code_insample)~=0
            count1 = count1 + 1;
        end
        if sum(hh_code(i)==hh_code([1:i-1 i+1:end]))~=0
            count2 = count2 + 1;
        end
    end
    if count1>0
        disp('At least one household occurs both in estimation sample and validation sample')
    end
    if count2>0
        disp('At least one household occurs more than once in the validation sample')
    end
    %____________end test______________________________________________________
else
    if N==2000
        if outofsample
            load hh_code_validation_sample
        else
            load hh_code_estimation_sample
        end
    else
        load hh_code_estimation_sample
        hh_code = TSS_notin(hh_code,[],N);
        hh_code = sort(hh_code,'ascend');
    end
end
%____________test__________________________________________________________
if length(unique(hh_code))~=N
    disp('Warning: the wrong number of consumers have been drawn')
end
%____________end test______________________________________________________

% logical index for location of sampled households:
hh_lix_1 = zeros(N0,1);          % within hh_codeN0
hh_lix_2 = zeros(2*N0,1);        % within hh_codeTN0
for i=1:N
    location1 = find(hh_codeN0==hh_code(i));
    location2 = find(hh_codeTN0==hh_code(i));
    hh_lix_1(location1) = true;
    hh_lix_2(location2) = true;
end
% regular (not logical) index
hh_ix_1 = find(hh_lix_1);       % hh_code = hh_codeN0(hh_ix_1)
hh_ix_2 = find(hh_lix_2);
%____________test__________________________________________________________
if ~isequal(hh_code,hh_codeN0(hh_ix_1)) ...
        | ~isequal(sum(hh_lix_1),N) ...
        | ~isequal(sum(hh_lix_2),NT) ...
        | ~isequal(repmat(hh_code,T,1),hh_codeTN0(hh_ix_2))
    disp('Warning: problem with hh_ix_1.')
end
%____________end test______________________________________________________

% index to access sampled households from ntjchars
ntjch_ix = [];
for i=1:NT
    endpoint = hh_ix_2(i)*J;            % there are J rows for each household / time period
    range = endpoint-J+1:endpoint;
    ntjch_ix = [ntjch_ix; range'];
end

% floor area, distance [30 nearest stores for each consumer]
if prelim_est
    ntjchars2 = ntjchars;
end
ntjchars = ntjchars(ntjch_ix,:); 
%____________test__________________________________________________________
for j=1:J
    ix_j = j:J:NT*J; % j-th store for each household/time period
    if ~isequal(ntjchars(ix_j,1),hh_codeTN0(hh_ix_2))
        disp('Problem with ntjchars.')
    end
end
%____________end test______________________________________________________
% chain number (1-16) first store ... 30th store, for each consumer
ntjfirmnum = ntjfirmnum(hh_ix_2,:); 
hhchars = hhchars(hh_ix_2,:); 
% prices in each of 8 categories for each of 30 stores, for each consumer
cat_prices = cat_prices(hh_ix_2,:);
% instrumental variables for price (same format as cat_prices)
cat_iv = cat_iv(hh_ix_2,:);
% expenditure in each of 8 categories for each of 30 stores, for each
% consumer
cat_spends(:,4:end) = 1e-1*cat_spends(:,4:end); % in units of 10 pounds per week
if prelim_est
    cat_spends2 = cat_spends;       % all observations - for use in estimation of simple choice model
end
cat_spends = cat_spends(hh_ix_2,:); % observations used in estimation of main model
if prelim_est
    cat_USIs2 = cat_USIs;           % all observations - for use in estimation of simple choice model
end
cat_USIs = cat_USIs(hh_ix_2,:);
store_visit = store_visit(hh_ix_2,4:end);
sameday = sameday(hh_ix_2,4:end);
number_stores = sameday(:,1);       % the number of stores visited (whether or not main store for any category)
same_day = sameday(:,2);            % 1 if the consumer visited two stores on the same day. 0 if not.            

week0 = cat_prices(:,2);
% observations that are in 'week 78' (middle week of sample period)
ix_w78 = find(week0 == week78);
period = cat_prices(:,3);
year = floor(week0/100);
week = week0-year*100;
qu1 = week<14;
qu2 = week>=14 & week<=26;  
qu3 = week>=27 & week<=39;
qu4 = week>=40;
y2 = year==2003;
y3 = year==2004;
y4 = year==2005;
quarter = [qu2 qu3 qu4 (qu4 & y2 | (qu1|qu2|qu3)& y3) ...
    (qu4 & y3 | (qu1|qu2|qu3)& y4)];

on = ones(NT,1);                    % vectors predefined for expansion
oj = ones(J,1);
ok = ones(K,1);
o2 = ones(2,1);
os = ones(S,1);
oc = ones(C,1);

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


% demographic variables: adults+children, income, time dummies
% income divided by 10, so it is in tens of pounds, like expenditure
z = [sum(hhchars(:,[10 11]),2) hhchars(:,end)/10 quarter];
% number of demographic characteristics (including time dummies)
L2 = size(z,2);
% labels for household characteristics with estimated coefficients 
hhchcl = cell(L2-1,2);
hhchcl{1,1} = 'Household characteristics';
hhchcl{1,2} = 'household size';
hhchcl{2,2} = 'Q2';
hhchcl{3,2} = 'Q3';
hhchcl{4,2} = 'Q4';
hhchcl{5,2} = 'Y2004';
hhchcl{6,2} = 'Y2005';

z_NTL2 = z;                         % NT x L2
z(:,2) = z(:,2)./z(:,1);            % change income to per capita income

x0 = ntjchars;                      % [Hhnumber week t Ncommon n salesarea easting northing dist]
ncommon = x0(:,4);                  % the first ncommon stores in i's choice set are there in t=1 and t=2
ncommon = ncommon(x0(:,5)==1,:);    % pick one row for each household (first store for each hh.)
x0(:,6) = log(x0(:,6));             % log of store size (salesarea)
x0(:,7) = x0(:,7)/1000;             % convert to metres. (dist = x0(:,9) is already in metres).
x0(:,8) = x0(:,8)/1000;
% log of storesize 
if prelim_est
    x02 = ntjchars2;
    x02(:,6) = log(x02(:,6));
end
storechars = zeros(NT,J,4);
for i=1:N
    for t=1:T
        ix_temp = x0(:,1)==hh_codeN0(hh_ix_1(i)) & x0(:,3)==t;
        storechars(i+N*(t-1),:,:) = x0(ix_temp,6:end);      % logsalesarea easting northing dist
    end
end

if prelim_est
    storechars2 = zeros(2*N0,J,2);                          % for preliminary estimation
    for i=1:N0
        for t=1:T
            ix_temp = x02(:,1)==hh_codeN0(i) & x02(:,3)==t;
            storechars2(i+N0*(t-1),:,:) = x02(ix_temp,[6 9]);   % logsalesarea dist
        end
    end
    clear x02
end

% store characteristics which enter first-order term of continuous utility
x = storechars(:,:,1);
% number of such characteristics
L1 = size(x,3);
stchcl = cell(L1,2); % labels for store characteristics
stchcl{1,1} = 'Store characteristics';
stchcl{1,2} = 'log store size';

% NxJxKxJxL3. store-PAIR characteristics, for discrete utility
% number of such characteristics [two-store or not, combined distance]
% [log_salesarea easting northing dist]
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

spcl = cell(4,2);
spcl{1,1} = 'Store pair fixed coefficients';
spcl{1,2} = 'choice contains two stores (0/1)';
spcl{2,2} = 'sum of distances';
spcl{3,1} = 'Store pair random coefficients';
spcl{3,2} = 'choice contains two stores (0/1)';
spcl{4,2} = 'sum of distances';

% NxJ. Entry i,j gives the number (from 1 to S) of the chain to which
% consumer i's j-th store belongs.
firmnum = ntjfirmnum(:,4:end);
frmcl = cell(16,1);         % string labels for firms (16 different firms, including one 'other' category)
% ALDI		1	DISC
frmcl{1} = 'Aldi';
% ASDA		2
frmcl{2} = 'Asda';
% BUDGEN	3	OTHER
frmcl{3} = 'Budgen';
% COOP		4	
frmcl{4} = 'Coop';
% ICELAND	5
frmcl{5} = 'Iceland';
% KWIKSAVE	6	OTHER
frmcl{6} = 'Kwiksave';
% LIDL		7	DISC
frmcl{7} = 'Lidl';
% MORRISONS	8
frmcl{8} = 'Morissons';
% MS		9
frmcl{9} = 'MS';
% NETTO		10	DISC
frmcl{10}  = 'Netto';
% OTHER		11
frmcl{11} = 'Other';
% SAFEWAY   12  OTHER
frmcl{12} = 'Safeway';
% SAINSBURY	13	
frmcl{13} = 'Sainsbury';
% SOMERFIELD14
frmcl{14} = 'Somerfield';
% TESCO	    15
frmcl{15} = 'Tesco';
% WAITROSE	16
frmcl{16} = 'Waitrose';

%  Aggregating to nine firms
frmcl9 = cell(9,1);         % string labels for firms aggregated to 9
firmnum9 = zeros(NT,J);     % array of firm codes
% ASDA                      1
frmcl9{1} = 'Asda';
ix_temp = firmnum==2;
firmnum9(ix_temp) = 1;
% DISC (ALDI, LIDL, NETTO)  2	
frmcl9{2} = 'Discounter';
ix_temp = (firmnum==1|firmnum==7|firmnum==10);
firmnum9(ix_temp) = 2;
% ICELAND                   3
frmcl9{3} = 'Iceland';
ix_temp = firmnum==5;
firmnum9(ix_temp) = 3;
% MORRISONS                 4
frmcl9{4} = 'Morrisons';
ix_temp = firmnum==8;
firmnum9(ix_temp) = 4;
% MS                        5
frmcl9{5} = 'MS';
ix_temp = firmnum==9;
firmnum9(ix_temp) = 5;
% OTHER (OTHER, BUDGEN, KWIKSAVE, SAFEWAY, COOP, SOMERFIELD) 6
frmcl9{6} = 'Other';
ix_temp = (firmnum==3|firmnum==4|firmnum==6|firmnum==11|firmnum==12|firmnum==14);
firmnum9(ix_temp) = 6;
% SAINSBURY                 7
frmcl9{7} = 'Sainsbury';
ix_temp = firmnum==13;
firmnum9(ix_temp) = 7;
% TESCO                     8
frmcl9{8} = 'Tesco';
ix_temp = firmnum==15;
firmnum9(ix_temp) = 8;
% WAITROSE                  9
frmcl9{9} = 'Waitrose';
ix_temp = firmnum==16;
firmnum9(ix_temp) = 9;

% 1. Bakery
% 2. Dairy
% 3. Drink
% 4. Dry
% 5. FV
% 6. HH
% 7. Meat
% 8. Milk
catcl = cell(8,1);
catcl{1} = 'Bakery';
catcl{2} = 'Dairy';
catcl{3} = 'Drink';
catcl{4} = 'Dry grocery';
catcl{5} = 'Fruit & vegetables';
catcl{6} = 'Household goods';
catcl{7} = 'Meat';
catcl{8} = 'Milk';

% NTxJxS. Logical array to pick stores that are in a particular chain. 
chain9 = zeros(NT,J,9);
for s=1:9
    chain9(:,:,s) = firmnum9==s;
end
% NxTxJxKxS. Expand for categories
chain9 = reshape(chain9,NT,J,1,9);
chain9 = chain9(:,:,ok,:);

% same with chains not aggregated
% NTxJxS. Logical array to pick stores that are in a particular chain. 
chain16 = zeros(NT,J,16);
for s=1:16
    chain16(:,:,s) = firmnum==s;
end
% NxTxJxKxS. Expand for categories
chain16 = reshape(chain16,NT,J,1,16);
chain16 = chain16(:,:,ok,:);

% random draws saved so the same draws can be used for different runs. 
if N==2000
    if first_time
        % For estimation sample
        nu = randn(N,J+1+1,K);
        nu1 = nu(:,1:J,:);
        nu2 = nu(:,J+1,:);
        nu3 = randn(NT,1);
        nupr = rand(N,1,1);
        nupr = sqrt(-2*log(nupr)); % rayleigh distribution
        nugamma = randn(N,1,1,2);
        U = rand(N,5);
        nu_const = randn(N,1);
        
        nu3_78 = randn(N,1);
        hhcode6000 = hh_codeTN0(hh_ix_2);
        hh_78 = hhcode6000(ix_w78);
        for i=1:length(hh_78)
            nu3_78(hh_78(i) == hh_code) = nu3(ix_w78(i));
        end
        
        
        
        save nu2000_estimation_sample nu1 nu2 nu3 nu3_78 nupr nugamma U nu_const
        
        % For validation sample
        nu = randn(N,J+1+1,K);
        nu1 = nu(:,1:J,:);
        nu2 = nu(:,J+1,:);
        nu3 = randn(NT,1);
        nupr = rand(N,1,1);
        nupr = sqrt(-2*log(nupr)); % rayleigh distribution
        nugamma = randn(N,1,1,2);
        U = rand(N,5);
        nu_const = randn(N,1);
        save nu2000_validation_sample nu1 nu2 nu3 nupr nugamma U nu_const
    else
        if ~outofsample
            load nu2000_estimation_sample
        else
            load nu2000_validation_sample
        end
    end
else
    % Random draws (one per household)
    nu = randn(N,J+1+1,K);
    nu1 = nu(:,1:J,:);
    nu2 = nu(:,J+1,:);
    nu3 = randn(NT,1);
    nupr = rand(N,1,1);
    nupr = sqrt(-2*log(nupr)); % rayleigh distribution
    nugamma = randn(N,1,1,2);
    U = rand(N,5);
    nu_const = randn(N,1);
    save nu_testing nu1 nu2 nu3 nupr nugamma U nu_const
end

rncl = cell(4,2); % labels
rncl{1,1} = 'Random terms';
rncl{1,2} = 'store/cat. specific';
rncl{2,2} = 'category specific';
rncl{3,2} = 'time-varying';
rncl{4,2} = 'fixed across cat./store';

% expenditures 
q0 = cat_spends(:,4:end);
usi0 = cat_USIs(:,4:end);
% place expenditures in an NxJxK array
Q0 = zeros(NT,J,K);
USI = zeros(NT,J,K);
for k=1:K
    Q0(:,:,k) = q0(:,k:K:J*K);
    USI(:,:,k) = usi0(:,k:K:J*K);
end

if prelim_est
    % For estimation of simple discrete choice model below:
    q02 = cat_spends2(:,4:end);
    usi02 = cat_USIs2(:,4:end);
    % place expenditures in an N0xJxK array
    Q02 = zeros(T*N0,J,K);
    USI2 = zeros(T*N0,J,K);
    for k=1:K
        Q02(:,:,k) = q02(:,k:K:J*K);
        USI2(:,:,k) = usi02(:,k:K:J*K);
    end
    clear q02 usi02
end

% prices 
p0 = cat_prices(:,4:end);
P = zeros(NT,J,K);
for k=1:K
    temp2 = p0(:,k:K:J*K);
    P(:,:,k) = temp2;
end
% instrumental variable for price
iv0 = cat_iv(:,4:end);
IV = zeros(NT,J,K);
for k=1:K
    IV(:,:,k) = iv0(:,k:K:J*K);
end

if prelim_est
    %______Compute probabilities based on store size and distance _____________
    ix = true(N0,1);  % create logical index for households where store of purchase *is* determined
    for i=1:2*N0
        usi_i = USI2(i,:,1);                        % USIs (for cat 1, but same for all categories)
        unique_usi = nonzeros(unique(usi_i));       % list of unique nonzero USIs
        number_usi = length(unique_usi);            % number of unique USIs
        if number_usi>2                             % [test]
            disp('Warning: there are three or more different USIs for a given consumer.')
            return
        end
        for l=1:number_usi;                         % number_usi = 1 or 2
            usi_l = usi_i==unique_usi(l);
            occurrences_usi = sum(usi_l);           % number of times the USI occurs
            if occurrences_usi>1                    % USI not matched uniquely to store
                ix(i) = false;
            end
        end
    end
    clear USI2
    
    Q0t = sum(Q02,3);                               % total expenditure by store
    clear Q02
    [temp,I0] = max(Q0t,[],2);                      % pick the store with biggest expenditure
    % hist(I0)                                      % to see distribution of distance rank of first store:
    
    % create a 0-1 array with a 1 only for the largest-expenditure store for
    % each household. This is the dependent variable in a simple discrete choice
    % model used to determine probabilities based on storesize and distance
    I = zeros(2*N0,J);
    for i=1:2*N0
        I(i,I0(i)) = 1;
    end
    
    % Estimate parameters of simple logit model predicting main store with
    % distance and log of storesize.
    options = optimoptions('fminunc','Display','iter','Algorithm','quasi-newton');
    simpleobj = @(beta)-TSS_simplechoicelikelihood(I(ix,:),storechars2(ix,:,1),storechars2(ix,:,2),beta);
    startval = ones(2,1);
    beta = fminunc(simpleobj,startval,options);
    %beta = fminunc(simpleobj,startval);
    clear storechars2
    save beta beta
end
load beta

%___________________ Assign purchases to a single store ___________________
% for some chains, the data do not tell us which of a firm's stores is
% visited. So far the entire purchased quantity is entered for every
% possible store.
% The following lines randomly assign the purchased quantity to one of the
% possible stores for each consumer and category where this happens
% (according to the probabilities estimated with the simple discrete choice
% model above).

count_zeros = 0;
for i0=1:N
    assigned = [];                                  % a given USI will be entered here when a specific store has been assigned for it
    chosen_store = [];                              % the number (from 1 to 30) of the store chosen will be entered here
    l = 1;                                          % counter: the l-th store chosen by i is the next to be assigned.
    for t=1:T
        for k=1:K
            i = (t-1)*N+i0;
            usi_i = USI(i,:,k);                         % USIs for cat k
            unique_usi = nonzeros(unique(usi_i));       % list of unique nonzero USIs for cat k
            number_usi = length(unique_usi);            % number of unique USIs for cat k. Should be 1.
            if number_usi>1
                disp('Warning: there are two or more different USIs within a category.')
                return
            elseif number_usi==1                        % if there is no USI for cat k, the category is not purchased
                usi = usi_i==unique_usi;                % (1x30) logical array for occurrences of the USI
                occurrences_usi = sum(usi);             % number of times the USI occurs
                usi = find(usi);                        % indices (from 1 to 30) for the occurrences of the USI
                if occurrences_usi>1                    % USI not matched uniquely to store
                    if ~ismember(unique_usi,assigned)   % if the USI has not yet been assigned to a specific store
                        % A different USI belonging to the same firm could already have been assigned to a
                        % specific store, therefore remove such stores from the set of possible stores to
                        % assign the current USI to. Also remove stores
                        % that do not occur in both time periods for the
                        % household, provided that does not leave out all
                        % stores with this usi. 
                        usi_new = usi;
                        for r=1:length(usi)
                            if ismember(usi_new(r),chosen_store) | ...
                                    (usi_new(r)>ncommon(i0) & min(usi_new)<ncommon(i0) )
                                usi_new(r) = 0;
                            end
                        end
                        usi_new = nonzeros(usi_new);
                        if isempty(usi_new)
                            usi_new = usi(1);
                        end
                        % Unconditional probabilities of visiting each store
                        prob = TSS_simplechoiceprob(storechars(i,:,1),storechars(i,:,2),beta);
                        % Pick the probabilities only for those stores with the
                        % given USI *and* which have not already been
                        % chosen.
                        prob = prob(usi_new);
                        % Conditional probability of each store among those
                        % with the same USI
                        prob = prob/sum(prob);
                        % Pick one of the stores with the given USI randomly, with
                        % probabilities given by prob
                        u = U(i0,l);
                        cumprob = cumsum(prob);
                        n=1;
                        chosen = 0;
                        while n<=length(usi_new) & chosen==0
                            if u<cumprob(n)
                                chosen = usi_new(n);
                                l = l+1;
                                assigned = [assigned; unique_usi];
                                chosen_store = [chosen_store; chosen];
                            end
                            n = n+1;
                        end
                    else                                    % find the index of the store the USI has been assigned to
                        ix = find(unique_usi==assigned);
                        chosen = chosen_store(ix);
                    end
                    % delete quantities except in chosen store
                    if chosen~=0
                        quantity_temp = Q0(i,chosen,k);
                        % set all to zero within the given USI
                        Q0(i,usi,k) = 0;
                        % except the chosen one
                        Q0(i,chosen,k) = quantity_temp;
                    else
                        Q0(i,usi,k) = 0;
                        count_zeros = count_zeros + 1;
                    end
                end
            end
        end
    end
end
% [test]
% maximum number of strictly positive expenditures in a given category
% should be 1:
for i=1:NT
    for k=1:K
        if sum(Q0(i,:,k)>0)>1
            disp('Warning: there are consumers with positive purchases in two or more stores for a given category.')
            return
        end
    end
end

%%%%%%% quantities in real terms (divided by prices) %%%%%%%%%%%%%%%%%%%%%%%%
Q = Q0./P;

% which store pair did each consumer choose?
storepairvisited = zeros(NT,C);
for i=1:NT
    Qi = squeeze(Q0(i,:,:));                % quantity by store and category
    Qipos = find(sum(Qi,2));                % indicator for whether store is used for any category
    if length(Qipos)~=1 & length(Qipos)~=2  % test: either one or two stores are visited
        disp('Warning: wrong number of stores visited.')
        %return
    end
    if length(Qipos)<number_stores(i)       % a store is visited which is not the main store for any category
        if number_stores(i)~=2              % test: this line can only be reached if number_stores(i)=2 
            disp('Warning: wrong number of stores visited.')
            %return
        end
        store_visit_i = store_visit(i,:);   % NTxJ
        if store_visit_i(Qipos)~=1          % test: make sure assigned store is recorded as visited
            disp('Warning: wrong number of stores visited.')
            %return
        end
        USI_ik = squeeze(USI(i,:,:));       % find stores visited which have a different USI than store Qipos
        USI_i = sum(USI_ik,2)~=0;
        store_visit_i(USI_i) = 0;           % set stores with same USI as the one in Qipos to zero
        store_visit_i = find(store_visit_i); 
        if length(store_visit_i)==0         % number_stores(i)=2, but the second store is not among J closest stores
            store_visit_i = Qipos;       
        elseif length(store_visit_i)>1      % second visited store is undetermined - must be assigned
            prob = TSS_simplechoiceprob(storechars(i,:,1),storechars(i,:,2),beta);
            prob = prob(store_visit_i);
            % Conditional probability of each store among those
            % with the same USI
            prob = prob/sum(prob);
            % Pick one of the stores with the given USI randomly, with
            % probabilities given by prob
            u = U(i-N*floor(i/N),5);
            cumprob = cumsum(prob);
            n=1;
            while n<=length(store_visit_i) & length(store_visit_i)~=1
                if u<cumprob(n)
                    store_visit_i = store_visit_i(n);
                end
                n = n+1;
            end
        end
        visited = sort([Qipos store_visit_i],'ascend');    % put store with lowest index first
    elseif length(Qipos)==1                 % one-stop shopper
        visited = [Qipos Qipos];
    elseif length(Qipos)==2                 % two-stop shopper
        visited = Qipos;
    end
    if length(visited)<2
        test = 1;
    end
    storepairvisited(i,:) = visited(1)==storeindex(:,1) & visited(2)==storeindex(:,2);
    if sum(storepairvisited(i,:),2)~=1
        display('Warning: store pair visited not correctly determined.')
        return
    end
end

% ________________________ Print descriptive statistics ___________________
disp('________________________________________________________')
disp('________________________________________________________')
disp(['Households:                                     ' num2str(N)])
disp(['Observations:                                   ' num2str(NT)])
no_2ss = sum(sum(storepairvisited(:,~onestop),2),1);
disp(['2-store visits:                                 ' num2str(no_2ss) ])
disp(['2-store visits, %:                                ' num2str(round(100*no_2ss/NT,1))])
vis11 = zeros(N,1);
vis22 = zeros(N,1);
vis12 = zeros(N,1);
for i=1:N
    vis11(i,:) = sum(storepairvisited(i,onestop),2)==1 & sum(storepairvisited(i+N,onestop),2)==1;
    vis22(i,:) = sum(storepairvisited(i,~onestop),2)==1 & sum(storepairvisited(i+N,~onestop),2)==1;
    vis12(i,:) = 1-vis11(i,:)-vis22(i,:);
end
xp_temp = reshape(xp,NT,J*J,2);
xp_temp = xp_temp(:,ix_JJ2C,:);
dist_temp0 = sum(xp_temp(:,:,2).*storepairvisited,2);
dist_temp = mean(dist_temp0,1);
disp(['Mean distance, all, km:                           '  num2str(round(dist_temp,1))])    
onestop_hh = sum(storepairvisited(:,onestop),2)>0;
dist_temp_1SS = mean(dist_temp0(onestop_hh,:),1);
disp(['Mean distance, 1-stop shoppers, km:               '  num2str(round(dist_temp_1SS,1))])   
dist_temp_2SS = mean(dist_temp0(~onestop_hh,:),1);
disp(['Mean distance, 2-stop shoppers, km:               '  num2str(round(dist_temp_2SS,1))])   
clear xp_temp dist_temp
spend_temp0 = 10*sum(sum(Q0,3),2);           % spending in pounds is divided by 10 in line 141.
spend_temp_all = mean(spend_temp0);
spend_temp_1SS = mean(spend_temp0(onestop_hh));
spend_temp_2SS = mean(spend_temp0(~onestop_hh));
disp(['Mean spending, all, pounds:                       ' num2str(round(spend_temp_all,1))])
disp(['Mean spending, 1-stop shoppers, pounds:           ' num2str(round(spend_temp_1SS,1))])
disp(['Mean spending, 2-stop shoppers, pounds:           ' num2str(round(spend_temp_2SS,1))])
for k=1:K
    catname = catcl{k};
    disp(['Mean spending on ' catname(1:3) ', pounds:                      ' num2str(round(10*mean(sum(Q0(:,:,k),2)),1))])
end
clear spend_temp0 spend_temp_all spend_temp_1SSS spend_temp_2SS
disp(['Mean income per cap., all, pounds:               ' num2str(round(10*mean(z(:,2)),1))])  % income divided by 10 in line 198
disp(['Mean income per cap., 1-stop shoppers, pounds:   ' num2str(round(10*mean(z(onestop_hh,2)),1))])  % income divided by 10 in line 198
disp(['Mean income per cap., 2-stop shoppers, pounds:   ' num2str(round(10*mean(z(~onestop_hh,2)),1))])  % income divided by 10 in line 198
disp(['Mean number of adults, all:                        ' num2str(round(mean(hhchars(:,10)),2))])
% disp(['Mean number of adults, 1-stop shoppers:   ' num2str(round(mean(hhchars(onestop_hh,10)),2))])
% disp(['Mean number of adults, 2-stop shoppers:   ' num2str(round(mean(hhchars(~onestop_hh,10)),2))])
disp(['Mean number of children, all:                      ' num2str(round(mean(hhchars(:,11)),2))])
% disp(['Mean number of children, 1-stop shoppers: ' num2str(round(mean(hhchars(onestop_hh,11)),2))])
% disp(['Mean number of children, 2-stop shoppers: ' num2str(round(mean(hhchars(~onestop_hh,11)),2))])
disp(['Mean age, all:                                    ' num2str(round(mean(hhchars(:,7)),1))])
% disp(['Mean age, 1-stop shoppers:               ' num2str(round(mean(hhchars(onestop_hh,7)),1))])
% disp(['Mean age, 2-stop shoppers:               ' num2str(round(mean(hhchars(~onestop_hh,7)),1))])
disp('________________________________________________________')


%________________ DEPENDENT VARIABLES _____________________________________
positiveamount = zeros(NT,C,2,K);                   % (a) Indicator for positive amount or not
quantity = zeros(NT,C,2,K);                         % (b) Quantity purchased
for i=1:NT
    c_i = find(storepairvisited(i,:));              % index (between 1 and C) of the store pair visited
    for k=1:K
        store_k = find(squeeze(Q(i,:,k)));          % index (between 1 and J) of the store used for k
        Q_k = Q(i,store_k,k);                       % quantity purchased of category k
        if store_k == storeindex(c_i,1)             % if store_k is the first store in store pair c
            positiveamount(i,c_i,1,k) = 1;
            quantity(i,c_i,1,k) = Q_k;
        elseif store_k == storeindex(c_i,2)         % if store_k is the second store in store pair c
            positiveamount(i,c_i,2,k) = 1;          % (when storeindex(c_i,1)=storeindex(c_i,2) [one-store 'store pair'],
            quantity(i,c_i,2,k) = Q_k;              % the quantity is placed in column 1 of 'quantity')
        end
    end
end


%______________ EXPLANATORY VARIABLES______________________________________
% (a) Price instrument
priceinstrument_cross_cat = zeros(NT,C,2,K,2);
priceinstrument = zeros(NT,C,2,K);
for i=1:NT
    for k=1:K
        priceinstrument(i,:,1,k) = IV(i,storeindex(:,1),k);
        priceinstrument(i,:,2,k) = IV(i,storeindex(:,2),k);
        % use the price IV for other products as additional price IV's
        comb_k = IV_comb{k};
        Lk = length(comb_k);
        for l=1:Lk
            IV_temp1 = IV(i,storeindex(:,1),comb_k(l));
            IV_temp2 = IV(i,storeindex(:,2),comb_k(l));
            IV_temp_mean = (IV_temp1 + IV_temp2 )/2;
            priceinstrument_cross_cat(i,:,1,k,l) = IV_temp_mean;
            priceinstrument_cross_cat(i,:,2,k,l) = IV_temp_mean;
        end
    end
end

% (b) Chain number associated with each store
% (b.i) Aggregated to 9 different chains
chainnumber9 = zeros(NT,C,2,K);
for i=1:NT
    chainnumber9(i,:,1,1) = firmnum9(i,storeindex(:,1));
    chainnumber9(i,:,2,1) = firmnum9(i,storeindex(:,2));
    chainnumber9(i,:,:,:) = chainnumber9(i,:,:,ok);
end
% (b.ii) 16 different chains
chainnumber16 = zeros(NT,C,2,K);
for i=1:NT
    chainnumber16(i,:,1,1) = firmnum(i,storeindex(:,1));
    chainnumber16(i,:,2,1) = firmnum(i,storeindex(:,2));
    chainnumber16(i,:,:,:) = chainnumber16(i,:,:,ok);
end
% (b.iii) chaindummies
chaindummies     = zeros(NT,C,2,K,9);
for s=1:9
        chaindummies(:,:,:,:,s) = chainnumber9==s;
end

% (c) Store characteristics
% there is one store characteristic: store size 
storecharacteristics = zeros(NT,C,2);
for i=1:NT
    storecharacteristics(i,:,1) = log_stsize(i,storeindex(:,1)); % store size
    storecharacteristics(i,:,2) = log_stsize(i,storeindex(:,2)); % store size
end
% NxCx2xK
storecharacteristics = storecharacteristics(:,:,:,ok);

% (d) Household characteristics
% household size (adults + children) , income , 5 time dummies ; in all, 7
% household characteristics
% N x C x 2 x K x 7
householdcharacteristics = reshape(z_NTL2,NT,1,1,1,L2);
householdcharacteristics = householdcharacteristics(:,oc,o2,ok,:);

priceinstruments = cat(5,priceinstrument,priceinstrument_cross_cat);

clear x0 z_NTL2 y2 y3 y4 U cat_iv cat_prices cat_spends cat_USIs dist dist_t...
    east east_t north north_t ...
    hhchars IV iv0 ix_temp ntjch_ix ntjchars ntjfirmnum priceinstrument ...
    priceinstrument_cross_cat Q q0 Q0 p0 quarter store_visit storechars ...
    USI usi0 firmnum firmnum9 hh_codeN0 hh_code hh_codeTN0 hh_ix_1...
    hh_ix_2 hh_lix_21...
    hh_lix_2 log_stsize nu qu1 qu2 qu3 qu4 same_day sameday temp2 temp week ...
    week0 year dist_comb dist_i dist_sum dist_t dist

%_____________________INSTRUMENT ARRAY ____________________________________
% NxCx2 [code (1-9) for each store]
ACi = chainnumber9(:,:,:,1); % 4th dimension contains 8 identical pages.
clear chainnumber9
% NxCx2
ACi2 = chainnumber16(:,:,:,1);
clear chainnumber16

% Concatenation done in multiple steps to limit memory demands:
A = cat(5,householdcharacteristics(:,:,:,:,[1 3:7]),priceinstruments);
clear householdcharacteristics priceinstruments
A = cat(5,A,storecharacteristics);
clear storecharacteristics
A = cat(5,A,chaindummies(:,:,:,:,1:end-1));
clear chaindummies
A = cat(5,A,ones(NT,C,2,K));
if cat_1st_instr
    onestop_expanded = reshape(onestop,1,C);
    A = cat(5,A,onestop_expanded(on,:,o2,ok));
end

% Moment conditions corresponding to cross-category price instruments for
% categories that are not interacted will be dropped. Create index for
% moments which will not be dropped:
instr_ix1 = true(size(A,5),K);
instr_ix2 = true(size(A,5)-5,K);
for k=1:K
    number_iv = length(IV_comb{k});
    instr_ix1(6+1+number_iv+1:6+3,k) = false;         % 6 hh chars, 3 price instr - the last two of which are cross-cat
    instr_ix2(1+1+number_iv+1:1+3,k) = false;         % 1 hh char , 3 price instr - the last two of which are cross-cat
end
L_Asp = 6;  % number of instruments for store pair choice
instr_ix3 = true(L_Asp,1);
instr_ix = [instr_ix1(:); instr_ix2(:); instr_ix3(:)];

L = size(A,5);
olA = ones(L,1);

% Nx6 [household size (adults + children), 5 time dummies]
AH = squeeze(A(:,1,1,1,1:6));
% NxCx2xKx2 [pricinstr,cross cat price IV]
AP = A(:,:,:,:,7:9);
% NxCx2x1 [storesize]
AS = squeeze(A(:,:,:,1,10));
% Cx1 [dummy for whether one-stop 'pair' or not]
A1i = onestop;
onestore = ~xp(:,:,:,1);            % 'one-store' pair (*not* two-store)
onestore = onestore(:,ix_JJ2C);
distance = xp(:,:,:,2);
distance = distance(:,ix_JJ2C);
% NxCx3
mn_priceIV = squeeze(mean(mean(AP(:,:,:,:,1),3),4));
per_cap_inc = z(:,2);
Asp = cat(3,distance,onestore,distance.^2,onestore.*distance,mn_priceIV,mn_priceIV./per_cap_inc(:,oc));
% test
if size(Asp,3) ~= L_Asp
   disp('Fix instr_ix3.') 
end
AP = reshape(AP,NT*C*2,K,3);
LP = size(AP,3) + size(AS,4);
AP = permute(AP,[2 3 1]);
AP0 = cell(K,1); % used a cell here because it makes the calculation of the objective function faster
for k=1:K
    AP0{k} = squeeze(AP(k,:,:));
end
AP = AP0;
AS = reshape(AS,NT*C*2,1);
AS = AS';

%___________ Cross-period and cross-category moment conditions ____________
if ~outofsample
    if new_moments
        Ncommon = zeros(N,J);           % indicator for the stores that are common to the consumer's choice set
        for i=1:N                       % across *all* periods
            Ncommon(i,1:ncommon(i)) = 1;
        end
        p0 = reshape(P,NT,J,1,K);
        p0 = cat(3,p0(:,storeindex(:,1),:,:),p0(:,storeindex(:,2),:,:));
        r1 = squeeze(sum(sum(p0.*quantity,3),2));       % NTxK
        r0 = sum(r1,2);
        q0 = squeeze(sum(sum(quantity,3),2));           % NTxK
        storeindex_C2J = zeros(C,2,J);
        for j=1:J
            for l=1:2
                storeindex_C2J(:,l,j) = storeindex(:,l)==j;
            end
        end
        % ______________ create test quantity start ___________________________
        storeindex_C2J = reshape(storeindex_C2J,1,C,2,1,J);
        d0_2 = zeros(NT,K,J);
        for j=1:J
            d0_2(:,:,j) = sum(sum(positiveamount.*storeindex_C2J(on,:,:,ok,j),2),3);
        end
        % ______________ create test quantity end _____________________________
        storeindex_C2J = reshape(storeindex_C2J,C*2,J);
        pa0 = reshape(positiveamount,NT,C*2,K);
        pa0 = permute(pa0,[1 3 2]);
        pa0 = reshape(pa0,NT*K,C*2);
        d0 = pa0*storeindex_C2J;
        d0 = reshape(d0,NT,K,J);
        % ______________ test _________________________________________________
        if ~isequal(d0_2,d0)    % the test simply compares two ways of calculating the same quantity
            disp('Warning: problem in calculation of store visits.')
        end
        % ______________ test end _____________________________________________
        os0 = sum(storepairvisited(:,A1i),2);       % NTx1
        dst0 = sum(storepairvisited.*Asp(:,:,1),2);
        % ______________ test _________________________________________________
        os0_2 = sum(storepairvisited.*onestore,2);       % the test simply compares two ways of calculating the same quantity
        dst0_2 = sum(storepairvisited.*distance,2);
        if ~isequal(os0,os0_2) | ~isequal(dst0,dst0_2)
            disp('Warning: problem in calculation of #stores & distance.')
        end
        % ______________ test end _____________________________________________
        I = ones(K,K);  % index to pick lower triangle (not including diagonal) of KxK matrix for each household
        It = triu(I);
        It = It(:);
        It = ~It;       % exclude the diagonal; we only want cross-terms
        R0 = 0;
        Q0 = 0;
        D0 = 0;
        OS0 = 0;
        DST0 = 0;
        R_in = 0;
        R_cr = 0;
        for t=2:T
            R0 = R0 + r0(period==t-1).*r0(period==t);              % Nx1
            Q0 = Q0 + q0(period==t-1,:).*q0(period==t,:);          % NxK
            D0 = D0 + d0(period==t-1,:,:).*d0(period==t,:,:);      % NxKxJ
            OS0 = OS0 + os0(period==t-1).*os0(period==t);           % Nx1
            DST0 = DST0 + dst0(period==t-1).*dst0(period==t) / 1000;
            % cross-category moments
            r_1 = r1(period==t-1,:);
            r_2 = r1(period==t,:);
            r_1t = reshape(r_1,N,1,K);
            r_2t = reshape(r_2,N,1,K);
            r_11 = r_1(:,:,ok).*r_1t(:,ok,:); % Nx(K*K)
            r_22 = r_2(:,:,ok).*r_2t(:,ok,:); % Nx(K*K)
            r_12 = r_1(:,:,ok).*r_2t(:,ok,:); % Nx(K*K)
            r_21 = r_2(:,:,ok).*r_1t(:,ok,:); % Nx(K*K)
            R_in = R_in + (r_11(:,It) + r_22(:,It) ) / 2;
            R_cr = R_cr + (r_12(:,It) + r_21(:,It) ) / 2;
        end
    end
    
    % ___________________ Initial weighting matrix ____________________________
    Wk = cell(K,1);     % initial weighting matrix for continuous moments
    for k=1:K
        comb_k = IV_comb{k};
        Lk = length(comb_k);
        W0 = zeros(L-2+Lk,L-2+Lk); % 2-Lk variables dropped from A
        for i=1:NT
            for j=1:2
                % drop components of A: cross_cat instruments
                % for categories where these are set to zero
                temp = squeeze(A(i,:,j,k,[1:7+Lk 10:end]))';
                W0 = W0 + temp*temp';
            end
        end
        Wk{k} = inv(W0/NT);
    end
    
    Wki = cell(K,1);    % initial weighting matrix for discrete moments
    for k=1:K
        comb_k = IV_comb{k};
        Lk = length(comb_k);
        W0 = zeros(L-7+Lk,L-7+Lk);  % 5+(2-Lk) variables dropped from A
        for i=1:NT
            for j=1:2
                % drop components of A: 5 time dummies, cross_cat instruments
                % for categories where these are set to zero
                temp = squeeze(A(i,:,j,k,[1 7:7+Lk 10:end]))';
                W0 = W0 + temp*temp';
            end
        end
        Wki{k} = inv(W0/NT);
    end
    
    Lb = size(Asp,3);
    olb = ones(Lb,1);
    Wb = zeros(Lb,Lb);
    for i=1:NT
        % for c=1:C
        temp = squeeze(Asp(i,:,:))';
        Wb = Wb + temp*temp';
        %end
    end
    Wb = Wb/(NT);
    Wb = inv(Wb);
    
    W = blkdiag(Wk{1},Wk{2},Wk{3},Wk{4},Wk{5},Wk{6},Wk{7},Wk{8},...
        Wki{1},Wki{2},Wki{3},Wki{4},Wki{5},Wki{6},Wki{7},Wki{8},...
        Wb);
    
    LL = size(W,1);
    if new_moments
        LL = LL + 7;
    end
end

% ________________________ Input arguments (structure) ____________________
inp.fix = 'empty'; % for use when calculating conditional derivatives
inp.instr_ix = instr_ix;
inp.period = period;
inp.nu_const = nu_const;
inp.storeindex = storeindex;
inp.quantity = quantity;
inp.positiveamount = positiveamount;
inp.storepairvisited = storepairvisited;
inp.IV_comb = IV_comb;
inp.AH = AH;
inp.AP = AP;
inp.AS = AS;
inp.ACi2 = ACi2;
inp.ACi = ACi;
inp.A1i = A1i;
inp.Asp = Asp;
inp.cat_comb = cat_comb;
inp.nupr = nupr;
inp.nugamma = reshape(nugamma,N,1,1,2);
inp.C = C;
inp.price = P;
inp.N = N;
inp.NT = NT;
inp.T = T;
inp.J = J;
inp.K = K;
inp.S = S;
inp.L = L;
inp.LP = LP;
inp.L1 = L1;
inp.L2 = L2;
inp.L3 = L3;
inp.Li = Linteract; 
inp.on = on;
inp.oj = oj;
inp.ok = ok;
inp.o2 = o2;
inp.os = os;
inp.olA = olA;
inp.ix2 = ix_JJ2C;
inp.x = x;
inp.xp = xp;
inp.z = z;
inp.nu1 = nu1;
inp.nu2 = nu2;
inp.nu3 = nu3;
inp.chain0 = squeeze(chain9(:,:,1,:));
inp.discount = 0; % for use in counterfactuals
inp.new_moments = new_moments;
inp.cat_1st_instr = cat_1st_instr;
if ~outofsample
    inp.nu3_78 = nu3_78;
    inp.W = W;
    inp.Wb = Wb;
    inp.Lb = Lb;
    inp.LL = LL;
    inp.olb = olb;
    if new_moments
        inp.Ncommon = Ncommon;
        inp.It = It;
        inp.storeindex_C2J = storeindex_C2J;
        inp.p0 = p0;
        inp.R0_star = R0;
        inp.Q0_star = Q0;
        inp.D0_star = D0;
        inp.OS0_star = OS0;
        inp.DST0_star = DST0;
        inp.R_in_star = R_in;
        inp.R_cr_star = R_cr;
    end
    inp2.chain = chain9;
    inp2.chain2 = chain16;
    inp2.frmcl = frmcl; % labels
    inp2.frmcl9 = frmcl9;
    inp2.catcl = catcl;
    inp2.stchcl = stchcl;
    inp2.hhchcl = hhchcl;
    inp2.rncl = rncl;
    inp2.spcl = spcl;
    inp2.A = A;
else
    inp2 = [];
end
