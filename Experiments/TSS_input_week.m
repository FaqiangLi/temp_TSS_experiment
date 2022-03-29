function inp_TS = TSS_input_week(t_vector)

% 1 <= TS <= 156                number of time periods

%___________________Settings_______________________________________________
onetrip_dist = false;           % use distance of visiting both stores on one trip                     
N = 2000;                       % number of households
TS = 156;

%___________________Import data____________________________________________
% import data from data directory
dir = 'priceanalysis_data/';

cat_prices0 = table2array(readtable([dir 'cat_prices.csv']));
hhchars0 = table2array(readtable([dir 'HHchars2.csv']));
ntjchars0 = table2array(readtable([dir 'NTJchars2.csv']));
ntjfirmnum0 = table2array(readtable([dir 'NTJfirmnum2.csv']));

hhchars0 = repmat(hhchars0,TS,1);      % hhchars are constant across time

%____________test__________________________________________________________
% check that there are not observations that are NaN
test = zeros(4,1);
test(1) = sum(isnan(cat_prices0(:)));
test(2) = sum(isnan(ntjfirmnum0(:)));
test(3) = sum(isnan(ntjchars0(:)));
test(4) = sum(isnan(hhchars0(:)));
if sum(test)~=0
    disp('Warning: some data array contains NaN entries.')
    %return
end
%____________end test______________________________________________________

hh_codeTN0 = cat_prices0(:,1);   % first column in cat_prices is househould number: household codes
hh_codeN0 = unique(hh_codeTN0,'stable');    % list of unique household codes
J = 30;                         % stores in each person's choice set
K = 8;                          % categories of groceries
C = 465;                      % no. of households * no. of time periods


%____________test__________________________________________________________
% check that households are the same, and in same order, as in estimation
load hh_code_estimation_sample
if sum(abs(hh_codeN0(:) - hh_code(:))) ~= 0
    disp('Check list of households.')
end

%____________end test______________________________________________________
period = cat_prices0(:,3);

CC = length(t_vector);
cc = 0;  % counter
for t=t_vector
    cc = cc + 1;
    if cc==1 | cc==CC |  mod(cc, 10) == 0 % print 1st and thereafter every 10 iterations
        disp(['Data generated for time period number ' num2str(cc) ' of ' num2str(CC)])
    end
    ix = period==t;
    ntjch_ix = ntjchars0(:,3)==t;
    % floor area, distance [30 nearest stores for each consumer]
    ntjchars = ntjchars0(ntjch_ix,:);
    % chain number (1-16) first store ... 30th store, for each consumer
    ntjfirmnum = ntjfirmnum0(ix,:);
    hhchars = hhchars0(ix,:);
    % prices in each of 8 categories for each of 30 stores, for each consumer
    cat_prices = cat_prices0(ix,:);

    week0 = cat_prices(:,2);
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
    
    oj = ones(J,1);
    ok = ones(K,1);
      
    temp = true(J,J);                   % store / store pair indices
    temp = triu(temp);
    ix_JJ2C = false(J*J,1);
    ix_JJ2C(temp) = true;
    J1 = (1:J)';
    storeindex = zeros(J*J,2);          % all combinations (j,j'), 1<=j,j'<=J
    storeindex(:,1) = kron(ones(J,1),J1);
    storeindex(:,2) = kron(J1,ones(J,1));
    storeindex = storeindex(ix_JJ2C,:); % combinations (j,j') s.th 1<=j<=j'<=J
    
    % demographic variables: adults+children, income, time dummies
    % income divided by 10, so it is in tens of pounds, like expenditure
    z = [sum(hhchars(:,[10 11]),2) hhchars(:,end)/10 quarter];
    % number of demographic characteristics (including time dummies)
    z(:,2) = z(:,2)./z(:,1); % change income to per capita income
    
    x0 = ntjchars;                      % [Hhnumber week t Ncommon n salesarea easting northing dist]
    x0(:,5) = log(x0(:,5));             % log of store size (salesarea)
    x0(:,6) = x0(:,6)/1000;             % convert to metres. (dist = x0(:,8) is already in metres).
    x0(:,7) = x0(:,7)/1000;
    % log of storesize
    storechars = zeros(N,J,4);
    for i=1:N
        ix_i = x0(:,1)==hh_code(i);
        storechars(i,:,:) = x0(ix_i,5:end);      % logsalesarea easting northing dist
    end
       
    % store characteristics which enter first-order term of continuous utility
    x = storechars(:,:,1);
    
    % NxJxKxJxL3. store-PAIR characteristics, for discrete utility
    % number of such characteristics [two-store or not, combined distance]
    % [log_salesarea easting northing dist]
    east = storechars(:,:,2);
    north = storechars(:,:,3);
    dist = storechars(:,:,4);
    L3 = 2;
    dist = reshape(dist,N,J,1);
    dist = dist(:,:,oj);
    dist_t = permute(dist,[1 3 2]);         % transpose
    north = reshape(north,N,J,1);
    north = north(:,:,oj);
    north_t = permute(north,[1 3 2]);       % transpose
    east = reshape(east,N,J,1);
    east = east(:,:,oj);
    east_t = permute(east,[1 3 2]);         % transpose
    
    xp_triangular = zeros(N,J,J,L3);       % distance home -> A -> B -> home
    for i=1:N
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
    
    xp_twotrips = zeros(N,J,J,L3);         % distance home -> A -> home -> B -> home
    for i=1:N
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
    
    % NxJ. Entry i,j gives the number (from 1 to S) of the chain to which
    % consumer i's j-th store belongs.
    firmnum = ntjfirmnum(:,4:end);    
    %  Aggregating to nine firms
    firmnum9 = zeros(N,J);     % array of firm codes
    % ASDA                      1
    ix_temp = firmnum==2;
    firmnum9(ix_temp) = 1;
    % DISC (ALDI, LIDL, NETTO)  2
    ix_temp = (firmnum==1|firmnum==7|firmnum==10);
    firmnum9(ix_temp) = 2;
    % ICELAND                   3
    ix_temp = firmnum==5;
    firmnum9(ix_temp) = 3;
    % MORRISONS                 4
    ix_temp = firmnum==8;
    firmnum9(ix_temp) = 4;
    % MS                        5
    ix_temp = firmnum==9;
    firmnum9(ix_temp) = 5;
    % OTHER (OTHER, BUDGEN, KWIKSAVE, SAFEWAY, COOP, SOMERFIELD) 6
    ix_temp = (firmnum==3|firmnum==4|firmnum==6|firmnum==11|firmnum==12|firmnum==14);
    firmnum9(ix_temp) = 6;
    % SAINSBURY                 7
    ix_temp = firmnum==13;
    firmnum9(ix_temp) = 7;
    % TESCO                     8
    ix_temp = firmnum==15;
    firmnum9(ix_temp) = 8;
    % WAITROSE                  9
    ix_temp = firmnum==16;
    firmnum9(ix_temp) = 9;
    
    % NTxJxS. Logical array to pick stores that are in a particular chain.
    chain9 = zeros(N,J,9);
    for s=1:9
        chain9(:,:,s) = firmnum9==s;
    end
    % NxTxJxKxS. Expand for categories
    chain9 = reshape(chain9,N,J,1,9);
    chain9 = chain9(:,:,ok,:);
    
    % same with chains not aggregated
    % NTxJxS. Logical array to pick stores that are in a particular chain.
    chain16 = zeros(N,J,16);
    for s=1:16
        chain16(:,:,s) = firmnum==s;
    end
    % NxTxJxKxS. Expand for categories
    chain16 = reshape(chain16,N,J,1,16);
    chain16 = chain16(:,:,ok,:);
    
    chainnumber16 = zeros(N,C,2,K);
    for i=1:N
        chainnumber16(i,:,1,1) = firmnum(i,storeindex(:,1));
        chainnumber16(i,:,2,1) = firmnum(i,storeindex(:,2));
        chainnumber16(i,:,:,:) = chainnumber16(i,:,:,ok);
    end
    ACi2 = chainnumber16(:,:,:,1);
    
    % prices
    p0 = cat_prices(:,4:end);
    P = zeros(N,J,K);
    for k=1:K
        temp2 = p0(:,k:K:J*K);
        P(:,:,k) = temp2;
    end
    
    % _______________Input arguments (nonscalar structure) ____________________
    inp_TS(cc).storeindex = storeindex;
    inp_TS(cc).price = P;
    inp_TS(cc).x = x;
    inp_TS(cc).xp = xp;
    inp_TS(cc).z = z;
    inp_TS(cc).chain = chain9; 
    inp_TS(cc).chain2 = chain16; 
    inp_TS(cc).ACi2 = ACi2;
    inp_TS(cc).A1i = storeindex(:,1)==storeindex(:,2);
end
