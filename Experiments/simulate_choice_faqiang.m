function [quant_pred,psamt_pred,prob_vis_storepair,Q,jbest] = simulate_choice_faqiang(theta,inp)
% This script is adapted from TSS_quantities. It is used for predicting
% purchase behaviors given household characteristics, store characteristics  
% and so on under certain parameter values.
% 
% This script predict household's choices given characteristics on
% 1.which store to visit/spend positive expenditure: psamt_pred,prob_vis_storepair
% 2.which category purchased in each store: quant_pred,jbest
% 3.what quantity for each category purchased: Q
% 
% Identifiers for observations:
% i,t   : individual and time, or a "purchase"
% c,j   : a store choice c (2*1), and a particular store j in c
% k     : a category k
% 
% Above choices data can be read from the output
% quant_pred         : predicted quantity for each i,t,c,j,k, size=NT,C,2,K
% psamt_pred         : predicted whether positive amount, value=0or1, size=NT,C,2,K
% prob_vis_storepair : predicted probability of visit a score pair c by i,t, size=NT,C,1
% Q                  : quantity for each category for each purchase-storechoice i-t-c, size=K,NTC
% jbest              : whether the first of a store choice c is chosen for a particular purchase of a category i-t-k, size=NT,C,K
%
% Arguments:
% inp  : a struct including all the non-choice data. 
% theta: a D*1 vector of all the structrual parameter to be estimated


% Read data from inp.. NoComment means scalar.. each * separates a
% dimension
ix2 = inp.ix2; % JJ * 1, vectorize upper tri of store pair matrix
price = inp.price; %  NT*J*K, price
C = inp.C; % store combination, C=30*29/2+30=465
NT = inp.NT; %
T = inp.T; % periods
J = inp.J; % size of choice set
K = inp.K; % size of category set
S = inp.S; % # firms
on = inp.on; % NT*1, on, oj, ok are all matrix repetitors
oj = inp.oj; % J*1
ok = inp.ok; % K*1
x = inp.x; % NT*J, storechars(:,:,1), i.e. store size,
xp = inp.xp; % NT*J*J*L3, distance matrix, NT*J*J*L3 = each i in t * store A * store B * (onestore & distance)
z = inp.z; % N*4, 4 chars of matrix hhchars
L1 = inp.L1; % number of store chars enterring into taste mu
L2 = inp.L2; % number of hh chars enterring into taste mu
L3 = inp.L3; %  L3=2, two distance related utility
Li = inp.Li; % number of interactive elasticity
chain0 = inp.chain0; % NT*J*S, logical array to pick stores that are in a particular chain.
cat_comb = inp.cat_comb; % cat_comb = cell(4,1), this is a struct of category combination.
CC = size(cat_comb,1); % CC=4, there are four category combination 

storeindex = inp.storeindex; % C*2

% All the individual unobserved characteristic
nupr = inp.nupr; % N*1*1
nupr = repmat(nupr,T,1); % NT*1*1, notice this is becasue nupr=\nu_i^{alpha} is individual specific
nu_const = repmat(inp.nu_const,T,1); % NT*1
nugamma = repmat(inp.nugamma,T,1); % NT*1*1*2, nugamma = randn(N,1,1,2)
nu1 = repmat(inp.nu1,T,1); % see below
nu2 = inp.nu2;
nu2 = nu2(:,oj,:);
nu2 = repmat(nu2,T,1);
nu3 = inp.nu3;


% Read all the parameters and rename them
D = length(theta);  % number of parameters
theta = reshape(theta,1,D);

% beta_0k
scaleparam = theta(1:K-1); 
scaleparam = reshape(scaleparam,1,K-1);
scaleparam = [scaleparam 1];
scaleparam0 = reshape(scaleparam,1,1,K);
constant0 = theta(K);

% beta_1 and beta_2 (or beta_x, beta_z)
betax = theta(K+1:K+L1);  % store characteristic
L0 = K+L1; % simplify expression for index
betaz = theta(L0+1:L0+L2-1); % consumer chars. 
betaz = reshape(betaz,L2-1,1);

% sigma (the spread parameters of \nu) 
sigma = theta(L0+L2:L0+L2+2); 
sigma = reshape(sigma,3,1);
const_sigma = theta(L0+L2+3);

% lambda (quadratic parameters)
L0 = L0+L2+3; 
lambda0 = reshape(theta(L0+1:L0+K),K,1); % second-order terms
s_int = theta(L0+K+1:L0+K+Li); % interaction terms

% alpha (price)
alpha0 = theta(L0+K+Li+1);
alpha1 = theta(L0+K+Li+2);

% gamma (shopping cost)
L0 = L0+K+Li+2; % 
gamma0 = theta(L0+1:L0+2); % shopping cost parameters, i.e. gamma11 and gamma21
gm_sg = theta(L0+3:L0+4);  % gamma12 and gamma22

% Firm-category fixed effect, in total, K*(S-1) of them
delta0 = theta(L0+5:end);  


% default value. Not relevant for simulation.
inp.discount=0;
inp.fix = 'empty';


%_______________components of discrete (store-pair) utility________________
xp0 = reshape(xp,NT,J,J,L3);
SCt = - ( gamma0(1) + gm_sg(1)*nugamma(:,oj,oj,1) ); % shopping cost of store choice (one stop or two stop)
SCd = - ( gamma0(2) + gm_sg(2)*nugamma(:,oj,oj,2) ); % shopping cost of distance
xpgamma = xp0(:,:,:,1) .* SCt  +  xp0(:,:,:,2) .* SCd; 
xpgamma = reshape(xpgamma,NT,J*J); % each entry denotes for this purchase (row) and store comb(col), the SC
xpgamma = xpgamma(:,ix2); % avoid repetition, for xpgamma_storeAB = xpgamma_store_BA

%_______________components of continuous utility first-order terms_________
% Store characteristics
xbeta = 0;
for l1=1:L1 % adding store char one by one
    xbeta = xbeta + betax(l1)*x(:,:,l1);
end
% Demographic characteristics
zbeta = 0;
counter = 1;
for l2=[1 3:L2]   % hhsize, time dummies
    zbeta = zbeta + betaz(counter)*z(:,l2);
    counter = counter + 1;
end
xzbeta = xbeta + zbeta(:,oj);

% Random draws
nusigma = const_sigma*nu_const + sigma(3)*nu3;  % random terms that do not vary by store and category
nusigma = nusigma(:,oj,ok) + sigma(1)*nu1 + sigma(2)*nu2;
a_inc = alpha1./z(:,2);
apr1 = -(alpha0 + a_inc).*nupr;
apr1 = apr1(:,oj,ok).*price;                  % price term, price is NT*J*K

% Taste, aggregate above up
% mu is NT*J*K
mu = apr1 + scaleparam0(on,oj,:).*(constant0 + xzbeta(:,:,ok) + nusigma);
firmcat = zeros(NT,J,K); % \mu_{fk} term. category-firm fixed effect
for k=1:K % looping over all category
    fc = [0 delta0((k-1)*(S-1)+1:k*(S-1))]; % fc = firm 0,1,2,3,4,5,6,7,8 for category k, and k moves next
    fc = reshape(fc,1,1,S); %
    firmcat(:,:,k) = sum(fc(on,oj,:).*chain0,3); % chain0 is firm indicator, chain0 is NT*J*S, firmcat(:,:,k) is NT*J
end

mu = firmcat + mu;  % mu is NT*J*K, equation (9)

%__________put second-order parameters in matrix form______________________
Lambda = zeros(K,K);
l = 1;
for t=1:CC                       % CC = size(cat_comb,1). l.26
    c_cat_t = cat_comb{t};      % e.g. cat_comb{1} = [3 5] for cats. 3 & 5
    no_cat = size(c_cat_t,2);   % number of categories that are interacted
    if no_cat>1
        combs = combnk(c_cat_t,2);
        C0 = size(combs,1);     % #pairs = #interaction params.
        linInd = sub2ind([K K],combs(:,1),combs(:,2)); % mind the sub2ind argument
        for cc=1:C0
            Lambda(linInd(cc)) = s_int(l); % s_int is entry of interactive taste
            l = l + 1;
        end
    end
end
% put the (k,k') interaction parameter in entry (k',k) too.
Lambda = Lambda + Lambda' + diag(lambda0);  % mind that these are all parameters

%_____________Continuous quantities and utilities__________________________
% first-order terms by store within each store pair
mu1 = mu(:,storeindex(:,1),:); % utility got from store A in each trip
mu2 = mu(:,storeindex(:,2),:); % utility got from store B , both mu1,mu2 are NT*C*K

% explanation: as know from previous discussion, mu1 is NT*C*K matrix.
% In this cube (of matrix), along the second dimension, mu1 just gives the
% utility of the store A in that purchase for that cat while mu2 just gives the utility
% of the store B in that purchase for that cat.

%_____________Counterfactual discount calculations________________________
% not relevant for simulation
discount = inp.discount; % discount = 0 during estimation
% ignore this for now
if sum(discount)~=0 % for counterfactual calculations, in TSS_input, discount=0
    frm = inp.frm;
    ACi2 = inp.ACi2;
    a1 = apr1(:,storeindex(:,1),:);
    a2 = apr1(:,storeindex(:,2),:);
    ld = length(discount);
    if ld==1
        % store pairs including only firms frm's stores
        temp = ACi2(:,:,1)==frm & ACi2(:,:,2)==frm;
        temp = temp(:,:,ok);
        temp1 = mu1 - discount*temp.*a1;
        temp2 = mu2 - discount*temp.*a2;
    else
        temp1 = zeros(NT,C,K);
        temp2 = temp1;
        for k=1:ld
            % store pairs including only firms frm's stores
            temp = ACi2(:,:,1)==frm & ACi2(:,:,2)==frm;
            %temp = temp(:,:,:,ok);
            temp1(:,:,k) = mu1(:,:,k) - discount(k)*temp.*a1(:,:,k);
            temp2(:,:,k) = mu2(:,:,k) - discount(k)*temp.*a2(:,:,k);
        end
    end
    mu1 = temp1;
    mu2 = temp2;
end


% Create jbest matrix
if strcmp(inp.fix, 'empty') || strcmp(inp.fix, 'save') || strcmp(inp.fix, 'c') % ignore this line of the code
    jbest = mu1 >= mu2;             % j best, as opposed to j' best, in c=(j,j')
    if strcmp(inp.fix, 'save')  % strcmp, compare string.
        save jbest jbest
    end
elseif strcmp(inp.fix, 'dc') % if calculating elasticities for fixed 'd', load saved jbest
    load jbest
end


%___From all above, mu is ready. Now consider quantity and shopping cost__


%_____________non-linear programming to determine quanitty______________

muC = mu2; % initial MuC as all utility from store B, for a nt-k
muC(jbest) = mu1(jbest); % whenever store A is better for a nt-k. If yes, replace muC with mu from store A
muC = permute(muC,[3 1 2]); % muC now K*NT*C, 
NTC = NT*C;
mu = reshape(muC,K,NTC); % reshape to better feed the QP_nonneg

Q = zeros(K,NTC); % quantity for each cat of purchase under store choice c
U = zeros(CC,NTC); % optimal utility, this is the w in the paper, CC is how many cat-combinations, in this case 4
for t=1:CC % looping for all possible complementarity (CC)
    
    ix = cat_comb{t}; % temporary category combination index
    [Q(ix,:),U(t,:)] = QP_nonneg(Lambda(ix,ix),mu(ix,:)); % quadratic programming, see script outside.
    
end
U = sum(U,1); % adding across category combination so we get total consumption utility
Q = reshape(Q,K,NT,C);
Q = permute(Q,[2 3 1]); % NT*C*K
U = reshape(U,1,NT,C);
U = permute(U,[2 3 1]); % NT*C*1

%__________________unconditional predictions_______________________________
v = U + xpgamma; % recall the quantity decision is made given c, so we have to add back shopping cost
ev = exp(v); % NT*C*1
evdenom = sum(ev,2); % logit denominator, it is added within nt, across all its c
evdenom = evdenom(:,ones(C,1)); 

% ignore the first line condition, it is passed
if strcmp(inp.fix, 'empty') || strcmp(inp.fix, 'save')
    pp = ev./evdenom;  % this is equation (18), NT*C*1
    prob_vis_storepair = pp;
    ppK = pp(:,:,ones(K,1)); % just replicate it for K times
    % store j (within c) chosen for category k
    if strcmp(inp.fix, 'save')
        save ppK ppK
    end
elseif strcmp(inp.fix, 'c') || strcmp(inp.fix, 'dc') % if calculating elasticities for fixed 'c', load saved ppK
    load ppK
end

% The  Q .* jbest .* ppK is what we want: at the chosen store, what is the
% quantity timing the probability of making this choice. As for  Q .* ~jbest .* ppK
% this is the predicted purchase quantity for the non-chosen store
quant_pred = cat(3,reshape(  Q .* jbest .* ppK , NT,C,1,K),...
    reshape(  Q .* ~jbest .*ppK,NT,C,1,K));
psamt_pred = cat(3,reshape( (Q>0) .* jbest .* ppK , NT,C,1,K),...
    reshape( (Q>0) .* ~jbest .*ppK , NT,C,1,K));
%positive amount prediction

