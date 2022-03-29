function [quant_pred,psamt_pred,prob_vis_storepair,Q,jbest] = TSS_quantities(theta,inp)

ix2 = inp.ix2;
price = inp.price;
C = inp.C;
NT = inp.NT;
T = inp.T;
J = inp.J;
K = inp.K;
S = inp.S;
on = inp.on;
oj = inp.oj;
ok = inp.ok;
x = inp.x;
xp = inp.xp;
z = inp.z;
L1 = inp.L1;
L2 = inp.L2;
L3 = inp.L3;
Li = inp.Li;
chain0 = inp.chain0;
cat_comb = inp.cat_comb;
CC = size(cat_comb,1);
storeindex = inp.storeindex;
nupr = inp.nupr;
nupr = repmat(nupr,T,1);
nu_const = repmat(inp.nu_const,T,1);
nugamma = repmat(inp.nugamma,T,1);
nu1 = repmat(inp.nu1,T,1);
nu2 = inp.nu2;
nu2 = nu2(:,oj,:);
nu2 = repmat(nu2,T,1);
nu3 = inp.nu3;

D = length(theta);
theta = reshape(theta,1,D);
scaleparam = theta(1:K-1);
scaleparam = reshape(scaleparam,1,K-1);
scaleparam = [scaleparam 1];
scaleparam0 = reshape(scaleparam,1,1,K);
constant0 = theta(K);
betax = theta(K+1:K+L1);  % store characteristic
L0 = K+L1; % simplify expression for index
betaz = theta(L0+1:L0+L2-1); % consumer chars. (income not included, hence
% length is L2-1)
betaz = reshape(betaz,L2-1,1);
sigma = theta(L0+L2:L0+L2+2); % random terms
sigma = reshape(sigma,3,1);
const_sigma = theta(L0+L2+3);
L0 = L0+L2+3; % simplify expression for index
lambda0 = reshape(theta(L0+1:L0+K),K,1); % second-order terms
s_int = theta(L0+K+1:L0+K+Li); % interaction terms
alpha0 = theta(L0+K+Li+1);
alpha1 = theta(L0+K+Li+2);
L0 = L0+K+Li+2; % simplify expression for index
gamma0 = theta(L0+1:L0+2); % shopping cost parameters
gm_sg = theta(L0+3:L0+4);  
delta0 = theta(L0+5:end);

%_______________components of discrete (store-pair) utility________________
xp0 = reshape(xp,NT,J,J,L3);
SCt = - ( gamma0(1) + gm_sg(1)*nugamma(:,oj,oj,1) );
SCd = - ( gamma0(2) + gm_sg(2)*nugamma(:,oj,oj,2) );
xpgamma = xp0(:,:,:,1) .* SCt  +  xp0(:,:,:,2) .* SCd;
xpgamma = reshape(xpgamma,NT,J*J);
xpgamma = xpgamma(:,ix2);

%_______________components of continuous utility first-order terms_________
% Store characteristics
xbeta = 0;
for l1=1:L1
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
apr1 = apr1(:,oj,ok).*price;                    % price term
mu = apr1 + scaleparam0(on,oj,:).*(constant0 + xzbeta(:,:,ok) + nusigma);
% NTxJxKxS
% for each j chain is 1 only for one of the S pages and 0 for the others,
% so summing does not lead to any duplication
firmcat = zeros(NT,J,K);
for k=1:K
    fc = [0 delta0((k-1)*(S-1)+1:k*(S-1))];
    fc = reshape(fc,1,1,S);
    firmcat(:,:,k) = sum(fc(on,oj,:).*chain0,3);
end
mu = firmcat + mu;

%__________put second-order parameters in matrix form______________________
Lambda = zeros(K,K);
l = 1;
for t=1:CC                       % CC = size(cat_comb,1). l.26
    c_cat_t = cat_comb{t};      % e.g. cat_comb{1} = [3 5] for cats. 3 & 5
    no_cat = size(c_cat_t,2);   % number of categories that are interacted
    if no_cat>1
        combs = combnk(c_cat_t,2);
        C0 = size(combs,1);     % #pairs = #interaction params.
        linInd = sub2ind([K K],combs(:,1),combs(:,2));
        for cc=1:C0
            Lambda(linInd(cc)) = s_int(l);
            l = l + 1;
        end
    end
end
% put the (k,k') interaction parameter in entry (k',k) too.
Lambda = Lambda + Lambda' + diag(lambda0);

%_____________Continuous quantities and utilities__________________________
% first-order terms by store within each store pair
mu1 = mu(:,storeindex(:,1),:);
mu2 = mu(:,storeindex(:,2),:);

%_____________Counterfactual discount calculations________________________
discount = inp.discount; % discount = 0 during estimation
if sum(discount)~=0 % for counterfactual calculations
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
%______________________________________________________________________
% set first-order term to zero unless greater than for the other store
if strcmp(inp.fix, 'empty') || strcmp(inp.fix, 'save') || strcmp(inp.fix, 'c')
    jbest = mu1 >= mu2;             % j best, as opposed to j' best, in c=(j,j')
    if strcmp(inp.fix, 'save')
        save jbest jbest
    end
elseif strcmp(inp.fix, 'dc') % if calculating elasticities for fixed 'd', load saved jbest
    load jbest
end

muC = mu2;
muC(jbest) = mu1(jbest);
muC = permute(muC,[3 1 2]);
NTC = NT*C;
mu = reshape(muC,K,NTC);

Q = zeros(K,NTC);
U = zeros(CC,NTC);
for t=1:CC
    ix = cat_comb{t};
    [Q(ix,:),U(t,:)] = QP_nonneg(Lambda(ix,ix),mu(ix,:));
end
U = sum(U,1);
Q = reshape(Q,K,NT,C);
Q = permute(Q,[2 3 1]);
U = reshape(U,1,NT,C);
U = permute(U,[2 3 1]);

%__________________unconditional predictions_______________________________
v = U + xpgamma;
ev = exp(v);
evdenom = sum(ev,2);
evdenom = evdenom(:,ones(C,1));

if strcmp(inp.fix, 'empty') || strcmp(inp.fix, 'save')
    pp = ev./evdenom;
    prob_vis_storepair = pp;
    ppK = pp(:,:,ones(K,1));
    % store j (within c) chosen for category k
    if strcmp(inp.fix, 'save')
        save ppK ppK
    end
elseif strcmp(inp.fix, 'c') || strcmp(inp.fix, 'dc') % if calculating elasticities for fixed 'c', load saved ppK
    load ppK
end

quant_pred = cat(3,reshape(  Q .* jbest .* ppK , NT,C,1,K),...
    reshape(  Q .* ~jbest .*ppK,NT,C,1,K));
psamt_pred = cat(3,reshape( (Q>0) .* jbest .* ppK , NT,C,1,K),...
    reshape( (Q>0) .* ~jbest .*ppK , NT,C,1,K));

