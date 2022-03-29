function [visits,profit,revenue] = TSS_profit_visits(theta,inp,mcost,ix_g,frm)

% TSS_profit_visits returns 3 Kx2 matrices of visits, profits, revenue for
% each category and each of two consumer groups for a given firm 'frm'.

C = inp.C;
K = inp.K;
NT = inp.NT;
S = 16;
J = inp.J;
ACi = inp.ACi2;
o2 = ones(2,1);
chain = inp.chain2;
storeindex = inp.storeindex;

[~,~,P,Q,jbest] = TSS_quantities(theta,inp);
% price
pr = inp.price;
price = zeros(NT,C,2,K);
price(:,:,1,:) = pr(:,storeindex(:,1),:);
price(:,:,2,:) = pr(:,storeindex(:,2),:);
% marginal cost
mc = reshape(mcost,1,1,K,S);
mc = mc(ones(NT,1),ones(J,1),:,:);
mc = sum(mc.*chain,4);
mc0 = mc;
mc = zeros(NT,C,2,8);
mc(:,:,1,:) = mc0(:,storeindex(:,1),:);
mc(:,:,2,:) = mc0(:,storeindex(:,2),:);
%markup
markup = price-mc;
% revenue (per unit)
rev = price;
% given the store pair, use this store or not for this category
% ("store visit")
st1 = reshape(jbest,NT,C,1,K);
st2 = reshape(~jbest,NT,C,1,K);
stvis = cat(3,st1,st2);
% quantity purchased in the store pair
Q2 = reshape(Q,NT,C,1,K);
Q2 = Q2(:,:,o2,:);
% profits from a store pair if visited
prof = markup.*Q2.*stvis;
rev = rev.*Q2.*stvis;
quant = Q2.*stvis;

Afrm = ACi==frm;
G = size(ix_g,2);
visits = zeros(K,G);
profit = zeros(K,G);
revenue = zeros(K,G);
for k=1:K
    for g=1:G;
        ix = ix_g(:,g);
        firm = Afrm(ix,:,:,:);
        pr = P(ix,:,:,:);
        visits(k,g) = sum(sum(sum(firm.*(quant(ix,:,:,k)>0),3) .*pr));
        profit(k,g) = sum(sum(sum(firm.*prof(ix,:,:,k),3) .*pr));
        revenue(k,g) = sum(sum(sum(firm.*rev(ix,:,:,k),3) .*pr));
    end
end



