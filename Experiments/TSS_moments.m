function [G1,G2,G3,G4] = TSS_moments(theta,inp)

new_moments = inp.new_moments;
cat_1st_instr = inp.cat_1st_instr;
C = inp.C;
K = inp.K;
NT = inp.NT;
LP = inp.LP;

% Nx6 [household size (adults + children), 5 time dummies]
AH = inp.AH;
% NxCx2xKx5 [pricinstr,priceinstr/(hhincome/hhsize),...
%             priceinstrument_cross_cat (ut to 2), priceinstrument_cross_store]
AP = inp.AP;
% NxCx2x2 [storesize,storepairdistance]
AS = inp.AS;
% NxCx2 [code (1-9) for each store]
ACi= inp.ACi;
if cat_1st_instr
    % Cx1 [dummy for whether one-stop 'pair' or not]
    A1i = inp.A1i;
end
% NxCx3 [storepairdistance,onestoredummy,...
% mean priceinstrument across stores in store pair and categories]
Asp = inp.Asp;

quantity = inp.quantity;
positiveamount = inp.positiveamount;
storepairvisited = inp.storepairvisited;
% [NxCx2x8  NxCx2x8    NxC]
[quant_pred,psamt_pred,prob_vis_storepair] = TSS_quantities(theta,inp);
% residuals
res1 = quantity - quant_pred;
res2 = positiveamount - psamt_pred;
res3 = storepairvisited - prob_vis_storepair;

temp = squeeze(sum(sum(res1,2),3));
G1a = zeros(6,8);
temp2 = squeeze(sum(sum(res2,2),3));
G2a = zeros(1,8);
for k=1:8 % household characteristics
    temp_b  = sum(AH.*temp(:,k*ones(6,1)))';
    G1a(:,k) = temp_b; % category quantity residuals
    temp_b2  = sum(AH(:,1).*temp2(:,k))';
    G2a(:,k) = temp_b2; % category 'visit' residuals
end

r1 = reshape(res1,NT*C*2,K);
r2 = reshape(res2,NT*C*2,K);
G1bc = zeros(LP,K);
G2bc = G1bc;
temp1 = [AS*r1 AS*r2];
for k=1:8
    APk = AP{k};
    G1bc(:,k) = [APk*r1(:,k); temp1(:,k)];
    G2bc(:,k) = [APk*r2(:,k); temp1(:,K+k)];
end

G1d = zeros(8,K); % chain dummies (1-8)
G2d = zeros(8,K);
for s=1:8
    itemp = (ACi(:)==s)';
    G1d(s,:) = itemp*r1;
    G2d(s,:) = itemp*r2;
end

if cat_1st_instr
    G1e = squeeze(sum(sum(sum(res1(:,A1i,:,:),1),2),3))'; % one-stop dummy
    G2e = squeeze(sum(sum(sum(res2(:,A1i,:,:),1),2),3))'; % one-stop dummy
else
    G1e = [];
    G2e = [];
end

G1f = sum(r1); % constant
G2f = sum(r2); % constant
G3a = reshape(Asp,NT*C,6)'*res3(:);

G1 = [G1a; G1bc; G1d; G1f; G1e];
G2 = [G2a; G2bc; G2d; G2f; G2e];
G1 = G1(:);
G2 = G2(:);
G3 = G3a(:);
 
if new_moments
    G4 = TSS_nonlin_moments(quant_pred,psamt_pred,prob_vis_storepair,inp);
    G4 = sum(G4,2);
else
    G4 = [];
end



