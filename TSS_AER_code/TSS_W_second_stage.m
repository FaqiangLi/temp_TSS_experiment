function W = TSS_W_second_stage(theta,inp,inp2)

K = inp.K;
N = inp.N;
W = inp.W;
A = inp2.A; 
Asp = inp.Asp;
clear inp2 % free memory
instr_ix = inp.instr_ix;
LL = inp.LL;
new_moments = inp.new_moments;
if new_moments
    Ncommon = inp.Ncommon;
end

quantity = inp.quantity;
positiveamount = inp.positiveamount;
storepairvisited = inp.storepairvisited;

[quant_pred,psamt_pred,prob_vis_storepair] = TSS_quantities(theta,inp);
res1 = quantity - quant_pred;
res2 = positiveamount - psamt_pred;
res3 = storepairvisited - prob_vis_storepair;
if new_moments
    G4 = TSS_nonlin_moments(quant_pred,psamt_pred,prob_vis_storepair,inp);
end

G1 = A.*res1(:,:,:,:,ones(size(A,5),1));
A = A(:,:,:,:,[1 7:end]); % drop time dummies
G2 = A.*res2(:,:,:,:,ones(size(A,5),1));
clear A
G3 = Asp.*res3(:,:,ones(size(Asp,3),1));
clear Asp
%__________________________________________________________________________
% check that GMM objective has the same value as in TSS_gmm(x,inp)
if ~new_moments
    g1 = squeeze(sum(sum(sum(G1,1),2),3))';
    g1 = g1(:);
    g2 = squeeze(sum(sum(sum(G2,1),2),3))';
    g2 = g2(:);
    g3 = squeeze(sum(sum(G3,1),2));
    g = [g1; g2; g3]/N;
    g = g(instr_ix);
    gmm1 = g'*W*g;
    gmm2 = TSS_gmm(theta,inp);
    if abs(gmm1-gmm2)>1e-12
        disp('Warning: check GMM calculation')
    end
end
%__________________________________________________________________________
W = zeros(LL,LL);
for i=1:N
    temp1 = [];
    temp2 = [];
    for k=1:K
        temp1 = [temp1; squeeze(sum(sum(G1(i,:,:,k,:) + G1(i+N,:,:,k,:),2),3))];
        temp2 = [temp2; squeeze(sum(sum(G2(i,:,:,k,:) + G2(i+N,:,:,k,:),2),3))];
    end
    temp3 = squeeze(sum(G3(i,:,:) + G3(i+N,:,:),2));
    temp = [temp1; temp2; temp3];
    temp = temp(instr_ix);
    if new_moments
        temp4 = G4(:,i);
    else
       temp4 = []; 
    end
    temp = [temp; temp4];
    W = W + temp*temp';
end
W = W/N;