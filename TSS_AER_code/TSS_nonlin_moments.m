function G4 = TSS_nonlin_moments(quant_pred,psamt_pred,prob_vis_storepair,inp)

N = inp.N;
T = inp.T;
C = inp.C;
ok = inp.ok;
J = inp.J;
period = inp.period;
It = inp.It;
storeindex_C2J = inp.storeindex_C2J;
p0 = inp.p0;
Ncommon = inp.Ncommon;
R0_star = inp.R0_star;
Q0_star = inp.Q0_star;
D0_star = inp.D0_star;
OS0_star = inp.OS0_star;
DST0_star = inp.DST0_star;
R_in_star = inp.R_in_star;
R_cr_star = inp.R_cr_star;
Asp = inp.Asp;
A1i = inp.A1i;
NT = inp.NT;
K = inp.K;

r1 = squeeze(sum(sum(p0.*quant_pred,3),2));        % NTxK
r0 = sum(r1,2);                                    % NTx1
q0 = squeeze(sum(sum(quant_pred,3),2));            % NTxK
pa0 = reshape(psamt_pred,NT,C*2,K);
pa0 = permute(pa0,[1 3 2]);
pa0 = reshape(pa0,NT*K,C*2);
d0 = pa0*storeindex_C2J;
d0 = reshape(d0,NT,K,J);
os0 = sum(prob_vis_storepair(:,A1i),2);            % NTx1
dst0 = sum(prob_vis_storepair.*Asp(:,:,1),2);
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

R0 = R0 - R0_star;
Q0 = sum(Q0 - Q0_star,2);
D0 = sum(squeeze(sum(D0_star - D0,2)).*Ncommon,2);
OS0 = OS0 - OS0_star;
DST0 = DST0 - DST0_star;
R_in = sum(R_in - R_in_star,2);
R_cr = sum(R_cr - R_cr_star,2);
G4 = [R0'; Q0'; D0'; OS0'; DST0'; R_in'; R_cr'];