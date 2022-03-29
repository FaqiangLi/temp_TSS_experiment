function f = TSS_predictions(theta,inp)

C = inp.C;
K = inp.K;
N = inp.NT;
z = inp.z;      % household characteristics
z = reshape(z,N,1,1,inp.L2);
ACi= inp.ACi;
A1i= inp.A1i;
Asp = inp.Asp;
quantity = inp.quantity;
positiveamount = inp.positiveamount;


% [NxCx2x8  NxCx2x8    NxC]
[quant_pred,psamt_pred] = TSS_quantities(theta,inp);

predictions = [];
% 1. quantities
obs = zeros(9,K);
temp1 = reshape(quantity,N*C*2,K);
pred = zeros(9,K);
temp2 = reshape(quant_pred,N*C*2,K);
for s=1:9
    itemp = ACi(:)==s;
    temp = squeeze(sum(temp1(itemp,:),1));
    obs(s,:) = temp;
    temp_2 = squeeze(sum(temp2(itemp,:),1));
    pred(s,:) = temp_2;
end
predictions = [predictions obs(:) pred(:)];

% 2. visitors
obs = zeros(9,K);
temp1 = reshape(positiveamount,N*C*2,K);
pred = zeros(9,K);
temp2 = reshape(psamt_pred,N*C*2,K);
for s=1:9
    itemp = ACi(:)==s;
    temp = squeeze(sum(temp1(itemp,:),1));
    obs(s,:) = temp;
    temp_2 = squeeze(sum(temp2(itemp,:),1));
    pred(s,:) = temp_2;
end
predictions = [predictions obs(:) pred(:)];

% 3. quantities one-stop shoppers
obs = zeros(9,K);
temp1 = reshape(quantity,N*C*2,K);
pred = zeros(9,K);
temp2 = reshape(quant_pred,N*C*2,K);
A1ie = reshape(A1i,1,C);
A1ie = A1ie(ones(N,1),:,ones(2,1));
for s=1:9
    itemp = ACi(:)==s  & A1ie(:)==1;
    temp = squeeze(sum(temp1(itemp,:),1));
    obs(s,:) = temp;
    temp_2 = squeeze(sum(temp2(itemp,:),1));
    pred(s,:) = temp_2;
end
predictions = [predictions obs(:) pred(:)];

% 4. visitors one stop shoppers
obs = zeros(9,K);
temp1 = reshape(positiveamount,N*C*2,K);
pred = zeros(9,K);
temp2 = reshape(psamt_pred,N*C*2,K);
A1ie = reshape(A1i,1,C);
A1ie = A1ie(ones(N,1),:,ones(2,1));
for s=1:9
    itemp = ACi(:)==s  & A1ie(:)==1;
    temp = squeeze(sum(temp1(itemp,:),1));
    obs(s,:) = temp;
    temp_2 = squeeze(sum(temp2(itemp,:),1));
    pred(s,:) = temp_2;
end
predictions = [predictions obs(:) pred(:)];


% 3i. quantities two-stop shoppers
obs = zeros(9,K);
temp1 = reshape(quantity,N*C*2,K);
pred = zeros(9,K);
temp2 = reshape(quant_pred,N*C*2,K);
A1ie = reshape(A1i,1,C);
A1ie = ~logical(A1ie(ones(N,1),:,ones(2,1)));
for s=1:9
    itemp = ACi(:)==s  & A1ie(:)==1;
    temp = squeeze(sum(temp1(itemp,:),1));
    obs(s,:) = temp;
    temp_2 = squeeze(sum(temp2(itemp,:),1));
    pred(s,:) = temp_2;
end
predictions = [predictions obs(:) pred(:)];

% 4i. visitors two-stop shoppers
obs = zeros(9,K);
temp1 = reshape(positiveamount,N*C*2,K);
pred = zeros(9,K);
temp2 = reshape(psamt_pred,N*C*2,K);
A1ie = reshape(A1i,1,C);
A1ie = ~logical(A1ie(ones(N,1),:,ones(2,1)));
for s=1:9
    itemp = ACi(:)==s  & A1ie(:)==1;
    temp = squeeze(sum(temp1(itemp,:),1));
    obs(s,:) = temp;
    temp_2 = squeeze(sum(temp2(itemp,:),1));
    pred(s,:) = temp_2;
end
predictions = [predictions obs(:) pred(:)];

% 5. average per capita income (in pounds)
inc = 1e+1*z(:,:,:,2);
inc = inc(:,ones(C,1),ones(2,1),ones(K,1));

obs = zeros(9,K);
temp1 = reshape(positiveamount.*inc,N*C*2,K);
pred = zeros(9,K);
temp2 = reshape(psamt_pred.*inc,N*C*2,K);
for s=1:9
    itemp = ACi(:)==s;
    temp = squeeze(sum(temp1(itemp,:),1));
    obs(s,:) = temp;
    temp_2 = squeeze(sum(temp2(itemp,:),1));
    pred(s,:) = temp_2;
end
predictions = [predictions obs(:)./predictions(:,3) pred(:)./predictions(:,4)];

% 6. Average household size
hhs = z(:,:,:,1);
hhs = hhs(:,ones(C,1),ones(2,1),ones(K,1));

obs = zeros(9,K);
temp1 = reshape(positiveamount.*hhs,N*C*2,K);
pred = zeros(9,K);
temp2 = reshape(psamt_pred.*hhs,N*C*2,K);
for s=1:9
    itemp = ACi(:)==s;
    temp = squeeze(sum(temp1(itemp,:),1));
    obs(s,:) = temp;
    temp_2 = squeeze(sum(temp2(itemp,:),1));
    pred(s,:) = temp_2;
end
predictions = [predictions obs(:)./predictions(:,3) pred(:)./predictions(:,4)];

% 7. Average distance travelled
dist = Asp(:,:,1);
dist = dist(:,:,ones(2,1),ones(K,1));

obs = zeros(9,K);
temp1 = reshape(positiveamount.*dist,N*C*2,K);
pred = zeros(9,K);
temp2 = reshape(psamt_pred.*dist,N*C*2,K);
for s=1:9
    itemp = ACi(:)==s;
    temp = squeeze(sum(temp1(itemp,:),1));
    obs(s,:) = temp;
    temp_2 = squeeze(sum(temp2(itemp,:),1));
    pred(s,:) = temp_2;
end
predictions = [predictions obs(:)./predictions(:,3) pred(:)./predictions(:,4)];

f = predictions;


