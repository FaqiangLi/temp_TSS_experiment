function [] = simulate_data_mainscript(parameter_customized,total_households,input_dir,output_dir,seed)
% parameter_customized : parameters used to simulate the model
% total_households     : the number of the candidate households to whom we
%                        use the model to give values to their endogeneous
%                        variables.  
% input_dir            : directory that stores the exogeneous variables
% output_dir           : directory taat stores both the generated
%                        endogeneous variables and the exogeneous variables
% seed                 : seed for random number generator    

%% Run to read and reformat the *exogeneous variables*
% The script is modified according to TSS_input and TSS_quantities

run clip_simulate_choice_probability.m

%% Simulate choice data

% Generate probabilistic behaviors
[quant_pred,psamt_pred,prob_vis_storepair,Q,jbest]=simulate_choice(theta,inp);

% Draw determinnistic behaviors from the above data (store choice, volumn choice)
device = rand(NT,1); % NT*1 random uniform number
cumsum_prob = cumsum(prob_vis_storepair,2);
choice_draw = sum(device > cumsum_prob,2) + 1; % used later to draw choices

%% Draw and Select ``healthy" households
% N_select denotes number of healthy households for estimation.
% The reason why I need to check healthiness here is that, I
% want the following code giving me those households that have full periods
% of purchase with positive amount in some category given all
% the pair of store choice.

N_select = total_households;

if size(Q)~=[NT C K]
    Q=permute(Q,[1 3 2]) % We want size(Q)=NT*C*K
end

% Now realize actual quantity by using choice_draw
Q_realized=zeros(NT,K);
for i=1:NT
    cc=choice_draw(i,1);
    Q_realized(i,:)=Q(i,cc,:);
end

% Tag households having three periods of non-zero purchase
positive_purchase=sum(Q_realized,2)>0;
id_pos = zeros(N,1);
for i=1:N
    all_three_purchase = positive_purchase(i,1)+positive_purchase(i+N,1)+positive_purchase(i+(T-1)*N,1);
    index = all_three_purchase==3;
    id_pos(i,1) = index;
end
check_pos_ix = length(find(id_pos));

% show in the diary
N_select
check_pos_ix
assert( N_select<=check_pos_ix , 'Not enough households having positive full periods purchase')

% Indices to select health households' data
position = find(id_pos);
position_select = randsample(position,N_select);
position_select_sorted = sort(position_select,1);
ix_n   = position_select_sorted(1:N_select);  
ix_nt = [ix_n;ix_n+N;ix_n+2*N];
ix_ntj_temp = zeros(N*T*J,1);
for ii = 1:N_select
    which_household = ix_n(ii);
    for tt = 1:T
        starting=(tt-1)*N*J + (which_household-1)*J+1;
        ix_ntj_temp(starting:1:starting+J-1,1) =1;
    end
end
ix_ntj=find(ix_ntj_temp);

% Allocate quantity to store A and store B
jbest=permute(jbest,[1 3 2]); % originally jbest is NT*C*K
jbest_realized=zeros(NT,K);
for i=1:NT
    cc=choice_draw(i,1);
    jbest_realized(i,:)=jbest(i,:,cc);
end
Q_storeAA=Q_realized.*jbest_realized;
Q_storeBB=Q_realized.*(1-jbest_realized);

% small test, making sure Q_storeAA and Q_storeBB are reasonable
Q_realized_sum=Q_storeAA+Q_storeBB;
small_test2=sum(sum(Q_realized_sum(ix_nt,:),2)>0,1)==N_select*T;

% Preparing matrix for expenditure related matrics
purchase_quant=zeros(NT,JK); % to construct cat_spend.csv
num_store=zeros(NT,1); % to construct sameday.csv component
which_store=zeros(NT,J); % to construct store_visit.csv
for i=1:NT
    c=choice_draw(i,1);
    temp=storeindex(c,:);
    AA=temp(1,1);
    BB=temp(1,2);
    num_store(i,:)= 1*(AA==BB)+2*(AA~=BB);
    which_store(i,AA)=1;
    which_store(i,BB)=1;
    
    if num_store(i,:)==1
        purchase_quant(i,K*(AA-1)+1:K*AA)=Q_storeAA(i,:)+Q_storeBB(i,:);
    else
        purchase_quant(i,K*(AA-1)+1:K*AA)=Q_storeAA(i,:);
        purchase_quant(i,K*(BB-1)+1:K*BB)=Q_storeBB(i,:);
    end
end

% small test, making sure no purchase is totally empty
temp_4test=sum(purchase_quant(ix_nt,:),2)>0;
small_test3=sum(temp_4test,1)==N_select*T;
if small_test3
    disp('No purchase records have all zero category')
end

% For sameday.csv
sameday=zeros(NT,5);
sameday(:,4)=num_store;
sameday(:,5)=0;

% For store_visit.csv , it is from which_store constructed above
store_visit=zeros(NT,3+J);
store_visit(:,4:end)=which_store;

% For cat_USI.csv 
cat_USIs=zeros(NT,3+JK);
% Generate choice set (it won't matter at this stage)
% they are just store index
choice_set=zeros(N,J);
for i=1:N
    choice_set(i,:)=randperm(1000,J); % store number is unique
end
% So 1000-1999 are the store index
choice_set=choice_set+999;
usi_temp=repmat(choice_set,T,1);
usi_temp=kron(usi_temp,ones(1,K));
% Only if positive amount is spent, we have usi non-zero
cat_USIs(:,4:end)=(purchase_quant>0).*usi_temp;

% For cat_spends.csv
cat_spends=zeros(NT,3+JK);
cat_spends(:,4:end)=purchase_quant.*cat_prices(:,4:end);

% Check isnan
test2 = zeros(9,1);
test2(1) = sum(isnan(cat_iv(:)));
test2(2) = sum(isnan(cat_prices(:)));
test2(3) = sum(isnan(cat_spends(:)));
test2(4) = sum(isnan(cat_USIs(:)));
test2(5) = sum(isnan(NTJfirmnum2(:)));
test2(6) = sum(isnan(NTjchars2(:)));
test2(7) = sum(isnan(Hhchars2(:)));
test2(8) = sum(isnan(store_visit(:)));
test2(9) = sum(isnan(sameday(:)));
if sum(test2)~=0
    disp('Warning: some data array contains NaN entries.')
    %return
end

% Reselect all data picked to be sent to the structural estimation
cat_prices = cat_prices(ix_nt,:);
NTjchars2 = NTjchars2(ix_ntj,:);
NTJfirmnum2 = NTJfirmnum2(ix_nt,:);
Hhchars2 = Hhchars2(ix_n,:);
cat_iv = cat_iv(ix_nt,:);
sameday=sameday(ix_nt,:);
store_visit=store_visit(ix_nt,:);
cat_USIs=cat_USIs(ix_nt,:);
cat_spends=cat_spends(ix_nt,:);

% supplementing identifiers for choice data.
% note that ix_nt is for selection and id_nt is for refilling real id. 
id_nt = cat_prices(:,1:3);
sameday(:,1:3)=id_nt;
store_visit(:,1:3)=id_nt;
cat_USIs(:,1:3)=id_nt;
cat_spends(:,1:3)=id_nt;


final_test=sum(sum(cat_spends(:,4:end),2)>0,1)==N_select*T;


%% write and save in txt
dir=ITAM_output_dir;
cd(dir)

writematrix(sameday)
writematrix(store_visit)
writematrix(Hhchars2)
writematrix(NTjchars2)
writematrix(NTJfirmnum2)
writematrix(cat_prices)
writematrix(cat_iv)
writematrix(cat_USIs)
writematrix(cat_spends)


display(['The data have been simualated. The seed is', num2str(seed_in_use)])
display('Now one can go back to working script to conduct the subsequent task.')

clear

return
