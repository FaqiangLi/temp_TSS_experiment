function [obsvis,vis,obsrev,revenue,table] = TSS_pairs(theta,inp)

% obsvis: 16x16 matrix of giving the number of visitors to each combination
% of firms
% vis: same as obsvis, but the model's predictions
% obsrev: 16x16 matrix like obsvis, but giving revenue generated rather
% than number of consumers
% revenue: same as obsrev, but the model's predictions
% table: 32x10 matrix giving, for each firm, the observed (above) and
% predicted (below) values for 1SS (left) and 2SS (right) of:
% (a) Revenue, % of total		
% (b) visits, % of total no. of households		
% (c) Average household size		
% (d) Average household income		
% (e) Average distance travelled to stores	
									
C = inp.C;
K = inp.K;
N = inp.NT;
S = 16;
z = inp.z;
z = reshape(z,N,1,1,inp.L2);
Asp = inp.Asp;
ACi = inp.ACi2;
storeindex = inp.storeindex;
storepairvisited = inp.storepairvisited;
quantity = inp.quantity;

[quant_pred,~,P]  = TSS_quantities(theta,inp); 

pr = inp.price;
prices = zeros(N,C,2,K);
prices(:,:,1,:) = pr(:,storeindex(:,1),:);
prices(:,:,2,:) = pr(:,storeindex(:,2),:);

% Revenue generated in each store in each store pair from each consumer
% (summed across categories)
% (NxCx2)
rev_o = sum(prices.*quantity,4);
rev_p = sum(prices.*quant_pred,4);

% householdsize
hhs = z(:,:,:,1);
% income
inc = 1e+1*z(:,:,:,2);
% distance to store pair
dist = Asp(:,:,1);

vis = zeros(S,S);
obsvis = vis;
revenue = vis;
obsrev = vis;

table = zeros(2*S,2*5);
% Firms - 1 store; Firms - 2 store
for s=1:S
    onefirm = ACi(:,:,1)==s & ACi(:,:,2)==s;
    twofirm =(ACi(:,:,1)==s & ACi(:,:,2)~=s)|...
             (ACi(:,:,1)~=s & ACi(:,:,2)==s);

    % revenue for firm s from consumer n in store pair c
    % NxC
    revo = (rev_o(:,:,1).*(ACi(:,:,1)==s))+(rev_o(:,:,2).*(ACi(:,:,2)==s));
    revp = (rev_p(:,:,1).*(ACi(:,:,1)==s))+(rev_p(:,:,2).*(ACi(:,:,2)==s));     
    
    rev1o = sum(sum(revo.*onefirm,2),1);
    rev1p = sum(sum(revp.*onefirm,2),1);
    rev2o = sum(sum(revo.*twofirm,2),1);
    rev2p = sum(sum(revp.*twofirm,2),1);
    
    table(2*s-1:2*s,1:2) = [ [rev1o;rev1p] [rev2o;rev2p] ];
    
    % visitors
    vi1o = sum(sum(storepairvisited.*onefirm,2),1);
    vi1p = sum(sum(P.*onefirm,2),1);
    vi2o = sum(sum(storepairvisited.*twofirm,2),1);
    vi2p = sum(sum(P.*twofirm,2),1);
    
    table(2*s-1:2*s,3:4) = [ [vi1o;vi1p] [vi2o;vi2p] ];
    
    for z=1:S
        temp = (ACi(:,:,1)==s &  ACi(:,:,2)==z) |  (ACi(:,:,1)==z &  ACi(:,:,2)==s); 
        obsvis(s,z) = sum(sum(storepairvisited.*temp));
        vis(s,z) = sum(sum(P.*temp));
        obsrev(s,z) = sum(sum(revo.*temp));
        revenue(s,z) = sum(sum(revp.*temp));
    end
    
    do1 = sum(sum(storepairvisited.*onefirm,2),1);
    dp1 = sum(sum(P.*onefirm,2),1);
    do2 = sum(sum(storepairvisited.*twofirm,2),1);
    dp2 = sum(sum(P.*twofirm,2),1);
    % average householdsize
    hh1o = sum(hhs.*sum(storepairvisited.*onefirm,2),1)/do1;
    hh1p = sum(hhs.*sum(P.*onefirm,2),1)/dp1;
    hh2o = sum(hhs.*sum(storepairvisited.*twofirm,2),1)/do2;
    hh2p = sum(hhs.*sum(P.*twofirm,2),1)/dp2;
    
    table(2*s-1:2*s,5:6) = [ [hh1o;hh1p] [hh2o;hh2p] ];
    % average income
    in1o = sum(inc.*sum(storepairvisited.*onefirm,2),1)/do1;
    in1p = sum(inc.*sum(P.*onefirm,2),1)/dp1;
    in2o = sum(inc.*sum(storepairvisited.*twofirm,2),1)/do2;
    in2p = sum(inc.*sum(P.*twofirm,2),1)/dp2;
    
    table(2*s-1:2*s,7:8) = [ [in1o;in1p] [in2o;in2p] ];
    % average distance travelled
    dst1o = sum(sum(dist.*storepairvisited.*onefirm,2),1)/do1;     
    dst1p = sum(sum(dist.*P.*onefirm,2),1)/dp1;      
    dst2o = sum(sum(dist.*storepairvisited.*twofirm,2),1)/do2;     
    dst2p = sum(sum(dist.*P.*twofirm,2),1)/dp2;
   
    table(2*s-1:2*s,9:10) = [ [dst1o;dst1p] [dst2o;dst2p] ];
end
odd = 1:2:2*S;
even= 2:2:2*S;
sm1 = sum(sum(table(odd,1:2)));
sm2 = sum(sum(table(even,1:2)));
table(odd,1:2) = 100*table(odd,1:2)/sm1;
table(even,1:2) = 100*table(even,1:2)/sm2;
table(:,3:4) = 100*table(:,3:4)/N;

