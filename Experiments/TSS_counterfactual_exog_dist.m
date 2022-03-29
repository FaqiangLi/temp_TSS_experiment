function f = TSS_counterfactual_exog_dist(theta,inp,mcost,pr_w,steplength)

print = inp.print;
K = inp.K;
S = 16;
chain = inp.chain2;
onestop = inp.A1i;
price = inp.price;
% shopping costs:
L0 = 34;
dist_i = 10;
pr = reshape(pr_w,8,16);

nugamma = squeeze(inp.nugamma);
gm_sg = theta(L0+3:L0+4);
shop_cost =  gm_sg(1)*nugamma(:,1) + dist_i.*(gm_sg(2)*nugamma(:,2)); % enters utility with a minus in front

[~,~,P] = TSS_quantities(theta,inp);

%criterion = 'shop_cost';
criterion = 'prob_1ss';
%cutoff = sum(sum(P(:,~onestop)))/2000;
if strcmp(criterion,'shop_cost') % shopping cost
    cutoff = 0.5;
    ix_A = shop_cost >= quantile(shop_cost,cutoff);
    ix_B = shop_cost < quantile(shop_cost,cutoff);
    ix_g = [ix_A ix_B];
elseif strcmp(criterion,'prob_1ss')% probability of visiting one store only
    ix_A = sum(P(:,onestop),2)>=0.5;
    ix_B = sum(P(:,onestop),2)<0.5;
    ix_g = [ix_A ix_B];
end

% A1s = 100*sum(sum(P(ix_A,onestop)))/sum(ix_A);
% A2s = 100*sum(sum(P(ix_A,~onestop)))/sum(ix_A);
% B1s = 100*sum(sum(P(ix_B,onestop)))/sum(ix_B);
% B2s = 100*sum(sum(P(ix_B,~onestop)))/sum(ix_B);
% TypeShare = [A1s B1s; A2s B2s]

A1s = 100*sum(sum(P(ix_A,onestop)))/sum(sum(P(:,onestop)));
A2s = 100*sum(sum(P(ix_A,~onestop)))/sum(sum(P(:,~onestop)));
B1s = 100*sum(sum(P(ix_B,onestop)))/sum(sum(P(:,onestop)));
B2s = 100*sum(sum(P(ix_B,~onestop)))/sum(sum(P(:,~onestop)));
TypeShare = [A1s B1s; A2s B2s]


%_______ Price increase for one category, all consumers ___________________
% A1: No. of visitors for the category in each group before and after price change.
% B: No. of visitors for each (all 8) category in each group before and after
% C3: Total profit (all cats) before and after; revenue for the category whose price changes.

T = 16*8; % number of combinations of firm and category
A1 = zeros(3,T);
B = zeros(3,T);
C3 = zeros(3,T);
%for cat = [3 5 6 7];
cat_vec = repmat(1:8,1,16);
frm_vec = kron([2 8 13 15 1 7 10 3 4 5 6 9 11 12 14 16],ones(1,8));
parfor t=1:T
    frm = frm_vec(t);
    cat = cat_vec(t);
%     disp(['Firm ' num2str(frm)])
%     disp(['category ' num2str(cat)])
    
    %______ Before _____________
    price0 = price(:,:,:,ones(S,1)).*chain;
    price0(:,:,cat,frm) = price0(:,:,cat,frm) - (steplength/2)*chain(:,:,cat,frm);
    price0 = sum(price0,4);
    inp0 = inp;
    inp0.price = price0;
    [visits_b,profit_b] = TSS_profit_visits(theta,inp0,mcost,ix_g,frm);
    
    %______ After ______________
    price0 = price(:,:,:,ones(S,1)).*chain;
    price0(:,:,cat,frm) = price0(:,:,cat,frm) + (steplength/2)*chain(:,:,cat,frm);
    price0 = sum(price0,4);
    inp0 = inp;
    inp0.price = price0;
    [visits_a,profit_a] = TSS_profit_visits(theta,inp0,mcost,ix_g,frm);
    %___________________________
%     A1(:,t) = [(visits_a(cat,1)-visits_b(cat,1));...
%         (visits_a(cat,2)-visits_b(cat,2))];
    A1(:,t) = [(visits_a(cat,1)-visits_b(cat,1))/visits_b(cat,1);...
        (visits_a(cat,2)-visits_b(cat,2)) / visits_b(cat,2);...
        ( sum(visits_a(cat,:)) - sum(visits_b(cat,:)) ) / sum(visits_b(cat,:)) ;...
        ] /(steplength/pr(cat,frm));
    d_vis = visits_a - visits_b;
    ix = true(K,1);
    ix(cat) = false;
    other = sum(d_vis(ix,:));
    ratio = other ./ d_vis(cat,:);
    B(:,t) = [ratio(1); ratio(2); sum(other) / sum(d_vis(cat,:)) ];
    %prof_b = sum(profit_b,1); 
    %d_prof = sum(profit_a - profit_b);
    d_prof = sum(profit_a - profit_b);

%     disp(['Total change in profit is             ' num2str(sum(d_prof)) ])
%     disp(['Change in profit for prob. 1ss type   ' num2str(d_prof(1)) ])
%     disp(['Change in profit for prob. 2ss type   ' num2str(d_prof(2)) ])
%     disp(['Total change in cat. visitors is      ' num2str(sum(d_vis(cat,:))) ])
%     disp(['Change in cat.vis. for prob. 1ss type ' num2str((visits_a(cat,1)-visits_b(cat,1))) ])
%     disp(['Change in cat.vis. for prob. 2ss type ' num2str((visits_a(cat,2)-visits_b(cat,2))) ])
    C3(:,t) = [d_prof(1); d_prof(2); d_prof(1) + d_prof(2)] / steplength;
end
Table = [A1; B; C3];
%          Big4                      Disc                       All
f = [[median(Table(:,1:4*8),2) median(Table(:,4*8+1:7*8),2) median(Table,2)];...
     [mean(Table(7,1:4*8)<0,2) mean(Table(7,4*8+1:7*8)<0,2) mean(Table(7,:)<0,2)];...
     [mean(Table(8,1:4*8)>0,2) mean(Table(8,4*8+1:7*8)>0,2) mean(Table(8,:)>0,2)]]
 
if print
    FileName = ['Output/' 'Exog_type_counterf','_',criterion,'_',datestr(now, 'ddmm_yyyy_HHMM'),'.xlsx'];
    xlswrite(FileName,f,1,'G4:I14')
    xlswrite(FileName,{'Big 4','Disc', 'All'},1,'G3:I3')
    xlswrite(FileName,{'Elasticity of 1ss type category visits w.r.t. price'},1,'A4')
    xlswrite(FileName,{'Elasticity of 2ss type category visits w.r.t. price'},1,'A5')
    xlswrite(FileName,{'Elasticity of all category visits w.r.t. price'},1,'A6')
    xlswrite(FileName,{'Diversion ratio for 1ss type visits'},1,'A7')
    xlswrite(FileName,{'Diversion ratio for 2ss type visits'},1,'A8')
    xlswrite(FileName,{'Diversion ratio for all visits'},1,'A9')
    xlswrite(FileName,{'Derivative of firm profit from 1ss w.r.t. price'},1,'A10')
    xlswrite(FileName,{'Derivative of firm profit from 2ss w.r.t. price'},1,'A11')
    xlswrite(FileName,{'Derivative of firm profit from both types w.r.t. price'},1,'A12')
    xlswrite(FileName,{'Proportion where the derivative of 1ss-type profit <0'},1,'A13')
    xlswrite(FileName,{'Proportion where the derivative of 2ss-type profit >0'},1,'A14')
    
    xlswrite(FileName,{['Finite-difference derivative step length (forward+backward): ' num2str(steplength)]},1,'I1')
    
    if strcmp(criterion,'shop_cost') % shopping cost
        xlswrite(FileName,{'Partition based on shopping cost'},1,'A1')
    elseif strcmp(criterion,'prob_1ss')
        xlswrite(FileName,{'Partition based on probability of one-stop shopping being >0.5 or <0.5'},1,'A1')
    end
end
 