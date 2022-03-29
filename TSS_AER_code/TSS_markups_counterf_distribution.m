function [markups,markups_deleg,markups_etc_full] = TSS_markups_counterf_distribution(theta,acov,inp,S,NS,steplength)

% TSS_markups_distribution calculates a distribution of median markups (and
% related quantities such as marginal cost), elasticities and results from
% counterfactual experiments, based on NS draws from the the asymptotic
% distribution of the estimator of theta, and with the median being across
% different weeks of the data.

% __________ Input arguments ______________________________________________
% theta:            estimated parameters (Dx1)
% acov:             estimated covariance matrix of estimator (DxD)
% inp:              input structure
% S:                9 or 16. The number of firms, depending on level of
%                   aggregation
% NS:               number of draws from the distribution of the estimator
% steplength:       step length for finite difference derivatives of demand
%                   with respect to price
% TS:               number of weeks


week = 78;
D = inp.D;
NT = inp.NT;
J = inp.J;
K = inp.K;
storeindex = inp.storeindex;
price = inp.price;
if S==9
    ACi = inp.ACi;
elseif S==16
    ACi = inp.ACi2;
end

quantity = inp.quantity;
pr = reshape(price,NT,J,1,K);
pr = cat(3,pr(:,storeindex(:,1),:,:),pr(:,storeindex(:,2),:,:));
revenue = pr.*quantity;
revenue_milk = revenue(:,:,:,8);
revenue_all = revenue;
rev_s_m = zeros(S,1);
rev_s_a = zeros(S,K);
for s=1:S
    temp = (ACi==s).*revenue_milk;
    rev_s_m(s) = sum(temp(:));
    ACi_K = ACi(:,:,:,ones(K,1));
    temp = (ACi_K==s).*revenue_all;
    rev_s_a(s,:) = sum(sum(sum(temp,1),2),3);
end
% for use in calculating a weighted median:
min_rev_m = min(rev_s_m);
if min_rev_m<10
    rev_s_m = round(rev_s_m);
elseif min_rev_m<100
    rev_s_m = round(rev_s_m/10);
end

min_rev_a = min(rev_s_a(:));
if min_rev_a<100
    rev_s_a = round(rev_s_a/10);
elseif min_rev_a<1000
    rev_s_a = round(rev_s_a/100);
end

chol_acov = chol(acov,'lower');
R = 1;
disp('Generating weekly data.')
inp1 = TSS_weekly_input_struct(inp,week);
inp1.print = false;
%ix_milk = K:K:S*K;
markups = [];
markups_deleg = [];
counterf = [];
markups_etc_full = [];
counter = 0;
tic
for ns=1:NS
    counter = counter + 1;
    disp(['Random draw number ' num2str(ns) ' of ' num2str(NS)])
    done = false;       % in case of numerical problems in calculation of markups, a new draw is attempted
    while ~done
        rng('shuffle')
        theta_ns = theta + chol_acov*randn(D,1);       % draw number ns from the distribution of the estimator
        
%         disp('NO DRAW - USING EST!') 
%         theta_ns = theta
        
        
        mu = TSS_multipledraws_markups(theta_ns,inp1,S,R,steplength);
        markups_ns = mu;
        pr_ns = mu(:,1); % 'weighted' prices
        mcost_ns = mu(:,2);
        %load mcost
        counterf_ns = TSS_counterfactual_exog_dist(theta_ns,inp1,mcost_ns,pr_ns,steplength);
        
        if all(isfinite(markups_ns(:)))
            done = true;
            mu_a = markups_ns(:,4);
            mu_m = markups_ns(8:8:end,4);
            markups_all = median(mu_a);
            markups_milk = median(mu_m);
            mu_m_exp = [];
            mu_a_exp = [];
            for s=1:S
                mu_m_exp = [mu_m_exp; repmat(mu_m(s),rev_s_m(s),1)];
                for k=1:K
                    mu_a_exp = [mu_a_exp; repmat(mu_a((s-1)*K+k),rev_s_a(s,k),1)];
                end
            end
            markups_all_wgt = median(mu_a_exp);
            markups_milk_wgt = median(mu_m_exp);
            
            mu_a = markups_ns(:,7);
            mu_m = markups_ns(8:8:end,7);
            markups_all_deleg = median(mu_a);
            markups_milk_deleg = median(mu_m);
            mu_m_exp = [];
            mu_a_exp = [];
            for s=1:S
                mu_m_exp = [mu_m_exp; repmat(mu_m(s),rev_s_m(s),1)];
                for k=1:K
                    mu_a_exp = [mu_a_exp; repmat(mu_a((s-1)*K+k),rev_s_a(s,k),1)];
                end
            end
            markups_all_wgt_deleg = median(mu_a_exp);
            markups_milk_wgt_deleg = median(mu_m_exp);
            
%             Big4 = [2 8 13 15];
%             Disc = [1 7 10];
            ix_b4 = [(2-1)*8+1:2*8 (8-1)*8+1:8*8 (13-1)*8+1:13*8 (15-1)*8+1:15*8];
            ix_ds = [(1-1)*8+1:1*8 (7-1)*8+1:7*8 (10-1)*8+1:10*8];
            
            pA = markups_ns(:,4)/100;
            pB = markups_ns(:,7)/100;
            pC = pB-pA;
            
            markup_means_ns = [ mean(pA);...
                                mean(pA(ix_b4)) ;...
                                mean(pA(ix_ds)) ;...
                                mean(pB);...
                                mean(pB(ix_b4));...
                                mean(pB(ix_ds));...
                                mean(pC);...
                                mean(pC(ix_b4));...
                                mean(pC(ix_ds))];
             
        else
            disp('A draw led to numerical problems - a new draw is attempted.')
        end
    end
    markups = [markups; [markups_milk markups_milk_wgt markups_all markups_all_wgt]];
    markups_deleg = [markups_deleg; [markups_milk_deleg markups_milk_wgt_deleg markups_all_deleg markups_all_wgt_deleg]];
    markups_etc_full = cat(3,markups_etc_full,markup_means_ns);
    counterf = cat(3,counterf,counterf_ns);
    if mod(ns, 10) == 0  % Save for every 10th draw. 
       time = toc;
       disp(['Time per draw is approximately ' num2str(time/counter) ' seconds.'])
       counter = 0;
       FileName = ['Output/','Markups_counterf_distribution',datestr(now, 'ddmm_yyyy_HHMM'),'.mat'];
       save(FileName, 'markups', 'markups_deleg','markups_etc_full','counterf')
       tic
    end
end

