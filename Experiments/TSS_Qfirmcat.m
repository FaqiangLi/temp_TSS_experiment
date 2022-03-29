function [Q,DQ,DQ_pr] = TSS_Qfirmcat(theta,inp,chain,steplength)

storeindex = inp.storeindex;
J = inp.J;
NT = inp.NT;
K = inp.K;
price = inp.price;
t = steplength;
S = size(chain,4);
os = ones(S,1);

chainCJ = reshape(chain,NT,J,1,K,S);
chainCJ = cat(3,chainCJ(:,storeindex(:,1),:,:,:), ...
    chainCJ(:,storeindex(:,2),:,:,:));
price0 = price(:,:,:,ones(S,1)).*chain;

% If calculating derivatives holding c or d&c fixed, save these choices
% here (at observed prices):
if strcmp(inp.fix, 'c') || strcmp(inp.fix, 'dc')
    inp0 = inp;
    inp0.fix = 'save';  %
    q1 = TSS_quantities(theta,inp0);
else
    q1 = TSS_quantities(theta,inp);
end

if any(~isfinite(q1(:)))
    Q = NaN;
    DQ = NaN;
    return
end
DQ = zeros(S*K,S*K);
DQ_pr = DQ;

q = q1(:,:,:,:,os);
% *Note:* store "j" is not the same store for two different consumers, and
% are therefore not in the same chain. On "page s" an individual's
% quantities at a given time in a given store will be set to zero unless
% they are from chain s. Summing over stores within page s we then get the
% sum over stores in chain s.
% NxTxJxKxS. Set to zero the stores on page s not belonging to chain s
q = q.*chainCJ;
clear chainCJ
% total quantity
q0 = reshape(sum(sum(sum(q,1),2),3),K,S);
% quantity purchased by consumers who only visit one firm (whether 1 or 2
% stores)

price_all = reshape(inp.price,NT,J,1,K);
price_all = cat(3,price_all(:,storeindex(:,1),:,:),...
    price_all(:,storeindex(:,2),:,:));

if NT>2000  % do not use parfor loop when data are too big
    for i=1:K*S
        if i==1 || i==K*S || mod(i, 30)==0 % print at 1st and thereafter every 100 iterations
            disp(['Derivative w.r.t to the price of firm-category combination '...
                num2str(i) ' out of ' num2str(K*S)])
        end
        s = ceil(i/K);
        k = i - K*(s-1);
        inp1 = inp;
        % Forwards
        price2 = price0;
        price2(:,:,k,s) = price0(:,:,k,s) + (t/2)*chain(:,:,k,s);
        price2 = sum(price2,4);
        inp1.price = price2;
        q2 = TSS_quantities(theta,inp1);
        q10 = zeros(K,S);
        q10_pr = zeros(K,S);
        for s2=1:S
            chainCJs = reshape(chain(:,:,:,s2),NT,J,1,K);
            chainCJs = cat(3,chainCJs(:,storeindex(:,1),:,:), ...
                chainCJs(:,storeindex(:,2),:,:));
            q_s = q2.*chainCJs;
            q10(:,s2) = reshape(sum(sum(sum(q_s,1),2),3),K,1);
            q10_pr(:,s2) = reshape(sum(sum(sum(q_s.*price_all,1),2),3),K,1);
        end
        % Backwards
        price2 = price0;
        price2(:,:,k,s) = price0(:,:,k,s) - (t/2)*chain(:,:,k,s);
        price2 = sum(price2,4);
        inp1.price = price2;
        q2 = TSS_quantities(theta,inp1);
        q00 = zeros(K,S);
        q00_pr = zeros(K,S);
        for s2=1:S
            chainCJs = reshape(chain(:,:,:,s2),NT,J,1,K);
            chainCJs = cat(3,chainCJs(:,storeindex(:,1),:,:), ...
                chainCJs(:,storeindex(:,2),:,:));
            q_s = q2.*chainCJs;
            q00(:,s2) = reshape(sum(sum(sum(q_s,1),2),3),K,1);
            q00_pr(:,s2) = reshape(sum(sum(sum(q_s.*price_all,1),2),3),K,1);
        end
        
        Dq0 = (q10-q00)/t;
        DQ(:,i) = Dq0(:);
        Dq0 = (q10_pr-q00_pr)/t;
        DQ_pr(:,i) = Dq0(:);
    end
else
    %disp('TSS_Qfirmcat line 89: using for loop. Change to parfor when more memory available!')
    disp('Calculating demand derivatives.')
    parfor i=1:K*S
        if i==1 || i==K*S || mod(i, 30)==0 % print at 1st and thereafter every 100 iterations
            disp(['Derivative w.r.t to the price of firm-category combination '...
                num2str(i) ' out of ' num2str(K*S)])
        end
        s = ceil(i/K);
        k = i - K*(s-1);
        inp1 = inp;
        % Forwards
        price2 = price0;
        price2(:,:,k,s) = price0(:,:,k,s) + (t/2)*chain(:,:,k,s);
        price2 = sum(price2,4);
        inp1.price = price2;
        q2 = TSS_quantities(theta,inp1);
        q10 = zeros(K,S);
        q10_pr = zeros(K,S);
        for s2=1:S
            chainCJs = reshape(chain(:,:,:,s2),NT,J,1,K);
            chainCJs = cat(3,chainCJs(:,storeindex(:,1),:,:), ...
                chainCJs(:,storeindex(:,2),:,:));
            q_s = q2.*chainCJs;
            q10(:,s2) = reshape(sum(sum(sum(q_s,1),2),3),K,1);
            q10_pr(:,s2) = reshape(sum(sum(sum(q_s.*price_all,1),2),3),K,1);
        end
        % Backwards
        price2 = price0;
        price2(:,:,k,s) = price0(:,:,k,s) - (t/2)*chain(:,:,k,s);
        price2 = sum(price2,4);
        inp1.price = price2;
        q2 = TSS_quantities(theta,inp1);
        q00 = zeros(K,S);
        q00_pr = zeros(K,S);
        for s2=1:S
            chainCJs = reshape(chain(:,:,:,s2),NT,J,1,K);
            chainCJs = cat(3,chainCJs(:,storeindex(:,1),:,:), ...
                chainCJs(:,storeindex(:,2),:,:));
            q_s = q2.*chainCJs;
            q00(:,s2) = reshape(sum(sum(sum(q_s,1),2),3),K,1);
            q00_pr(:,s2) = reshape(sum(sum(sum(q_s.*price_all,1),2),3),K,1);
        end
        
        Dq0 = (q10-q00)/t;
        DQ(:,i) = Dq0(:);
        Dq0 = (q10_pr-q00_pr)/t;
        DQ_pr(:,i) = Dq0(:);
    end
end
Q =  q0(:);
