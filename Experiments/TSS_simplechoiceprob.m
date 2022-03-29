function f = TSS_simplechoiceprob(storesize,distance,beta)

% returns the NxJ matrix of choice probabilities: (n,j) gives the
% probability of choosing alternative j from among the J alternatives
% available to household n. 

J = size(storesize,2);
% N0xJxK
u = beta(1)*storesize + beta(2)*distance;
eu = exp(u);
eu_denom = sum(eu,2);
eu_denom = eu_denom(:,ones(J,1),:);
f = eu./eu_denom;
