function [markups,el] = TSS_multipledraws_markups(theta,inp,S,~,steplength)

% TSS_multipledraws_markups: returns a (K*S)x7 matrix
% [pr mcost markup percmu mcostDeleg markupDeleg percmuDeleg]
% of markups and perc.markups for supermarket pricing and delegated pricing
% S=9 when firms are aggregated and S=16 otherwise.

K = inp.K;
omega = kron(eye(S),ones(K,K)); % ownership matrix
if S==9
    chain = inp.chain; 
elseif S==16
    chain = inp.chain2;
end

[Q,DQ,DQ_pr] = TSS_Qfirmcat(theta,inp,chain,steplength);
%save eletc_testing Q DQ DQ_pr
% disp('warning: loading saved derivatives')
% load eletc_testing

Omega = omega.*(DQ');
omegaDeleg = eye(S*K); % 'ownership matrix' for single-category 'firms'
OmegaDeleg = omegaDeleg.*(DQ');
sO = size(Omega,1);
sOd = size(OmegaDeleg,1);
sQ = size(Q,1);
if sO == sQ && sOd == sQ
    mc = Omega\(Q + sum((omega.*DQ_pr'),2));
    mc_deleg = OmegaDeleg\(Q + sum(omegaDeleg.*DQ_pr',2));
else
    disp('Error: dimensions of Omega and Q do not match') % not sure why this error might occur, but it has.
    markups = NaN;
    return
end

p_w = diag(DQ_pr')./diag(DQ');

markup = p_w - mc;
markupDeleg = p_w - mc_deleg;
percmu = 100*markup./p_w;
percmuDeleg = 100*markupDeleg./p_w;

markups = [p_w mc markup percmu mc_deleg markupDeleg percmuDeleg] ;
if nargout>1 
    p = p_w(:,ones(K*S,1))';
    q = Q(:,ones(K*S,1));
    el = p.*DQ./q;
end
%save mcost_testing mcost
