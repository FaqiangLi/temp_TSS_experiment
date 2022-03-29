function [qstar,vstar] = QP_nonneg(A,B)
% QP_NONNEG quadratic programming with nonnegativity constraints.
% Calculates the solutions and optimal values to the
% problem max_{q>=0} u(q) where  u(q) = q'b - 0.5q'Aq. The matrix B has b
% corresponding to different problems as its columns.
% If B is KxN, qstar is KxN and vstar is 1xN.
% Notes: using linsolve or Cramer's rule (analytical solution) for k=1,2,3
% instead of mldivide (\) does not increase speed (and may be numerically
% inaccurate)

[K,N] = size(B);
R = 2^K;
q = zeros(K,N,R);
r = 1;
for k=1:K
    % combinations of k's with interior (positive) solution
    combination = combnk(1:K,k);
    L = size(combination,1);
    for l=1:L
        ix = combination(l,:);
        q(ix,:,r) = A(ix,ix)\B(ix,:); 
        r = r+1;
    end
end
negatives_r = min(q,[],1) < 0 ;
v = zeros(1,N,R);
for r=1:R
    qr = q(:,:,r);
    v(:,:,r) = sum((B - 0.5*A*qr).*qr,1);
end
% if the optimum in case r violates a constraint, set v(r) = -1e30, to
% ensure that r does not become the chosen outcome.
v(negatives_r) = -1e30;
% highest utility (v*) and best outcome (r*)
[vstar,rstar] = max(v,[],3);
ix = sub2ind([N,R],1:N,rstar);
qstar = q(:,ix);

%__________ test: comparison with quadprog ________________________________
% options = optimoptions('quadprog','Algorithm','trust-region-reflective','Display','off','TolX',1e-10,'TolFun',1e-10);
% for n=1:N
%     startval = qstar(:,n) + 0.1*randn(K,1);
%     [x,fval] = quadprog(A,-B(:,n),[],[],[],[],zeros(K,1),[],startval,options);
%     dist = max(abs([-fval; x]-[vstar(:,n); qstar(:,n)]));
%     if  dist > 1e-5
%         disp('warning: quadprog gives different solution from QP_nonneg')
%         return
%     end
% end
%__________________________________________________________________________




