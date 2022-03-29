function acov = TSS_acov(theta,inp,inp2)

N = inp.N;
W = inp.W;
LL = inp.LL;
D = length(theta);
instr_ix = inp.instr_ix;
new_moments = inp.new_moments;

%_____________derivatives of moments_______________________________________
stepfactor = 5e-4;
dG = zeros(LL,D);
parfor pd=1:D
    if pd==1 | pd==D |  mod(pd, 10)==0 % print 1st and thereafter every 10 iterations
        disp(['derivative w.r.t. parameter number '...
            num2str(pd) ' out of ' num2str(D)])
    end
    dGtemp = zeros(LL,2);
    for t=1:2
        theta2 = theta;
        step = ((-1)^t) * stepfactor * abs(theta2(pd));
        theta2(pd) = theta2(pd) + step;
        [G1,G2,G3,G4] = TSS_moments(theta2,inp);
        g = [G1; G2; G3];
        g = g(instr_ix);
        if new_moments
            g = [g; sum(G4,2)];
        end
        dGtemp(:,t) = g;
    end
    dG(:,pd) = (dGtemp(:,2)-dGtemp(:,1))./(2*step);
    if length(nonzeros(dG(:,pd))) == 0
        test = 0;
    end
end
save dG_testing dG
% disp('Warning: loading saved derivatives of moments instead of computing new ones.')
% load dG_testing

%_______________Covariance matrices of moments_____________________________
Lambda = TSS_W_second_stage(theta,inp,inp2);
% See Wooldridge, Cross Section and Panel, 2nd ed., p. 527 for asymptotic
% covariance matrix, and Stern (1997): Simulation-Based Estimation, p.
% 2029-2030, eq. 3.18, for correcting for simulation. 
dG = dG/N;                          % Wooldridge, eq. (14.15), p. 528
A = dG'*W*dG;                       % Wooldridge, eq. (14.14), p. 527
iA = inv(A);
B = dG'*W*Lambda*W*dG;              % Wooldridge, eq. (14.14), p. 527
acov = iA*B*iA/N;                   % Wooldridge, p. 527, above eq. (14.14)
disp('Covariance matrix has been computed.')
