function [f,g] = TSS_gmm(theta,inp)

new_moments = inp.new_moments;
W = inp.W;
N = inp.N;
instr_ix = inp.instr_ix;            % logical index to drop unused cross-cat moment conditions

[G1,G2,G3,G4] = TSS_moments(theta,inp);
g = [G1; G2; G3]/N;
g = g(instr_ix);
if new_moments
    g = [g; G4/N];
end
if all(isfinite(g(:))) 
    f = g'*W*g;
else
    f = 1e10;
end