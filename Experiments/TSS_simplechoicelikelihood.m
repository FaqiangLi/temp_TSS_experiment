function f = TSS_simplechoicelikelihood(I,storesize,distance,beta)


p = TSS_simplechoiceprob(storesize,distance,beta);

f = I.*log(p);
f = sum(sum(f));
if isnan(f)|isinf(f)
    f=-1e+5
end
