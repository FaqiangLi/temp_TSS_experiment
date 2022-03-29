function B2 = TSS_notin(A,B1,N2)

% Given a set A (with elements listed in an N0x1 array), and a subset B1 of
% A (with elements listed in an N1x1 array, N1<N0), draw a new subset B2 of A
% of size N2<=N0-N1 such that the intersection of B1 and B2 is empty. 

ix_notB1 = ~ismember(A,B1);     % ix_notB1 has the same size as first argument in ismember.
AnotB1 = A(ix_notB1);
N0N1 = length(AnotB1);
ix_B2 = randperm(N0N1,N2);
B2 = AnotB1(ix_B2);
