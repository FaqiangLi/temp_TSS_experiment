function [f,g] = yearindices(yer)

n = length(yer);
y_l = zeros(n,1);
for j=1:n-1
    y_l(j) = yer(j)~=yer(j+1);
end
y_l(n) = true;
y_i = cumsum(y_l)-y_l+1;
f = y_l;
g = y_i;