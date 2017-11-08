function x = discreternd(p,n)
if nargin == 1
    n = 1;
end
k = length(p);
p = reshape(p,k,1);
x = sum(repmat(rand(1,n),k,1)> repmat(cumsum(p)/sum(p),1,n),1)+1;