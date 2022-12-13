function p = bary(x)
n = length(x);
x= x(:);
for k =1:n
    p(k) = 1/(prod(x(k)-x([1:k-1,k+1:n])));
end 

