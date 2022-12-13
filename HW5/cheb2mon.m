function A = cheb2mon(N)

%constructs A to switch basis 
syms x;
T = [1,x];

for i = 3:N+1
    t = 2*x*T(i-1)-T(i-2);
    T = [T,t];
    
end
B = [];
co = [];
for i = 1:N+1
    c = coeffs(T(i), 'All');
    c = fliplr(c);
    m = length(c);
    if m < N+1
        for i = m+1:N+1
            c(i) = 0;
        end
    end 
    c = c';
    B = [B, c];
    
end


A = zeros(N+1, N+1);
for i = 1: N+1
    for j = 1:N+1
        A(i,j) = B(i,j);
    end
end


