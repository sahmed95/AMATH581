N = 16;
xx = -1.01:.005:1.01; clf
for i = 1:2
    if i==1
        s = "equispaced points";
        x = -1 + 2*(0:N)/N; 
    end
    if i==2
        s = "Chebyshev points"; 
        x = cos(pi*(0:N)/N);
    end
    
    subplot(2,2,i)
    u = 1./(1+16*x.^2);
    uu = 1./(1+16*xx.^2);
    a = bary(x);
    p = polyfit(x,u,N);% interpolation
    mean(diff(a-p))
    pp = polyval(p,xx); % evaluation of interpolant

end
