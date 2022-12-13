syms x;
f = x*sin(3*x)-exp(x);
g = diff(f);
x0 = -1.6;
X = [x0];
iterations_N = 1;
value = double(subs(f,x0));

while abs(value)>=10^-6
    x1 = x0 - value/double(subs(g,x0));
    diff = x1 - x0;
    x0 = x1;
    X = [X,x0];
    iterations_N = iterations_N + 1;
    value = double(subs(f,x0));
end
X;

save A1.dat  X  -ascii


xl = -0.7;
xr = -0.4;
xmid = (xl+xr)/2;
g = double(subs(f,xmid));
iterations_b = 1;
X_bisection = [xmid];
while abs(g)>=10^(-6)
   iterations_b = iterations_b + 1;
   l = double(subs(f,xl));
   r = double(subs(f,xr));

   if (l > 0 & g>0)| (l<0 & g<0)
      xl = xmid;
   else
      xr = xmid;
   end
   xmid = (xl+xr)/2;
   X_bisection = [X_bisection,xmid];
   g = double(subs(f,xmid));
end 
iterations = [iterations_N,iterations_b];
X_bisection;


save A2.dat X_bisection -ascii
save A3.dat iterations  -ascii


A = [1,2;-1,1]; B = [2,0;0,2]; C = [2,0,-3;0,0,-1]; D = [1,2;2,3;-1,0]; x=[1;0];y=[0;1];z=[1;2;-1];

M = A+B;

save A4.dat M -ascii

M = (3*x) - (4*y);
save A5.dat M -ascii
M = A*x;
save A6.dat M -ascii
M = B*(x-y);
save A7.dat M -ascii
M = (D*x);
save A8.dat M -ascii
M = (D*y) + z;
save A9.dat M -ascii
M = A*B;
save A10.dat M -ascii
M = B*C; 
save A11.dat M -ascii
M = C*D;
save A12.dat M -ascii

