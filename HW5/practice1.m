%Problem 1

N = 20; 
[D,x] = cheb(N); 
DD = D*D;
D_2 = DD(2:N, 2:N);
D_1 = D(2:N, 2:N); 
I = eye(N-1);

for i = 2:N
    I(i-1,i-1) = exp(x(i));
end

A = D_2 + 4*D_1 + I;
rhs = sin(8*x(2:N)); 
u = A\rhs;

u = [0;u;0];

A1 = u(11); 


save A1.dat A1 -ascii

%Problem 2
% p14.m - solve nonlinear BVP u_xx = exp(u), u(-1)=u(1)=0
%         (compare p13.m)

  N = 16;
  [D,x] = cheb(N); D2 = D^2; D2 = D2(2:N,2:N);
  u = zeros(N-1,1);
  change = 1; it = 0;
  C =[change];
  while change > 10^(-15)% fixed-point iteration
    J = D2 - diag(exp(u)); %Jacobian
    F = D2*u-exp(u);
    unew = u - J^(-1)*F;
    %res = D2*unew - exp(unew);
    %obj = norm(res); 
    %grad = J'*res;
    change = norm(unew-u,inf);
    C= [C;change];
    u = unew; 
    it = it+1;
  end
  C
  
   u = [0;u;0];
  clf, subplot('position',[.1 .4 .8 .5])
  plot(x,u,'.','markersize',16)
  xx = -1:.01:1;
  uu = polyval(polyfit(x,u,N),xx);
  line(xx,uu), grid on
  title(sprintf('no. steps = %d      u(0) =%18.14f',it,u(N/2+1)))
  
  
  A2 = C(2);
  A3 = C(3);

  
  save A2.dat A2 -ascii
  save A3.dat A3 -ascii
  
  %Problem 3
  
  %Solving ut = uxx + e^u
  
  N = 20; 
  [D,x] = cheb(N); 
  DD = D*D; 
  D_2 = DD(2:N, 2:N);
  I = eye(N-1); 
  D_xx = kron(I, D_2);
  dt = 0.0001; 
  tspan = 0:dt:3.55; 
  %[t,y] = ode23s(D_xx-diag(exp(
  u0 = zeros(N-1,1);
  [t,u] = ode23s(@(t,u) D_2*u + exp(u), tspan, u0);
  m = length(tspan);
  z = zeros(m,1);
  u = [z, u, z];
  [I,J] = find(tspan==3.5);
  A4 = u(J,11);
  save A4.dat A4 -ascii
  t5 = 35443;
  A5 = tspan(t5);
  save A5.dat A5 -ascii
  
  
 
 