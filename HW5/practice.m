% p14.m - solve nonlinear BVP u_xx = exp(u), u(-1)=u(1)=0
%         (compare p13.m)

  N = 16;
  [D,x] = cheb(N); D2 = D^2; D2 = D2(2:N,2:N);
  u = zeros(N-1,1);
  %u = D2\exp(u);
  change = 1; it = 0;
  C =[change];
  while change > 10^(-14)          % fixed-point iteration
    J = D2 - diag(exp(u)); %Jacobian
    F = D2*u-exp(u);
    unew = u - J^(-1)*F;
    res = D2*unew - exp(unew);
    obj = norm(res); 
    change = obj;
    C= [C;change];
    u = unew; 
    it = it+1;
  end
  A2 = C(2);
  A3 = C(3);
  
  save A2.dat A2 -ascii
  save A3.dat A3 -ascii
   u = [0;u;0];
  clf, subplot('position',[.1 .4 .8 .5])
  plot(x,u,'.','markersize',16)
  xx = -1:.01:1;
  uu = polyval(polyfit(x,u,N),xx);
  line(xx,uu), grid on
  title(sprintf('no. steps = %d      u(0) =%18.14f',it,u(N/2+1)))