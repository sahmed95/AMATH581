function [error,t, y] = myObj(beta, params)
%myObj objective function for shooting 

xp = params.xp;
x0 = params.x0;
n0 = params.n0;
L = params.L;

[t,y] = ode45('newrhs', xp, x0, [], n0, beta);

error = abs(y(end,1)*sqrt(n0*L^2-beta)+y(end,2)); 

end

