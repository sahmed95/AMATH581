clear all; close all; clc;
col = ['r', 'b', 'g', 'c', 'm', 'k'];


n0 = 1; % n0 as in book
A = 1;    % initial slope at x = -1 
L =4; 
xp = -L:0.1:L; % span of computational domain ; 
betaIn = [0, 1.8, 4.9, 6, 8.5]; 
params.L = L;
params.xp = xp; 
params.n0 = n0; 
modenum = 5; 
betas = zeros(modenum,1);
iterations = zeros(modenum,1);

for modes = 1:modenum
    x0 = [A A*sqrt(n0*L^2-betaIn(modes))];
    params.x0 = x0;
    fobj = @(beta)myObj(beta, params);
    x0 = [A A*sqrt(n0*L^2-betaIn(modes))];
    params.x0 = x0;
    [beta,fval,exitflag,output] = fminunc(fobj, betaIn(modes)); % run our optimizatiom
    betas(modes) = beta;
    iterations(modes) = getfield(output, 'iterations');
    [error,t,y] = fobj(beta); % get the y
    
end
save A13.dat betas -ascii
save A14.dat iterations -ascii
