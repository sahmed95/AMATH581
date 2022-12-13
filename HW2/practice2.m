clear all; close all; clc;
n0 = 1;
A = 1;
xp = -4:0.1:4;
betaIn = [0,2,4,6,8];
num_modes = 5;

for i=1:num_modes
    betain = betaIn(i)
    x0 = [A A*sqrt(n0*16-betain)];
    fobj = (@betain)ode45(rhs(t,y,betain),xp,x0);
    beta = fminunc(fobj, betaIn(i))
    beta
end