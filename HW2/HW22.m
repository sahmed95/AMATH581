%Direct solve
clear all; close all; clc;
L = 4;
x0 = -L;
xf = L;
h = 0.1;
xspan = -L:0.1:L;
N = (xf-x0)/h; % number of steps 
N = N-1;
assert(N-round(N) == 0); % make sure number of steps is an integer.
e = ones(N,1);
A = (1/h^2)*spdiags([e -2*e e], -1:1, N, N); 
A = -A;
xp = xspan(:,2:80)';
xp = xp.^2;
B = spdiags( [xp], 0, N, N);
A = A +B;
A(1,1) = (1/h^2)*(2/3)+xspan(1)^2;

A(1,2) = -(1/h^2)*(2/3);
A(N,N) = (1/h^2)*(2/3)+xspan(81)^2;

A(N,N-1) = -(1/h^2)*(2/3);
A
[V,D] = eig(full(A));
[d,ind] = sort(diag(abs(D)));
Ds = D(ind,ind);
Vs = V(:,ind);
eigenvalues = [];
for i =1:5
    eigenvalues(i,:) = Ds(i,i);
end
eig1 = Vs(:,1);
eig2 = Vs(:,2);
eig3 = Vs(:,3);
eig4 = Vs(:,4);
eig5 = Vs(:,5);
K = [];
for i=1:5
    k = 1/(2+3*h*sqrt(L^2-eigenvalues(i)));
    K(i)=k;
end
eig1 = [(eig1(1)+eig1(3))*K(1); eig1];
eig2 = [(eig2(1)+eig2(3))*K(2); eig2];
eig3 = [(eig3(1)+eig3(3))*K(3); eig3];
eig4 = [(eig4(1)+eig4(3))*K(4); eig4];
eig5 = [(eig5(1)+eig5(3))*K(5); eig5];

eig1 = [eig1; (eig1(78)+eig1(79))*K(1)];
eig2 = [eig2; (eig2(78)+eig2(79))*K(2)];
eig3 = [eig3; (eig3(78)+eig3(79))*K(3)];
eig4 = [eig4; (eig4(78)+eig4(79))*K(4)];
eig5 = [eig5; (eig5(78)+eig5(79))*K(5)];
    
e1 = abs(eig1./sqrt(trapz(xspan,eig1.^2)));
e2 = abs(eig2./sqrt(trapz(xspan,eig2.^2)));
e3 = abs(eig3./sqrt(trapz(xspan,eig3.^2)));
e4 = abs(eig4./sqrt(trapz(xspan,eig4.^2)));
e5 = abs(eig5./sqrt(trapz(xspan,eig5.^2)));
save A7.dat e1 -ascii
save A8.dat e2 -ascii
save A9.dat e3 -ascii
save A10.dat e4 -ascii
save A11.dat e5 -ascii
save A12.dat eigenvalues -ascii

