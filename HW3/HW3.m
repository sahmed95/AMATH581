%1a) Dx without periodic condition
n = 8;
span = [-1,1];
l = linspace(-1,1,n+1);
l = l(1:n);
h = mean(diff(l));

%Dx
e = ones(n,1);
A1 = spdiags([-e e], [-1 1], n,n)/(2*h);
A1 = full(A1);
save A1.dat A1 -ascii

%Dx periodic 
A2 = spdiags([-e e], [-1 1], n,n);
A2(1,end) = -1;
A2(end,1) = 1;
A2 = A2/(2*h);
A2 = full(A2);

save A2.dat A2 -ascii
%Dxx 
A3 = (1/h^2)*spdiags([e -2*e e], -1:1, n, n); 
A3 = full(A3);
save A3.dat A3 -ascii

%Dxx periodic
A4 = spdiags([e -2*e e], -1:1, n, n);
A4(1,end) = 1;
A4(end,1)= 1;
A4 = A4*(1/h^2);
A4 = full(A4); 
save A4.dat A4 -ascii

%Dx 
P = spdiags([-e e], [-1 1], n,n);
I = speye(n);
A5 = kron(I,P)/(2*h);
A5 = full(A5); 
save A5.dat A5 -ascii
%Dx periodic
P = spdiags([-e e], [-1 1], n,n);
P(1,end) = -1;
P(end, 1) = 1;
A6 = kron(I,P)/(2*h);
A6 = full(A6);
save A6.dat A6 -ascii
%Dy 
P = spdiags([-e e], [-1 1], n,n);
A7 = kron(P,I)/(2*h);
A7 = full(A7);
save A7.dat A7 -ascii
%Dy periodic
P = spdiags([-e e], [-1 1], n,n);
P(1,end) = -1; 
P(end,1) =1;
A8 =  kron(P,I)/(2*h);
A8 = full(A8); 

save A8.dat A8 -ascii

%Dxx 
P = spdiags([e -2*e e], -1:1, n, n); 
A9 = kron(I,P)/(h^2);
A9 = full(A9);

save A9.dat A9 -ascii
%Dxx periodic
P = spdiags([e -2*e e], -1:1, n, n); 
P(1,end) = 1;
P(end,1) = 1;
A10 = kron(I,P)/(h^2);
A10 = full(A10);
save A10.dat A10 -ascii
%Dyy 
P = spdiags([e -2*e e], -1:1, n, n); 
A11 = kron(P,I)/(h^2);
A11 =full(A11);
size(A11)
save A11.dat A11 -ascii

%Dyy periodic
P = spdiags([e -2*e e], -1:1, n, n); 
P(1,end) = 1;
P(end,1) = 1;
A12 = kron(P,I)/(h^2);
A12 = full(A12);

save A12.dat A12 -ascii
%Laplacian
T = spdiags([e -4*e e], -1:1, n,n);
S = spdiags([e e], [-1 1], n,n); 
A13 = (kron(I, T) + kron(S,I))/(h^2);
A13 = full(A13);
save A13.dat A13 -ascii 

%Laplacian periodic 
T = spdiags([e -4*e e], -1:1, n,n);
S = spdiags([e e], [-1 1], n,n); 
T(1,end)=1;
T(end,1) =1; 
R = sparse(n,n);
R(n,1) = 1;
R(1,n) = 1; 
A14 = (kron(I,T)+kron(S,I)+kron(R,I))/(h^2);
A14 = full(A14);
save A14.dat A14 -ascii

%2a)
n = 10;
tolerance = 1*10^-8;
b = ones(n, 1);
l = 1.16*b;
k = 0.16*b;
A11 = spdiags([-l b k], -1:1, n, n);
M = diag(A11);
M = spdiags([M], 0, n,n);
x = zeros(n,1); 

r = b-A11*x;
error = norm(r);
iter = 0;
i = [iter];
E = [norm(r)];
z = M\r;
while error>tolerance
    iter = iter+1;
    i = [i;iter];
    x = x+z; 
    r = b-A11*x;
    error = norm(r);
    E =[E;error];
    z = M\r;
end
E=E';
save A15.dat E -ascii


n =30;

b = ones(n, 1);
l = 1.16*b;
k = 0.16*b;
A11 = spdiags([-l b k], -1:1, n, n);
M = diag(A11);
M = spdiags([M], 0, n,n);
x = zeros(n,1); 
iter = 0;
i = [iter];
r = b-A11*x;
error = norm(r);
E =[error];

z = M\r;
while error>tolerance
    iter = iter +1; 
    i = [i; iter];
    x = x+z; 
    r = b-A11*x;
    error = norm(r);
    E = [E; error];
    z = M\r;
end
E=E';

save A16.dat E -ascii


n = 10;
tolerance = 1*10^-8;
b = ones(n, 1);
l = 1.16*b;
k = 0.16*b;
A11 = spdiags([-l b k], -1:1, n, n);
M = tril(A11);
x = zeros(n,1); 
iter = 0;
i =[iter];
r = b-A11*x;
error = norm(r);
E=[error];


z = M\r;
while error>tolerance
    iter = iter +1;
    i = [i; iter];
    x = x+z; 
    r = b-A11*x;
    error = norm(r);
    E =[E; error];
    z = M\r;
end
E=E';
save A17.dat E -ascii

n =30;

b = ones(n, 1);
l = 1.16*b;
k = 0.16*b;
A11 = spdiags([-l b k], -1:1, n, n);
M = tril(A11);
x = zeros(n,1); 
iter = 0; 
i = [iter];
r = b-A11*x;
error = norm(r);
E = [error];


z = M\r;
while error>tolerance
    iter = iter +1;
    i = [i; iter];
    x = x+z; 
    r = b-A11*x;
    error = norm(r);
    E = [E; error];
    z = M\r;
end
E =E';
save A18.dat E -ascii

