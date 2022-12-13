%1) Solve using Forward Euler

%value of kappa
k = 0.1; 
%time and space step
dt = 0.05;
dx = 0.05;
xspan = 0:dx:1; 
tspan = 0:dt:1;
m = length(xspan);
n = length(tspan);
gamma = dt/(dx)^2;
%initial condition
init = @(x) sin(2*pi*x);

%creating discretized A for Forward Euler using periodic boundary condition
e = ones(m-1,1); 
A = spdiags([e -2*e e], -1:1, m-1,m-1);
A(1,end) = 1;
A(end,1) = 1;

U = zeros(m-1,n); 
for i = 1:m-1
    U(i,1) = init(xspan(i));
end


for i = 2:n
    U(:,i) =U(:,i-1)+ k*gamma*A*U(:,i-1);
    
end

U = [U;U(1,:)];
A1 = U(:,end);
for i = 1:n
    plot(U(:,i))
    hold on
    getframe
    continue
end
save A1.dat A1 -ascii

%b) Solve using Backward Euler (implicit)
%value of kappa
k = 0.1; 
%time and space step
dt = 0.05;
dx = 0.1;
xspan = 0:dx:1; 
tspan = 0:dt:1;
m = length(xspan);
n = length(tspan);
gamma = dt/dx^2;

%initial condition
init = @(x) sin(2*pi*x);

%creating discretized A for Forward Euler using periodic boundary condition
e = ones(m-1,1); 
A = spdiags([e -2*e e], -1:1, m-1,m-1);
A(1,end) = 1;
A(end,1) = 1;

U = zeros(m-1,n); 
for i = 1:m-1
    U(i,1) = init(xspan(i));
end

I = speye(m-1); 
A = I - (k*gamma*A);
inverse = inv(A);

for i = 2:n
    U(:,i) = inverse*U(:,i-1);   
    
end
U = [U;U(1,:)];
A2 = U(:,end);
save A2.dat A2 -ascii
%changing time and space step 
dt = 0.01;
dx = 0.01;
xspan = 0:dx:1; 
tspan = 0:dt:1;
m = length(xspan);
n = length(tspan);
gamma = dt/dx^2;

%initial condition
init = @(x) sin(2*pi*x);

%creating discretized A for Forward Euler using periodic boundary condition
e = ones(m-1,1); 
A = spdiags([e -2*e e], -1:1, m-1,m-1);
A(1,end) = 1;
A(end,1) = 1;
U = zeros(m-1,n); 
for i = 1:m-1
    U(i,1) = init(xspan(i));
end

I = speye(m-1); 
A = I - (k*gamma*A);
inverse = inv(A);

for i = 2:n
    U(:,i) = inverse*U(:,i-1);   
    
end
U = [U;U(1,:)];
A3 = U(:,end);
save A3.dat A3 -ascii

%Solve using Crank Nicolson method

dt = 0.05;
dx = 0.1;
xspan = 0:dx:1; 
tspan = 0:dt:1;
m = length(xspan);
n = length(tspan);
gamma = dt/(2*dx^2);

%initial condition
init = @(x) sin(2*pi*x);

%creating discretized A for Forward Euler using periodic boundary condition
e = ones(m-1,1); 
A = spdiags([e -2*e e], -1:1, m-1,m-1);

A(1,end) = 1;
A(end,1) = 1;

U = zeros(m-1,n); 
for i = 1:m-1
    U(i,1) = init(xspan(i));
end

I = speye(m-1);
B = I - (k*gamma*A);
inverse = inv(B);

for i = 2:n
    U(:,i) = inverse*(k*gamma*A*U(:,i-1)+U(:,i-1));
end
U = [U; U(1,:)];
A4 = U(:,end);
save A4.dat A4 -ascii

%changing time and space step
dt = 0.01;
dx = 0.01;
xspan = 0:dx:1; 
tspan = 0:dt:1;
m = length(xspan);
n = length(tspan);
gamma = dt/(2*dx^2);

%initial condition
init = @(x) sin(2*pi*x);

%creating discretized A for Forward Euler using periodic boundary condition
e = ones(m-1,1); 
A = spdiags([e -2*e e], -1:1, m-1,m-1);

A(1,end) = 1;
A(end,1) = 1;

U = zeros(m-1,n); 
for i = 1:m-1
    U(i,1) = init(xspan(i));
end

I = speye(m-1);
B = I - (k*gamma*A);
inverse = inv(B);

for i = 2:n
    U(:,i) = inverse*(k*gamma*A*U(:,i-1)+U(:,i-1));
end
U = [U;U(1,:)];
A5 = U(:,end);

save A5.dat A5 -ascii

