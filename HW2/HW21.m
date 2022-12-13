%Shooting method 
clear all;
close all;
L = 4; %domain value
maxiter = 1000; %max iterations
xp = -L:0.1:L; %span of computational domain
tolerance = 10^(-4); %tolerance
A = 1;    
n0 = 1; %K
beta_initial = 0;
num_modes = 5;
eigenvalues=zeros(5,1);
for i = 1:num_modes
    beta = beta_initial;
    step = 1; 
    error = 1 ;
    iter = 0;
    x0 = [A A*sqrt(n0*L^2-beta)];
    while abs(error) > tolerance && iter < maxiter
        iter = iter + 1;
        [t,y] = ode45(@(t,y) rhs(t,y,beta), xp,x0);
        error = y(end,1)*sqrt(n0*L^2-beta)+y(end,2);
        eig = y(:,1);
        e1 = eig./sqrt(trapz(xp,eig.^2));
        if abs(error) < tolerance
            break
        elseif (-1)^(i+1)*error > 0
            beta = beta + step;
            step = step/2;
        else
            beta = beta - step;
            
            
     
        end
        
    end
    eigenvalues(i,:) = beta;
    eigenfunctions(:,i) = e1;
    beta_initial = beta +0.1;
   
   

end 

eigenvalues = abs(eigenvalues);
eigen1 = abs(eigenfunctions(:,1));
eigen2 = abs(eigenfunctions(:,2));
eigen3 = abs(eigenfunctions(:,3));
eigen4 = abs(eigenfunctions(:,4));
eigen5 = abs(eigenfunctions(:,5));

save A1.dat eigen1 -ascii
save A2.dat eigen2 -ascii
save A3.dat eigen3 -ascii
save A4.dat eigen4 -ascii
save A5.dat eigen5 -ascii
save A6.dat eigenvalues -ascii

