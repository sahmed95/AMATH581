a = [1; -2; 3; 2; 1; -1];
A = cheb2mon(5); 
c = A*a;
save A8.dat c -ascii

b = [1; -2; 3; 2; 1; -1]; 
B = inv(A); 
a_n = B*b;
save A9.dat a_n -ascii