y0 = 1;  lam = -15; dt = 0.2;  T = 5; N = T/dt; 

ts = linspace(0, T, N); yTrue = y0*exp(lam*ts); z = lam*dt;

yFE = 0*ts; yFE(1)= y0; 
for i = 2:length(ts)   
    yFE(i) = yFE(i-1)*(1+z);  
end

figure()
plot(ts, yTrue, 'Linewidth', 2);
hold on;
plot(ts, yFE, 'Linewidth', 2);
legend('true solution', 'FE estimate')