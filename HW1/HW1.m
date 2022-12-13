%Problem 1 a (Euler) 


Error = []; 
dtlist = [];%list of timestep
for k = [2:8]
    dt = [2^(-k)];
    dtlist =[dtlist,dt];
    y_est = E(pi/sqrt(2),dt,[0,5]); %list of y estimates for particular timestep
    y_true = L(pi/sqrt(2),dt,[0,5]); %list of true y for particular timestep
    e = mean(abs(y_est - y_true)); %mean error for a particular timestep
    Error = [Error, e]; %list of errors
end

last = E(pi/sqrt(2), dtlist(7),[0,5]); %list of y estimates for 2^-8



plot(log(dtlist), log(Error)), xlabel('log(dt)'), ylabel('log(Error)'), title('Euler method');
n = polyfit(log(dtlist), log(Error), 1);
slope = n(1);

save A1.dat last -ascii
save A2.dat Error -ascii
save A3.dat slope -ascii

%Problem 1 b (Heun)
Error_H = []; 
dtlist_H = [];%list of timestep
for k = [2:8]
    dt_H = [2^(-k)];
    dtlist_H =[dtlist_H,dt_H];
    y_estH = H(pi/sqrt(2),dt_H,[0,5]); %list of y estimates for particular timestep
    y_true = L(pi/sqrt(2),dt_H,[0,5]); %list of true y for particular timestep
    e = mean(abs(y_estH - y_true)); %mean error for a particular timestep
    Error_H = [Error_H, e]; %list of errors
end

last_H = H(pi/sqrt(2), dtlist_H(7),[0,5]); %list of y estimates for 2^-8


plot(log(dtlist_H), log(Error_H)), xlabel('log(dt)'), ylabel('log(Error)'), title('Euler method');
n = polyfit(log(dtlist_H), log(Error_H), 1);
slope_H = n(1);

save A4.dat last_H -ascii
save A5.dat Error_H -ascii
save A6.dat slope_H -ascii

%Problem 2 a) Van der Pol oscillator (we reduce the order of the
%differential equation and solve a system instead)

tspan = [0: 0.5: 32]; %timespan

epsilon = [0.1, 1, 20]; %list of epsilon values
Y=[]; %initializing list of y(t)
for i = 1:length(epsilon)
    [t1,y1] = ode45(@(t,y) odefun(t,y,epsilon(i)),tspan, [sqrt(3) 1]);
    Y =[Y,y1(:,1)];
end

save A7.dat Y -ascii

%Problem 2 b)

tspan = [0,32];
TOL = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9,1e-10]; %list of tolerances

%calculation for ode45
timestep1 = []; %list of timesteps
for i =1:length(TOL)
    options = odeset('AbsTol',TOL(i),'RelTol',TOL(i));
    [T,Y] = ode45(@(t,y) odefun(t,y,1), tspan,[2,pi^2],options);
    dt = mean(diff(T)); %finding the timestep for a particular tolerance
    timestep1 = [timestep1, dt];
end


plot(log(timestep1), log(TOL)), xlabel('log(dt)'), ylabel('log(TOL)'), title('Tolerance');
slope1 = polyfit(log(timestep1), log(TOL), 1);
slope1 = slope1(1);

%calculation for ode23
timestep2 = []; 
for i =1:length(TOL)
    options = odeset('AbsTol',TOL(i),'RelTol',TOL(i));
    [T,Y] = ode23(@(t,y) odefun(t,y,1), tspan,[2,pi^2],options);
    dt = mean(diff(T));
    timestep2 = [timestep2, dt];
end


plot(log(timestep2), log(TOL)), xlabel('log(dt)'), ylabel('log(TOL)'), title('Tolerance');
slope2 = polyfit(log(timestep2), log(TOL), 1);
slope2 = slope2(1);


%caluclation for ode113
timestep3 = [];
for i =1:length(TOL)
    options = odeset('AbsTol',TOL(i),'RelTol',TOL(i));
    [T,Y] = ode113(@(t,y) odefun(t,y,1), tspan,[2,pi^2],options);
    dt = mean(diff(T));
    timestep3 = [timestep3, dt];
end


plot(log(timestep3), log(TOL)), xlabel('log(dt)'), ylabel('log(TOL)'), title('Tolerance');
slope3 = polyfit(log(timestep3), log(TOL), 1);
slope3 = slope3(1);

save A8.dat slope1 -ascii
save A9.dat slope2 -ascii
save A10.dat slope3 -ascii


