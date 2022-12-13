function [y_est] = E(y0, dt, interval)
    syms y(t)
    syms f(t,y)

    first = interval(1); 
    last = interval(2);
    y_est =[y0]; %list of estimate values using Euler
    for t = first:dt:last-dt
        f = -3*y0*sin(t); %function
        func = double(subs(f,t));
        y1 = y0 + dt*func;
        y0 = y1;
        y_est = [y_est,y0];
 
    end
y_est = y_est';    
end
