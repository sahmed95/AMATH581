function [y_estH] = H(y0, dt, interval)
    syms y(t)
    syms f(t,y)

    first = interval(1); 
    last = interval(2);
    y_estH =[y0]; %list of estimate values using Heun's method
    for t = first:dt:last-dt
        f = -3*y0*sin(t); %function
        fun = double(subs(f,t));
        yz = y0 + dt*fun;
        syms g(x);
        g = -3*yz*sin(x);
        afun = double(subs(g,t+dt));
        y1 = y0 + (dt/2)*(fun+afun);
        y0 = y1;
        y_estH = [y_estH,y0];
 
    end
y_estH = y_estH';    
end
