function [y_true] = L(y0,dt,interval)
    first = interval(1);
    last = interval(2);
    syms y(s)
    y = (pi*exp(3*cos(s)-3))/sqrt(2);
    y_true = []; %list of true y values
    for s = first:dt:last
        value = double(subs(y,s));
        y_true = [y_true, value];
        
 
    end
y_true = y_true';    
end
    