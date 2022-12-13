function [rhs] = shoot2(xspan, x, dummy, n0, beta)
%shoot2 builds RHS
rhs = [y(2); (x^2-beta)*y(1)]; 
end

