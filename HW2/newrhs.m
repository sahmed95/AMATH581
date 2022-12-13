function [y] = rhs(x,y,xp, x0,beta)

y1 = y(2);
y2 = (x^2-beta)*y(1);

[y] = [y1; y2];
end