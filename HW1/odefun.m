function dydt = odefun(t,y, epsilon)
dydt = zeros(2,1);
u1 = y(1);
u2 = y(2);
dydt(1) = u2;
dydt(2) = -u1 -epsilon*(u1^2-1)*u2;
end



