function output = TwoMasses(t,y)
%Mass 1
xpE1 = y(1);
ypE1 = y(2);
vxE1 = y(3);
vyE1 = y(4);

%Mass 2
xpE2 = y(5);
ypE2 = y(6);
vxE2 = y(7);
vyE2 = y(8);

G = 6.67e-11; %Newton's gravitational constant
Me= 1e15; %kg

xpd = xpE1-xpE2; %magnitude of distance separating two masses in x direction
ypd = ypE1-ypE2; %magnitude of distance separating two masses in y direction



%The differential equations
output = [ vxE1, vyE1, (-G*Me*xpd)/(sqrt(xpd.^2+ypd.^2)^3), (-G*Me*ypd)/(sqrt(xpd.^2+ypd.^2)^3),vxE2, vyE2, (-G*Me*-xpd)/(sqrt(xpd.^2+ypd.^2)^3), (-G*Me*-ypd)/(sqrt(xpd.^2+ypd.^2)^3)]'; %differential equations right here. Answers returned in same columns.
