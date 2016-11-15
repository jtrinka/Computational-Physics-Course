%% Orbiting Masses
clear; clc;

d=1e6; %Distance between two masses
r=d/2; %half way to masses
m=1e15; %mass of two masses
G=6.67e-11; %Newtons universal law of G
v = sqrt((G*m*r)/(d^2)); %speed found using Gm^2/d^2=v^2*m/r
T = (2*pi*r)/v; %period found using circumference of circle 2pi*r/v
omega= 2*pi/T; %angular frequency
tmin=0;
tmax=3600*24*360;
n=100;
t = linspace(tmin,tmax,n);
%analytic solutions
xm1 = @(t) r*cos(omega*t); %solution is cos for x since we are starting at x=r. Can't include sine because we are not starting at any point up or down
ym1 = @(t) r*sin(omega*t); %solution is sine since are are starting at y=0. Can't include cos.
xm2 = @(t) -r*cos(omega*t);
ym2 = @(t) -r*sin(omega*t);

subplot(1,2,1)
plot(xm1(t),ym1(t))
hold on
plot(xm2(t),ym2(t),'ro')
xlabel('Horizontal Distance (m)')
ylabel('Vertical Distance (m)')
title('Analytical Solutions')

%Initial velocity for mass 1 is v since initially, we want mass 1 going
%straight up. Initial velocity for mass is -v since initially, we want mass
%2 going straight down.
[t,y] = ode45('TwoMasses',[tmin,tmax],[r,0,0,v,-r,0,0,-v],odeset('RelTol',10^-6));

subplot(1,2,2)
plot(y(:,1),y(:,2))
hold on
plot(y(:,5),y(:,6),'ro')

xlabel('Horizontal Distance (m)')
ylabel('Vertical Distance (m)')
title('Numerical Solutions')