%% Electrong in a Magnetic Field
clear;clc;
q=1.602e-19; %Coulombs
Bz = 3; %Teslas
m = 9.1e-31; %kg
omega = sqrt(((q*Bz)/m)^2);
vx0=1; %m/s
vy0=0; %m/s
ax0=(-q*vy0*Bz)/m; %m/s^2
ay0=(q*vx0*Bz)/m; %m/s^2
x0 = 0; %initial position x
y0 = 0; %initial position y
%solving for constants
c1x=vx0; 
c2x=ax0/omega;
c3x=x0+c2x/omega;
c1y=vy0;
c2y=ay0/omega;
c3y=y0+c2y/omega;
tmin=0;
tmax=1e-10;
n=10000;
t=linspace(tmin,tmax,n)';
%analytic solution
solx=@(t) (c1x/omega)*sin(omega*t)-(c2x/omega)*cos(omega*t)+c3x;
soly=@(t) (c1y/omega)*sin(omega*t)-(c2y/omega)*cos(omega*t)+c3y;


plot(solx(t),soly(t))
hold on
[t,y]=ode45('Electron',[tmin,tmax],[x0,y0,0,vx0,vy0,0],odeset('RelTol',10^-6));
%figure
plot(y(:,1),y(:,2),'ro')
xlabel('Horizontal Distance (m)')
ylabel('Vertical Distance (m)')
title('Analytic and Numerical Solutions')

Energy = .5*m*(y(:,4).^2+y(:,5).^2);
figure
plot(t,Energy(:,1))
xlabel('Time (s)')
ylabel('Energy (J)')
title('Conservation of Energy')

%To check  trajectory if magnetic field is a linear function x while
%pointing in x direction, turn off hold on, turn on figure command, and
%change Bz to Bz*xp.
