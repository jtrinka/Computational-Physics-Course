%% Schrodinger Equation Analysis
clear; clc;
xmin=0; %meters
xmax=1e-14; %meters
%initial conditions
psi0=1;
dpsidx0=0; 
E01=-1e-10; %Joules
E02=-1e-11; %Joules

%With psi0=1 and dpsidx0 = 0, there is an eigen-energy level between
%E01=-1e-14 joules and E02=-1e-11 joules. There is another between
%E01=-1e-11 joules and E02=-1e-10 joules.
[t1,y1]=ode45('Schrodinger1SpatialD',[xmin,xmax],[psi0,dpsidx0,E01],odeset('RelTol',10^-6));
[t2,y2]=ode45('Schrodinger1SpatialD',[xmin,xmax],[psi0,dpsidx0,E02],odeset('RelTol',10^-6));

subplot(1,2,1)
plot(t1,y1(:,1))
xlabel('Distance (m)')
ylabel('Wave Function')
title(['Search for Eigen-Energy Levels E01 = ',num2str(E01)])
subplot(1,2,2)
plot(t2,y2(:,1))
xlabel('Distance (m)')
ylabel('Wave Function')
title(['Search for Eigen-Energy Levels E02 = ',num2str(E02)])

% To change potential well, simply change V in Schrodinger1SpatialD to
% V(x)-400exp(-x/.3fm). 400 is in MeV so change to joules with 400*1.60218e-13, .3 is in
% femtometers so changes to meters with .3e-15.
