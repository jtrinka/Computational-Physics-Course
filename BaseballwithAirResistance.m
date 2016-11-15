%% A Baseball with Air Resistance Euler Cromer Adaptive Time Step
clear;clc;
m=.145; %kg
g=9.81; %m/s^2
rho = 1.225; %Density of air is kg/m^3
Cd = .3; % Drag coefficient
A = pi*(.038)^2; %Cross section area of a baseball in meters ^2
ys(1)=1.8796000000000002; %meters
xs(1)=0; %meters
vs(1)=44.704; %m/s
yd(1)=1.8796000000000002; %meters
xd(1)=0; %meters
yt(1)=1.8796000000000002; %meters
xt(1)=0; %meters
vd(1)=44.704; %m/s
theta = pi/6; %rad 
vys(1)=vs(1)*sin(theta); %meters per second
vxs(1)=vs(1)*cos(theta); %meters per second
vyd(1)=vd(1)*sin(theta); %meters per second
vxd(1)=vd(1)*cos(theta); %meters per second
% tmin=0; %seconds
% tmax=4; %seconds
% n=100000; %interior points
% t = linspace(tmin,tmax,n)'; %time
t(1)=0;
dt =4e-5;%t(2)-t(1); %delta t
tolg=10^-2;
toll=10^-8;
i=1;
n=60000;
%for i=1:1:n
 while yt(i)>0   
    %Euler-Cromer first single step
    
    %We want to compare the double step to the second single step so don't
    %save the first step.
    ayfs = -g-(vys(i)*.5*rho*vs(i)*Cd*A)/m;
    axfs = -(vxs(i)*.5*rho*vs(i)*Cd*A)/m;
    vyfs = vys(i)+ayfs*dt;
    vxfs = vxs(i)+axfs*dt;
    yfs = ys(i)+vyfs*dt;
    xfs = xs(i)+vxfs*dt;
    vfs = sqrt(vyfs^2+vxfs^2);
    
    
    %Euler-Cromer second single step
    
    
    ays(i)= -g-(vyfs*.5*rho*vfs*Cd*A)/m;
    axs(i) = -(vxfs*.5*rho*vfs*Cd*A)/m;
    vys(i+1) = vyfs+ays(i)*dt;
    vxs(i+1) = vxfs+axs(i)*dt;
    ys(i+1) = yfs+vys(i+1)*dt;
    xs(i+1) = xfs+vxs(i+1)*dt;
    vs(i+1) = sqrt(vys(i+1)^2+vxs(i+1)^2);
    
    %Euler-Cromer Double Step
    
    
    
    ayd(i) = -g-(vyd(i)*.5*rho*vd(i)*Cd*A)/m;
    axd(i) = -(vxd(i)*.5*rho*vd(i)*Cd*A)/m;
    vyd(i+1) = vyd(i)+ayd(i)*dt*2;
    vxd(i+1) = vxd(i)+axd(i)*dt*2;
    yd(i+1) = yd(i)+vyd(i+1)*dt*2;
    xd(i+1) = xd(i)+vxd(i+1)*dt*2;
    vd(i+1) = sqrt(vyd(i+1)^2+vxd(i+1)^2);
    
    
    diff = sqrt((xd(i+1)-xs(i+1))^2+(yd(i+1)-ys(i+1))^2);
    
    
    t(i+1)=t(i)+2*dt;
    if diff>tolg
        
        yt(i+1)=ys(i+1);
        xt(i+1)=xs(i+1);
        dt=dt/2;
        
    elseif diff>toll
        
        yt(i+1)=ys(i+1);
        xt(i+1)=xs(i+1);
        
    else
        
        yt(i+1)=yd(i+1);
        xt(i+1)=xd(i+1);
        
        dt=dt*1.5;
        
        
        
        
    end
i=i+1;

end


xs=xs';
ys=ys';
xd=xd';
yd=yd';
xt=xt';
yt=yt';


% plot(xs(:,1),ys(:,1))
% hold on
% plot(xd(:,1),yd(:,1))


figure
plot(xt(:,1),yt(:,1))
xlabel('Horizontal Distance (m)')
ylabel('Vertical Distance (m)')
title('Adaptive Time Step Euler-Cromer Method with Air Resistance')