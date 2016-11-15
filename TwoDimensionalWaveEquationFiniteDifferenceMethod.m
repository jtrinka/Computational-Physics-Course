%% Wave Equation Simulation
clear; clc;
%Analytic solution
syms x y t c a
z = cos(x*pi)*cos(y*(pi/3))*cos(c*t); %My guess at the analytic solution
zt = diff(z,t);
ztt = diff(zt,t);
zx = diff(z,x);
zxx = diff(zx,x);
zy = diff(z,y);
zyy = diff(zy,y);
lapz = zxx+zyy;
%Solve for c. 
cparam=solve(ztt==25*lapz,c);
%If we plug first entry into z, we get 0. Choosing the second or third entry yields same solution at center of the pool.
solz=subs(z,c,cparam(3,1));
%Solving for the solution when x and y are fixed to be at the center of the pool.
centerz= subs(solz,[x,y],[2,3]);
%Solving the period of oscillation by using T=2pi/omega where omega is the angular frequency and T is the period.
Tanalytic = 6*pi/(5*pi*sqrt(10));  

%Numerical Solution
%Define the domain
xmin = 0; xmax = 4; numx = 40;
ymin = 0; ymax = 6; numy = 40;
x=linspace(xmin,xmax,numx+2);
dx = x(2)-x(1);
y=linspace(ymin,ymax,numx+2);
dy = y(2)-y(1);
[x,y]=meshgrid(x,y);
%Define the initial conditions
IC = @(x,y) cos(pi.*x).*cos((pi/3).*y);
%Speed of propagation
m=@(x,y)5;% %5+.5.*x+.5.*y;
%Set up time interval, time vector, and dt.
ntimeit = 1112;
t = zeros(ntimeit,1);
t(1)=0;
dt = .0005;
%Pre-allocate space.
Z=zeros(ntimeit,length(x(:,1)),length(y(:,1)));
record = zeros(length(t),1);
%Initial Conditions
Z(1,:,:)= IC(x,y);
Z(2,:,:)=IC(x,y);
%Solution
for n=2:ntimeit
    for i=2:length(x(:,1))-1
        for j=2:length(y(:,1))-1
            
            
   %Interior Point Solution
   
   Z(n+1,i,j) = -Z(n-1,i,j)+2*Z(n,i,j)+((m(i,j)^2*((Z(n,i+1,j)-2*Z(n,i,j)+Z(n,i-1,j))/dx^2)+m(i,j)^2*((Z(n,i,j+1)-2*Z(n,i,j)+Z(n,i,j-1))/dy^2)))*dt^2;
       
   
   %left boundary solution
   Z(n+1,1,:)=Z(n+1,2,:);
   
   %right boundary solution
   Z(n+1,42,:) = Z(n+1,41,:);
   
   %bottom boundary solution
   Z(n+1,:,1)=Z(n+1,:,2);
   
   %top boundary solution
   
   Z(n+1,:,42) = Z(n+1,:,41);

            
        end
    end
    %Increment time vector
    t(n+1)=t(n)+dt;
    
    %Determining the period of oscillation for numerical solution 
    if Z(n+1,22,22)<=-1
        
        record(n+1)=n;
        
    end

   
   %Plot the numerical solution at each time step
    surf(reshape(Z(n+1,:,:),[42,42]));
    
            
            axis([0,length(x(:,1)),0,length(y(:,1)),-1,1])
            xlabel('Horizontal Distance (m)')
            ylabel('Vertical Distance (m)')
            zlabel('Depth (m)')
            title('Numerical Wave Equation')
            
 
    drawnow;
end
%Numerical period of oscillation using Time=n*dt.
Tnumerical = 942*dt;
%Absolute value difference between the analytic and numerical periods of
%oscillation
Tdiff = abs(Tanalytic-Tnumerical)
%The depth of the pool after one second at the center of the pool.
Depth=Z(end,22,22)