%% Traffic Flow Simulation

clear;clc;
vmax=35;
rhomax=1/3;
xmin=-100;
xmax=100;
nptx=100;
x=linspace(xmin,xmax,nptx+2);
dx=x(2)-x(1);
tmax=2;
dt=.001;
nptt=tmax/dt;
a=vmax/rhomax;
rho=zeros(nptt,length(x));
rho(1,1:nptx/2) = rhomax;
rho(1,nptx+1:end) = 0;
ConservationofCars = zeros(nptt,1);
ConservationofCars(1,1) = sum(rho(1,:));
DiffRegion1 = zeros(nptt,1);
DiffRegion2 = zeros(nptt,1);
DiffRegion3 = zeros(nptt,1);
for n=1:nptt
    for i=2:length(x)-1
        

        
        
        
 rho(n+1,i) = rho(n,i)+dt*-(...
     vmax*(rho(n,i+1)-rho(n,i-1))/2/dx*(1-2*(rho(n,i)/rhomax)));
        
        
        
    end
    
   rho(n+1,1)=rho(n+1,2);
    rho(n+1,end)=rho(n+1,end-1);
    ConservationofCars(n+1,1)=sum(rho(n+1,:));
    
DiffRegion1(n+1) = abs(rhomax-rho(n+1,2));

DiffRegion2(n+1) = abs((.5*(1)*rhomax)-rho(n+1,nptx/2+1));

DiffRegion3(n+1) = abs(0-rho(n+1,nptx-1));    
     plot(x,rho(n+1,:))
     axis([xmin,xmax,0,rhomax])
     xlabel ('Distance (m)')
     ylabel('Density of Cars (\rho)')
     title(n)
     
     drawnow;
    
end

figure
plot(0:nptt,ConservationofCars)
xlabel('Time Steps')
ylabel('Total Density of Cars')
title('Conservation of Cars')

%check solution at three regions

figure
plot(0:nptt,DiffRegion1)
xlabel('Time Step')
ylabel('Difference between the Numerical and Analytic Solution')
title('Difference between Numerical and Analytic Solution at x=-96')
figure
plot(0:nptt,DiffRegion2)
xlabel('Time Step')
ylabel('Difference between the Numerical and Analytic Solution')
title('Difference between Numerical and Analytic Solution at x=0')
figure
plot(0:nptt,DiffRegion3)
xlabel('Time Step')
ylabel('Difference between the Numerical and Analytic Solution')
title('Difference between Numerical and Analytic Solution at x=96')
