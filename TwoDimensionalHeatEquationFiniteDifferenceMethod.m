%% Two-Dimensional Heat Equation
clear; clc;
%Analytic Solution
syms t
centerx = 1.5;
centery = 2.5;
k1 = pi/3;
k2 = pi/5;
alpha = 5.2e-7;
c = alpha*(-30*k1^2-30*k2^2)/30;
T = 30*exp(c*t)*sin(k1*centerx)*sin(k2*centery)+40;
t50 = double(solve(T==50,t));
t50days = t50/60/60/24;

%Numerical Solution
%Define the domain
xmin = 0; xmax = 3; numx = 15;
ymin = 0; ymax = 5; numy = 15;
x=linspace(xmin,xmax,numx+2);
dx = x(2)-x(1);
y=linspace(ymin,ymax,numx+2);
dy = y(2)-y(1);
[x,y]=meshgrid(x,y);
%Initial Time
t(1)=0;
%Delta t
dt=1/(alpha*((1/dx^2)+(1/dy^2)))%2000;
%Set up T matrix
T=zeros(10,length(x(:,1)),length(y(:,1))); %Pre allocate space in T(n,length,length) to pre allocate space for n number of time steps
%Initial Condition set up
IC = @(x,y)40+30*sin((pi/3)*x).*sin((pi/5)*y);%70; % 
T(1,:,:) = IC(x,y);
%Set up time stepper
n=1;
%While loop 
while T(n,9,9)>50
    for i=2:length(x)-1
        for j=2:length(y)-1
         
  %Boundary conditions
T(n,1,:)=40; %left boundary
T(n,length(x(:,1)),:)=40; %right boundary
T(n,:,1)=40; %bottom boundary
T(n,:,length(y(:,1))) = 40; %top boundary    

            %Solution at interior points
   T(n+1,i,j) = T(n,i,j)+dt*alpha*(((T(n,i-1,j)-2*T(n,i,j)+T(n,i+1,j))/dx^2)+((T(n,i,j-1)-2*T(n,i,j)+T(n,i,j+1))/dy^2));     
  
            
        end
        
    end
    
    %Time vector
    t(n+1)=t(n)+dt;
    
    %Plot
     surf(reshape(T(n+1,:,:),[length(x),length(y)]),'Edgecolor','none');
     xlabel('Horizontal Distance (m)')
     ylabel('Vertical Distance (m)')
     zlabel('Temperature (Degrees Fahrenheit)')
     title('Two-Dimensional Heat Equation')
      colormap jet

    axis([0,length(x(:,1)),0,length(y(:,1)),40,70])
 
  
    %Time step
   n=n+1;
   %Animate
    drawnow;
    
end
%Ending Temperature
FinalTemp=T(end,9,9)
%Ending time
Finaltimeindays = double(t(end)/60/60/24)
%Analytic and Numerical Time Difference
Difference=double(abs(Finaltimeindays-t50days)) %1/25 difference in days between analytic and numerical



