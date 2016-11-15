%% Microbe Growth
clear; clc;
xmin = -20;
xmax = 20;
tmin = 0;
tmax = 6*60*60; %6 hours
npt= 100;
%Build domain
x = linspace(xmin,xmax,npt+2);
dx = x(2)-x(1);
dt=30;
%Number of Time Steps
ntimesteps = tmax/dt;
%Constants converted into seconds
k=1/300;
r=1/600;
%Constants calculated from algebra
p = (dt*k)/(2*dx^2);
q = dt*r;
%Initial condition is a Gaussian
IC = @(x) .16*exp(-x.^2);
%Space pre-allocation
U = zeros(length(x),ntimesteps);
%Initial Condition
U(:,1) = IC(x);
%Dirichlet Boundary Conditions on first step
U(1,1) = 0;
U(end,1) =0;
%Space pre-allocation
A = zeros(npt,npt);
b = zeros(npt,1);
%A never changes so just set it up outside of time stepper
for i=1:npt-1
A(i,i) = 1+2*p-(q/2);
A(i+1,i) =-p;
A(i,i+1) =-p;
end
A(end,end)=1+2*p;
%Take A inverse
Ainverse=inv(A);
%Time step
for n=1:ntimesteps
    %Use b to find interior point solutions so start at space index 2.
    b(1) = (1-2*p+(q/2))*U(2,n)+p*U(3,n);
    %Go till final space index-1
    b(end) = (1-2*p+(q/2))*U(end-1,n)+p*U(end-2,n);
    for i=2:npt
        %Loop over all the other interior points
   b(i) = p*U(i+2,n)+(1-2*p+(q/2))*U(i+1,n)+p*U(i,n);
    end
    
    
    
    
   %Only put interior points here and left multiply
    U(2:end-1,n+1)=Ainverse*b;

    
     plot(x,U(:,n+1))
     xlabel('meters')
     ylabel('Distribution of Microbes')
     title(n)
   axis([-20,20,0,4.5e13])
   drawnow
    
    
    
end
    

