%Spring Mass Simulation
clear; clc; clf; close all;
g = 9.81;
m = 1/2;
x = 9/50;
k = (m*g)/x;
tmin=0;
tmax=2;
n=100000;
t= linspace(tmin,tmax,n);
t=t';
h= t(2)-t(1);
aprime = @(t,a,b) b;
bprime = @(t,a,b) -(k/m).*a;
a(1) = -.07;
b(1) = 0;
energy(1) = .5*m*b(1).^2+.5*k*a(1).^2;
tic
for i=1:1:length(t)-1
   %k1 at start
   k1a =  aprime(t(i),a(i),b(i));
   k1b =  bprime(t(i),a(i),b(i));
   %k2 at a half step
   k2a = aprime(t(i)+h/2,a(i) + (h/2)*k1a, b(i)+(h/2)*k1b);
   k2b = bprime(t(i)+h/2,a(i) + (h/2)*k1a, b(i)+(h/2)*k1b);
   %k3 at a half step
   k3a = aprime(t(i)+h/2,a(i) + (h/2)*k2a, b(i)+(h/2)*k2b);
   k3b = bprime(t(i)+h/2,a(i) + (h/2)*k2a, b(i)+(h/2)*k2b);
   %k4 at a full step
   k4a = aprime(t(i)+h,a(i)+h*k3a,b(i)+h*k3b);
   k4b = bprime(t(i)+h,a(i)+h*k3a,b(i)+h*k3b);
   
   a(i+1) = a(i)+(h/6)*(k1a+2*k2a+2*k3a+k4a);
   b(i+1) = b(i)+(h/6)*(k1b+2*k2b+2*k3b+k4b);
   
   energy(i+1)=(1/2)*m*b(i+1).^2+.5*k*a(i+1).^2;
    
    
    
    
    
end
toc
sol = @(t) -.07.*cos(sqrt(k/m).*t);
subplot(1,2,1)
plot(t,a)
xlabel('Time (seconds)')
ylabel('Distance from Equilibrium (meters)')
title('Numerical Solution to the Spring-Mass System')
subplot(1,2,2)
plot(t,sol(t))
xlabel('Time (seconds)')
ylabel('Distance from Equilibrium (meters)')
title('Analytical Solution to the Spring-Mass System')
figure
plot(t,a)
xlabel('Time (seconds)')
ylabel('Distance from Equilibrium (meters)')
title('Analytical and Numerical Solutions to the Spring-Mass System')
hold on
plot(t,sol(t))
a=a';
figure
plot(t,energy)
xlabel('Time (seconds)')
ylabel('Total Mechanical Energy of the System (Joules)')
title('Mechanical Energy Conservation of the Spring-Mass System')
error=norm(sol(t)-a)






