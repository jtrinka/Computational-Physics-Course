function [derivs] = Schrodinger1SpatialD(t,y)
psi=y(1);
dpsi=y(2);
E = y(3); %Energy level

m = 1.6750e-27;%kg
h = 6.626e-34;%J s
hbar = h/(2*pi); %constant

%potential energy function
if(t>1.5e-15)
    
    V=0;
    
else
    
    V=-150*1.60218e-13; %MeV to joules


end

%Second order differential equation
d2psidx2 = (((V-E)*2*m)/(hbar^2))*psi;

[derivs]=[dpsi,d2psidx2,0]';
end

