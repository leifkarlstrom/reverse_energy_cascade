%     This script solves the 1D unsteady heat equation with source 
%     implicitly with a Crank-Nicolson method 
   
%     Copyright (C) May 2017  Leif Karlstrom leif@uoregon.edu
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
year=3600*24*365; %one year in seconds
h=.1 *1e3; %grid spacing, m
k=100 *year; %time step, years
L = 30 *1e3; %total length in km
N = L/h ; %number of nodes = domain length in km/spacing
M = 10e6*year/k; %total time / timestep
t=0:k:M*k; %time vector
x = 0:h:L; %distance vector
qm=.02 ; %mantle heat flux, W/m2 
q0=.045 ; %surface heat flux
hr=1e4 ;%radioactive production scale, m
rho=2600 ;%density, kg/m3
kc=2 ;%thermal conductivity, W/m K
cp=1200 ;%specific heat capacity, J/kg K
Kappa = kc/(rho*cp); %thermal diffusivity
Lf=310000; %latent heat fusion J/Kg
Tliq=1100; %magma liquidus in C
Tsurf = 25; %surface temperature in C
A=4.25*10^7; %parameters for viscosity law
Bmu=8.31; %gas constant
G=141*10^3; %activation energy
Keff=1e10;%Effective elastic modulus for Deborah number calculation
Eff=.01; %efficiency factor for heating of intrusions
QvolV=1e-3; %long term average magma flux in m3/yr
T=1e4*year; %period of magma injection
Qvol=QvolV*pi*(sin(2*pi*t/T)); %specify flux as positive part of sinusoidal influx 
Qvol(Qvol<0)=0;
IntrWid=20 *1e3; %width of intruded region
IntrTop = L-IntrWid; %top of intuded region
IntrBot=L; %bottom of intruded region
itop=find(x==IntrTop);
ibot=find(x==IntrBot);
%intrusion volume, assume .5 m width with length equal to intruded region
V0=IntrWid^2*.5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solve problem with Crank-Nicolson method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = Kappa*k/(2*h^2); e = ones(N,1);
B = spdiags([e -2*e e], -1:1, N, N); 
%Boundary conditions
B(N,N)=-1; B(1,1)=0; B(1,2)=0;
%assemble matrices
B1 = (eye(N) - v*B); B2 = (eye(N) + v*B);
%solution and parameter grids
u = zeros(N,M); Mu = zeros(N,M); De = zeros(N,M);
%Set initial conditions
H0=(q0-qm)/(rho*hr);
for i = 0:N-1%steady state geotherm
        u(i+1,1) = Tsurf + qm*x(i+1)/kc + rho*H0*hr^2/kc *(1-exp(-x(i+1)/hr)) ;
end
   Mu(:,1) = A.*exp(G./(Bmu*(u(:,1)+273))); %viscosity and Deborah number
   De(:,1) = Keff*V0./(QvolV*1e9/year * Mu(:,1));
%define initial heat production due to diking, distributed over region in crust 
Hmagma=(Qvol(1)*1e9/(year*L^2*IntrWid))*Eff*(cp*(Tliq-u(:,1))+Lf)... 
.*(1+tanh(10*(x(1:end-1)'-IntrTop)/L))/2; %amplitude of heat flux
F=k*(H0/cp*exp(-x(1:end-1)'/hr) + Hmagma); F(1)=0; %heating from radioactivity and magma
for m = 1:M-1%internal points 
    if m>1 %heat production 
Hmagma=(Qvol(m).*1e9./(year*L^2*IntrWid))*Eff*(cp*(Tliq-u(:,m))+Lf)... 
.*(1+tanh(10*(x(1:end-1)'-IntrTop)/L))/2;
F=k*(H0/cp*exp(-x(1:end-1)'/hr) + Hmagma);
%correct last source term to include nonzero boundary heat flux
F(end)=F(end) + 2*v*h*(qm + rho*H0*hr*(exp(-x(end)/hr))) /kc; 
    end
%solve the linear system
    u(:,m+1) = B1\(B2*u(:,m) + F); 
    u(1,m+1)= Tsurf;% inject the top BC
%calculate viscosity and Deborah number
    Mu(:,m+1) = A.*exp(G./(Bmu*(u(:,m+1)+273))); %temp in K
    De(:,m+1) = (QvolV.*1e9./year * Mu(:,m+1))./(Keff.*V0); 
    if De(itop,m+1)<1 
        break %have reached De=1 at top of intruded region, end timestep       
    end        
end
[~,idx]=min(abs(De(itop,:)-1));
disp(['Time until midpoint of intruded region has De=1 is ' num2str(t(idx)/year) ' years'])

    