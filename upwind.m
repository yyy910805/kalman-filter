function Tnp1 = upwind (Tn, t)

%UPWIND  Uses an upwind difference scheme to solve 
%        the advection-diffusion equation
%
%   Global parameters:
%       Courant_Number:  u  * dt / dx
%       Diffusivity:     mu * dt / (dx^2)
%
%  where u is the advective speed; mu is the physical diffusivity
%  parameter; dx is the grip point distance; and dt is the time
%  interval.
%
%  Notice that stability requires:  
%
%      dt <=  0.5 * dx^2 / ( 1 + 0.5 * R )
%
%  where R = Courant_Number / Diffusivity
%
%Ricardo Todling, 24 Jan 98    Initial code.

global  Courant_Number Diffusivity

C = Courant_Number;
s = Diffusivity;
jdim = size(Tn,1);
Tnp1 = zeros(jdim,1);
%
%  Loop over the central domain
%
for j = 2:jdim-1
    Tnp1(j) = (s+C)*Tn(j-1) + (1-2*s-C)*Tn(j) + s*Tn(j+1);
end 
%
%  Update periodic boundaries
%
j = 1;
Tnp1(j)  = (s+C)*Tn(jdim) + (1-2*s-C)*Tn(j) + s*Tn(j+1);
j = jdim;
Tnp1(j)  = (s+C)*Tn(j-1)  + (1-2*s-C)*Tn(j) + s*Tn(1); 