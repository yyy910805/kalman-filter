% code for 1D linear shallow water waves,
% written as first-order system of PDEs for
% q = H*u (flux) and h = surface height
% (H = water depth, u = water velocity)
% unstaggered grid
% for Cascadia, set xmax = 300km

function out = wave1d_assimilate(q00,h00,H,x,order,plot_solution,operator)

% q00, h00 = initial conditions of flux and wave height
% H = bathymetry (water depth)
% x = x grid
% order = order of accuracy of difference operators
% plot_solution = true to plot solution at all time steps
% operator = 'staggered' or 'unstaggered'
% out = data structure containing solution, parameters, etc.

% by default, do not plot solution
if nargin==1, order = 8; plot_solution = false; end
if nargin==2, plot_solution = false; end
N = length(x)-1;
L = x(end);
dx = x(2)-x(1);

% SBP grid and differentiation matrix
% D = first derivative
% x = coordinates
% Hinv = corner entry of diagonal norm (for SAT penalties)
% A = artificial dissipation operator
% order = interior order of accuracy
%[D,~,~,Hinv,A] = SBPoperators(N,L,order); % sparse arrays

switch operator
  case 'staggered'
    addpath cascadia/tsunami1d/staggered
    [xp,xm,Pp,Pm,Qp,Qm] = sbp_half(order,N,dx,true);
    Dp = inv(Pp)*Qp; Dm = inv(Pm)*Qm; Hinv = 1/Pp(1,1);
    Np = N+1; Nm = N+2;
  case 'unstaggered'
    addpath cascadia/tsunami1d/SBPoperators
    xm = x; xp = x;
    Np = N+1; Nm = N+1;
    [D,Hinv,~,A,~,~,~] = SBPoperatorsRussian(N,dx,order);
    Dp = D; Dm = D;
end

% material properties (stored in data structure M)
% (no specific units assumed, so either nondimensionalize
% or use SI or other self-consistent set of units)
%
% M.H = water depth
% M.g = gravitational acceleration
% M.c = wave speed

M.g = 9.8e-3; % km/s^2
M.H = interp1(x,H,xm'); % water depth (m)
M.c = sqrt(M.g*M.H); % wave speed (km/s)

% total simulation time
tmax = 4*L/min(M.c);

% time step
CFL = 0.5; dt = CFL*dx/max(M.c);

% SAT penalties
SAT = [M.c(1) M.c(end)]*Hinv; % penalty relaxation rate

% Runge-Kutta coefficients (4th order, low storage method)
RK.nstage = 5;
RK.A = [ 0d0, -567301805773d0/1357537059087d0, -2404267990393d0/2016746695238d0, ...
	 -3550918686646d0/2091501179385d0, -1275806237668d0/842570457699d0 ];
RK.B = [ 1432997174477d0/9575080441755d0, 5161836677717d0/13612068292357d0, ...
	 1720146321549d0/2090206949498d0, 3134564353537d0/4481467310338d0, ...
	 2277821191437d0/14882151754819d0 ];
RK.C = [ 0d0, 1432997174477d0/9575080441755d0, 2526269341429d0/6820363962896d0, ...
	 2006345519317d0/3224310063776d0, 2802321613138d0/2924317926251d0 ];

% solution is computed for following fields:
% q = flux
% h = surface displacement

% rates arrays
q0 = interp1(x,q00,xm');
h0 = interp1(x,h00,xp');
Dq = zeros(Nm,1); % dq/dt
Dh = zeros(Np,1); % dh/dt

% fields arrays
q=q0; h=h0;

% number of time steps
nt = ceil(tmax/dt);
% adjust dt to finish exactly at t=tmax
dt = tmax/nt;

% storage arrays for all time steps
out.q = zeros(Nm,nt+1);
out.h = zeros(Np,nt+1);

% save initial conditions
out.q(:,1) = q0;
out.h(:,1) = h0;

% save space and time vectors
out.t = (0:nt)*dt;
out.xp = xp; out.xm = xm;

% loop over time steps

t = 0; % start at t=0

for m=1:nt

  t0 = t;
  
  for k=1:RK.nstage

    % time
    t = t0+RK.C(k)*dt;

    % scale rates
    Dq = RK.A(k)*Dq;
    Dh = RK.A(k)*Dh;

    % set rates
    [Dq_new Dh_new] = RHS(t,q,h,M,Dp,Dm,SAT);

    % artificial dissipation
    Cdiss = 0; % set to 0 for no artificial dissipation
    if Cdiss>0
      Dq_new = Dq_new+Cdiss*A*q;
      Dh_new = Dh_new+Cdiss*A*h;
    end      

    % add rates to old rates
    Dq = Dq+Dq_new;
    Dh = Dh+Dh_new;

    % update fields
    q = q+dt*RK.B(k)*Dq;
    h = h+dt*RK.B(k)*Dh;

  end

  t = t0+dt;

  % save fields
  out.q(:,m+1) = q;
  out.h(:,m+1) = h;

  if plot_solution
    if mod(m,1)==0
      plot(xp,h)
      title(t)
      drawnow
    end
  else
    if mod(m,100)==0
      disp(['t = ' num2str(t) ' out of ' num2str(tmax)])
    end
  end
    
end

% return other variables as well

out.M = M;
out.L = L;
out.N = N;
out.dt = dt;
out.nt = nt;