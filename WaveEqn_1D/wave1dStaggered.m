% code for 1D linear shallow water waves,
% written as first-order system of PDEs for
% q = H*u (flux) and h = surface height
% (H = water depth, u = water velocity)
% staggered grid

function out=wave1dStaggered(N,order,plot_solution)

addpath sbplib
addpath staggered_acoustics_1d

% N = number of grid points
% order = order of accuracy of difference operators
% plot_solution = true to plot solution at all time steps
% out = data structure containing solution, parameters, etc.

% by default, do not plot solution
if nargin==1, order = 8; plot_solution = false; end
if nargin==2, plot_solution = false; end

L = 1000; % length of domain (km)

% SBP grid and differentiation matrix
% D = first derivative
% x = coordinates
% Hinv = corner entry of diagonal norm (for SAT penalties)
% A = artificial dissipation operator
% order = interior order of accuracy

% Get SBP operators
ops = sbp.D1StaggeredUpwind(N+1, {0, L}, order);
dx = ops.h;

H_primal = ops.H_primal;
H_dual = ops.H_dual;
H_primal_inv = inv(H_primal);
H_dual_inv = inv(H_dual);

D1_primal = ops.D1_primal;
D1_dual = ops.D1_dual;

% Create staggered grid using sbplib
[g_primal, g_dual] = grid.primalDual1D(N+1, {0, L});
x_primal = g_primal.points();
x_dual = g_dual.points();

% x_primal: "primal" grid. Completely equidistant.
% x_dual: "dual" grid. Completely equidistant except near boundaries.
% Primal grid has N+1 points, dual has N+2.

% material properties (stored in data structure M)
% (no specific units assumed, so either nondimensionalize
% or use SI or other self-consistent set of units)
% M.L = Domain length
% M.H = water depth
% M.g = gravitational acceleration
% M.c = wave speed

M.L = L;
M.g = 10e-3; % km/s^2
% M.H = linspace(1,7,N+1)'; % water depth (m)
% Create function handle for H that matches the above
M.H = @(x) 1 + 6/L*x;

% High-frequency coefficient
% M.H = @(x) 1 + 0.9*sin(50*2*pi*x/L);

M.c = @(x) sqrt(M.g*M.H(x)); % wave speed (km/s)

% Evaluate functions
M.c_primal = M.c(x_primal);

% total simulation time
tmax = 4*L/min(M.c_primal);

% time step
CFL = 0.5; dt = CFL*dx/max(M.c_primal);

% SAT penalties
SAT_primal = [M.c(0) M.c(L)]*H_primal_inv(1,1); % penalty relaxation rate
SAT_dual = [M.c(0) M.c(L)]*H_dual_inv(1,1); % penalty relaxation rate

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

% Let q live on primal grid and h on dual
% initial conditions
q0_fun = @(x) 0*x;
q0 = q0_fun(x_primal);

h0_fun = @(x) exp(-0.5*((x-L/2)/(0.1*L)).^2);
h0 = h0_fun(x_dual);

% rates arrays
Dq = zeros(size(q0)); % dq/dt
Dh = zeros(size(h0)); % dh/dt

% fields arrays
q=q0; h=h0;

% number of time steps
nt = ceil(tmax/dt);
% adjust dt to finish exactly at t=tmax
dt = tmax/nt;

% storage arrays for all time steps
out.q = zeros( length(q) ,nt+1);
out.h = zeros( length(h) ,nt+1);

% save initial conditions
out.q(:,1) = q0;
out.h(:,1) = h0;

% save space and time vectors
out.t = (0:nt)*dt;
out.x_primal = x_primal;
out.x_dual = x_dual;

% loop over time steps

t = 0; % start at t=0

%=== Setup plot ======
figure;
handle_h = plot(x_dual,h,'b-o');
title(['t = ' num2str(t)])
ylim([-1 1])
hold on
handle_q = plot(x_primal,q,'r-o');
hold off
legend({'h','q'})
%======================

for m=1:nt

  t0 = t;
  
  for k=1:RK.nstage

    % time
    t = t0+RK.C(k)*dt;

    % scale rates
    Dq = RK.A(k)*Dq;
    Dh = RK.A(k)*Dh;

    % set rates
    [Dq_new, Dh_new] = RHS_staggered(q,h,M,D1_primal,D1_dual,SAT_primal,SAT_dual);   

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
    if mod(m,2)==0
       handle_h.YData = h;
       handle_q.YData = q;
       title(['t = ' num2str(t)])
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