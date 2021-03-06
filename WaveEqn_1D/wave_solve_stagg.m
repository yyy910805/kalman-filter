% using staggered grid, and T operator
% for Cascadia
% solves for true solution of waveform over time given initial wave height
% and ocean bathymetry

function [u, dt, T] = wave_solve_stagg(filename, order, graph)

%{
@param filename: name of parameter file
@param order: Runge-Kutta order
@param graph: 0 for no plot, 1 for dynamic solution plot, 2 for initial and
final plots
%}

%%
% read paramaters from file
[xmin, xmax, tmin, tmax, dx, CFL, ~, ~, ~, ~] = readParams(filename);

% import water depth data
irreg = load('cascadia.mat');
% flip the axis to make coast on the left
x_irreg = -flipud(irreg.x_SF);
x_irreg = x_irreg - x_irreg(1);
y_irreg = flipud(irreg.y_SF); % irregularly sampled y-position [km]

% initial wave height at the moment of rupture on irregular grid
h_irreg = flipud(irreg.uy_os);

% get rid of overlapping locations in the x-position
h_irreg(diff(x_irreg)==0) = [];
y_irreg(diff(x_irreg)==0) = [];
x_irreg(diff(x_irreg)==0) = [];

% interpolation on regular grid
x = xmin:dx:xmax;
nx = length(x);
y = -interp1(x_irreg,y_irreg,x'); % sea floor depth [km]
h_reg = interp1(x_irreg,h_irreg,x'); % wave height

% compute weighting to remove the body waves in Gabe's simulation
weight = zeros(nx,1);
for i = 1:nx
    easting = xmin + dx*(i-1);
    if (easting > 300)
        weight(i,1) = sin(3.14/2*(xmax - easting)/(xmax-300))^10; 
    else
        weight(i,1) = 1;
    end
end
h0 = h_reg.*weight;

%%
% struct to store some essential variables
M.g = 10e-3; % [km/s^2]
M.H = y; % water depth [km]
M.c = sqrt(M.g*M.H); % wave velocity [km/s]

% set grid for t
dt = CFL*dx/max(M.c);
t = tmin:dt:tmax;
nt = length(t);

% true solution u, each column is a time step, stores q and h
u = zeros(2*nx+1,nt);

%%
% Parameters
m = length(x);
xvec = x;
c = @(x) interp1(xvec,M.c,x); % wave speed (km/s)

%==== Assemble coefficient matrices with parameters ====
% This should never change as long as we are 
% running the shallow water equations!
A = {@(x) 1./(c(x).^2), @(x)0*x;...
	 @(x)0*x, @(x) 0*x + 1};
B = {@(x)0*x, @(x)0*x + 1;...
	 @(x)0*x + 1, @(x)0*x};
%=====================================================

% Forcing
F = [];

% Domain
xlim = {xmin, xmax};

% Boundary data (empty for homogeneous BC)
g0 = [];
gL = []; 

% Create discretization object
discr = TsunamiDiscretization(m, order, xlim, F, A, B, g0, gL);

% Grid points
x_primal = discr.diffOp.grid_primal.points();
x_dual = discr.diffOp.grid_dual.points();

% Initial data
q0_fun = @(x) 0*x;
h0_fun = @(x) interp1(xvec,h0,x);
discr.v0 = [q0_fun(x_primal); h0_fun(x_dual)];
u(:,1) = discr.v0;

% get propagation matrix T
[RK,alpha] = rk_init(dt); % initialize RK coefficients
D = discr.D;
T = eye(length(D)) + dt*D + dt^2*D^2/2 + dt^3*D^3/6 + dt^4*D^4/24 + ...
    RK.B(1)*RK.B(2)*RK.B(3)*RK.B(4)*RK.B(5)*D^5;

% solve for true solution iteratively
for i = 2:nt 
    u(:,i) = T*u(:,i-1);
    if graph == 1
        plot(x_dual,u(nx+1:2*nx+1,i));
        drawnow
    end
end

if graph == 2
    plot(x, u(1+nx:2*nx,1),'Linewidth',1),xlabel('distance'),ylabel('height');
    hold on;
    plot(x, u(1+nx:2*nx,nt),'Linewidth',1);
    legend('initial','final');
    legend('Location','northeastoutside');
end

end