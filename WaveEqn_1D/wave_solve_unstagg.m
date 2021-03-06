% using unstaggered grid, and T operator
% for Cascadia
% solves for true solution of waveform over time given initial wave height
% and ocean bathymetry

function [u, dt, T] = wave_solve_unstagg(filename, order, graph)

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
u = zeros(2*nx,nt);
u((nx+1):2*nx,1) = h0; % initial condition for wave height

%%
% SBP operators
[D,Hinv,~,A,~,~,~] = SBPoperatorsRussian(nx-1,dx,order);
cdiss = 0.25; % artificial dissipation term
SAT = [M.c(1) M.c(end)]*Hinv; % boundary conditions
[RK,alpha] = rk_init(dt); % initialize RK coefficients

T = compute_T_matrix(nx-1,alpha,M,D,SAT,A,cdiss);
T = sparse(T);

% solve for true solution iteratively
for i = 2:nt
    u(:,i) = T*u(:,i-1);
    if graph == 1
        plot(x,u(1+nx:2*nx,i));
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