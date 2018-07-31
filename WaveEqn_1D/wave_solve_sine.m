% using unstaggered grid, and T operator
% for a general sine wave, uniform ocean depth
% examines relation between wavelength and station spacing
% solves for true solution of waveform over time given initial wave height

function [u, dt, T] = wave_solve_sine(filename, lambda, order, graph)

%{
@param filename: name of parameter file
@param lambda: wavelength of since wave
@param order: Runge-Kutta order
@param graph: 0 for no plot, 1 for dynamic solution plot, 2 for initial and
final plots
%}

%%
% read paramaters from file
[xmin, xmax, tmin, tmax, dx, CFL, ~, ~, ~, ~] = readParams(filename);

% set up initial wave height and ocean depth on the grid
x = xmin:dx:xmax;
nx = length(x);
y = ones(nx,1)*2; % sea floor depth [km]
h = @(s) sin(2*pi/lambda*s).*(s >= (xmax-lambda)/2 & s <= (xmax+lambda)/2); % initial wave height
h0 = h(x)';

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