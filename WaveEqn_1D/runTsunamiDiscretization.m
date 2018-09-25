function u = runTsunamiDiscretization(m,order,plot_flag)
default_arg('m',101);
default_arg('order',4)
default_arg('plot_flag',true)

% Parameters
L = 1000;
g = 10e-3; % km/s^2
H = @(x) 1 + 6/L*x;

% High-frequency coefficient
% H = @(x) 1 + 0.9*sin(50*2*pi*x/L);
c = @(x) sqrt(g*H(x)); % wave speed (km/s)

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
xlim = {0, L};

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
h0_fun = @(x) exp(-0.5*((x-L/2)/(0.1*L)).^2);
discr.v0 = [q0_fun(x_primal); h0_fun(x_dual)];

% Evaluate wave speed
c_primal = c(x_primal);

% Total simulation time
T = 4*L/min(c_primal);

% Solve
if plot_flag
	animationSpeed = 5000;
	noname.animate(discr, animationSpeed, T);
	e = [];
	u = [];
else
	[ts,N] = discr.getTimestepper([],T);
	fprintf('Computing o = %d, m = %d.  \n',order,m(1));
	ts.stepN(N,true);
	[u, ~] = ts.getV;
end

% Return meshsize
% h = (xr-xl)/(m-1);