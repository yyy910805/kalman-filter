% unstaggered grid
% obtains all assimilation -> propagation results for Cascadia
% Ensemble Kalman filter

function out = optinterp_cascadia(filename, order, step)
%{
@param filename: name of file that contains parameters
@param order: order of SBP operator
@param step: interval of assimilation run (e.g. step=10 if we assimilate
             every 10*dt)
@return out: true solution, x, t, observation points, assimilated result 
             matrix (nx by nt by runs)
%}

%%
addpath WaveEqn_1D
% read paramaters from file
[xmin, xmax, tmin, tmax, dx, ~, Ld, xsd, esd, freq] = readParams(filename);

% first get the true solution, for comparison purpose later
% u stores true solution (2*nx by nt)
% T is the propagation matrix
[u, dt, T] = wave_solve_unstagg(filename, order, 0);

% set up the grid
x = xmin:dx:xmax;
t = tmin:dt:tmax;
nx = length(x);
nt = length(t);

% value of max wave height at coast, and arrival time index
[val,idx] = max(u(nx+1,:));
% only run assimilation up to when the highest wave arrives
%runs = ceil(idx/step);
runs = floor(nt/step) - 1;

%%
%{
set observation(measurement) matrix that comes from true state plus a
random noise ~N(0,esd^2)
number of observations is smaller than or equal to number of states
depending on the frequency of observation
%}
obs = freq:freq:nx;
m = length(obs);
H = zeros(m,nx);
H(:,obs(:)) = eye(m);
z = H*u((nx+1:2*nx),:) + esd*randn(m,nt);

% set observation error covariance matrix
R = esd^2*eye(m);

% initial prior estimate error covariance matrix (diff between true state
% and a priori state estimate)
P = xsd^2*gaussian(nx, dx, Ld); % Gaussian covariance matrix
% Constant Kalman gain matrix
K = P*H'*pinv(H*P*H' + R);

% FIXME: intial guess for wave height (need pressure time seris)
h0 = zeros(nx,1);
q0 = zeros(nx,1);

% all assimilation run results
all_res = zeros(nx,nt,runs);

% initial condition
X = [q0;h0];
    
% assimilation for all runs
tic
for i = 1:runs
    % create vector v to store assimilation result for current assimilation
    v = zeros(2*nx,nt);
    
    % propagate this number of steps
    for k = 1:step
        X = T*X;
    end
    h = X(1+nx:2*nx);

    % update 
    h = h + K*(z(:,i*step+1) - H*h);
    X(1+nx:2*nx) = h;
    v(1+nx:2*nx,i*step) = h;
    % propagate and save
    for j = i*step+1:nt
        v(:,j) = T*v(:,j-1);
    end
    % only save the wave height in final result
    all_res(:,:,i) = v(nx+1:2*nx,:);
    
end
toc

% save output
out.real = u(nx+1:2*nx,:);
out.all_res = all_res;
out.obs = obs;
out.x = x;
out.t = t;
end