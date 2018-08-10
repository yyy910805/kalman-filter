% unstaggered grid
% obtains all assimilation -> propagation results for Cascadia
% Ensemble Kalman filter

function out = enkf_cascadia(filename, order, assim_step, fcst_step)
%{
@param filename: name of file that contains parameters
@param order: order of SBP operator
@param assim_step: interval of assimilation run (e.g. assim_step=10 if we 
                   assimilate every 10*dt)
@param fcst_step: interval of forecast run (e.g. fcst_step=100 if we
                  forecast every 100*dt)
@return out: true solution, x, t, observation points, assimilated result 
             matrix (nx by nt by runs), and evolution of standard deviation
%}

%%
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

% number of assimilations to do
assim_runs = floor(nt/assim_step) - 1;

% number of forecasts to do
fcst_runs = floor(nt/fcst_step) - 1;

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

% all forecast results
all_res = zeros(nx,nt,fcst_runs);

% FIXME: intial guess for wave height (need pressure time seris)
N = 100;
mu = zeros(1,nx);
mu_h = zeros(1,nx);
mu_q = zeros(nx,1);
% set initial ensemble
h0 = mvnrnd(mu_h,P,N)';
q0 = repmat(mu_q,1,N);
X = [q0;h0];
std = zeros(assim_runs + 1,nx);
std(1,:) = sqrt(diag(P));

tic
% assimilation for all runs
for i = 1:assim_runs
    % propagate this number of steps
    for k = 1:assim_step
        X = T*X;
    end
    q = X(1:nx,:);
    h = X(1+nx:2*nx,:);
    
    % reconstruct covariance matrix
    P_ = cov(X');
    P = P_(1+nx:2*nx,1+nx:2*nx);
    std(i+1,:) = sqrt(diag(P));
    % update
    h = h + P*H'*linsolve((H*P*H' + R),(repmat(z(:,assim_step*i+1),1,N) + ...
        mvnrnd(mu*H',R,N)' - H*h));
    X(1+nx:2*nx,:) = h;
    
    if mod(i*assim_step,fcst_step) == 0
        % rth forecast run
        r = i*assim_step/fcst_step;   
        % create vector v to store current forecast result
        v = zeros(2*nx,nt);
        % the updated result becomes the initial condition
        v(1:nx,i*assim_step) = mean(q,2);
        v(1+nx:2*nx,i*assim_step) = mean(h,2)';
    
        % propagate and save
        for j = i*assim_step+1:nt
            v(:,j) = T*v(:,j-1);
        end
        all_res(:,i*assim_step:nt,r) = v(nx+1:2*nx,i*assim_step:nt);
        if r > 1
            all_res(:,1:i*assim_step-1,r) = all_res(:,1:i*assim_step-1,r-1);
        end
    end
end
toc

% save all output
out.real = u(nx+1:2*nx,:);
out.all_res = all_res;
out.obs = obs;
out.std = std;
out.x = x;
out.dx = dx;
out.t = t;
out.dt = dt;
out.assim_runs = assim_runs;
out.fcst_runs = fcst_runs;
end