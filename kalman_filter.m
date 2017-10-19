function v = kalman_filter(xmin, xmax, tmin, tmax, dx, dt, c, esd, freq)
%
% @param xmin: minimum value of x
% @param xmax: maximum value of x in the period
% @param tmin: minimum value of time
% @param tmax: maximum time over which the assimilation is done
% @param dx: step size of x
% @param dt: step size of t
% @param c: speed of propagation
% @param Ld: spatial decorrelation length
% @param xsd: initial estimate error standard deviation
% @param esd: measurement error standard deviation
% @param freq: frequency of observation/measurement
%
% @return v: assimilated data stored in a matrix, with each row representing
%            one time step, and each column a location in space
%

% first get the true solution, for comparison purpose later
u = linear_1d(xmin, xmax, tmin, tmax, dx, dt, c);

% setting up the grid
x = xmin:dx:xmax;
t = tmin:dt:tmax;
nx = length(x);
nt = length(t);

% create matrix v to store assimilation results
v = zeros(nt,nx);
v(1,:) = zeros(1,nx); % set initial guess to zero

% set matrix A for the process
A = (1-c*dt/dx)*eye(nx);
A(1,nx) = c*dt/dx;
for i = 2:nx
    A(i,i-1) = c*dt/dx;
end

% set observation(measurement) matrix that comes from true state plus a
% random noise ~N(0,esd^2)
% number of observations is smaller than or equal to number of states
% depending on the frequency of observation
obs = 1:freq:nx;
m = length(obs);
H = zeros(m,nx);
H(:,obs(:)) = eye(m);
z = H*u' + esd*randn(m,nt);
% set observation error covariance matrix
R = esd^2*eye(m);

% initial prior estimate error covariance matrix (diff between true state
% and a priori state estimate)
P = correl(nx, dx, Ld); % Gaussian correlation matrix
P_ = xsd^2*P; % covariance matrix

% Kalman Filter assimilation through time
for n = 2:nt
    % time update
    x_ = A*v(n-1,:)';
    P_ = A*P_*A';
    
    % measurement update
    K = P_*H'*pinv(H*P_*H' + R);
    x_hat = x_ + K*(z(:,n) - H*x_);
    P_ = (eye(nx) - K*H)*P_;
    
    % store assimilated solution in v
    v(n,:) = x_hat';
    
    plot(x,u(n,:),x,v(n,:))
    drawnow
    
end

subplot(2,1,1);
plot(x, v(1,:),'Linewidth',1),xlabel('x'),ylabel('solution');
hold on;
plot(x, u(1,:),'Linewidth',1),legend('initial guess','true initial condition');

subplot(2,1,2);
plot(x, v(nt,:),'Linewidth',1),xlabel('x'),ylabel('solution');
hold on;
plot(x, u(nt,:),'Linewidth',1),legend('final assimilation result','true final state');
legend('Location','northwest');

figure
pcolor(x,t,v),shading flat,colorbar;

end