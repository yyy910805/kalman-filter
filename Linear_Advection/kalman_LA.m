function v = kalman_LA(filename, graph, filter)

%{
@param filename: name of file that contains parameters
@param graph: 0 for no graph, 1 to display only dynamic graph,
              2 to display only initial and final comparison graph
@param filter: 1 to use dynamic filter, 2 to set Kalman gain matrix fixed
@return v: assimilated result matrix
%}

% read paramaters from file
[xmin, xmax, tmin, tmax, dx, dt, c, Ld, xsd, esd, freq] = readParams(filename);

% first get the true solution, for comparison purpose later
[u, g] = bc1linear_1d(xmin, xmax, tmin, tmax, dx, dt, c, 0);

% set up the grid
x = xmin:dx:(xmax-dx);
t = tmin:dt:tmax;
nx = length(x);
nt = length(t);

% create matrix v to store assimilation results
v = zeros(nt,nx);
% set initial guess
v(1,:) = zeros(1,nx);
% set boundary condition
v(:,1) = g;

% set process matrix A for periodic boundary condition
%{
A = (1-c*dt/dx)*eye(nx);
A(1,nx) = c*dt/dx;
for i = 2:nx
    A(i,i-1) = c*dt/dx;
end
%}

% set process matrix A for non-periodic boundary condition
e0 = zeros(nx,1);
e0(1) = 1;
M = eye(nx);
M(1,1) = -1;
for i = 2:nx
    M(i,i-1) = -1;
end
M = c*M/dx;

A = eye(nx) - M*dt - 2*dt/dx*e0*(e0');

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
z = H*u' + esd*randn(m,nt);
% set observation error covariance matrix
R = esd^2*eye(m);

% initial prior estimate error covariance matrix (diff between true state
% and a priori state estimate)
P = gaussian(nx, dx, Ld); % Gaussian correlation matrix
P_ = xsd^2*P; % covariance matrix

if filter == 2
    K = P_*H'*pinv(H*P_*H' + R);
end
tic
% Kalman Filter assimilation through time
for n = 2:nt
    % time update
    x_ = A*v(n-1,:)' + 2*e0*dt/dx*g;
    if filter == 1
        P_ = A*P_*A';
        % measurement update
        K = P_*H'*pinv(H*P_*H' + R);
        P_ = (eye(nx) - K*H)*P_;
    end
    x_hat = x_ + K*(z(:,n) - H*x_);
    
    % store assimilated solution in v
    v(n,:) = x_hat';
    
    if graph == 1
        plot(x,u(n,:),x,v(n,:),x,sqrt(diag(P_)),x(obs),0*obs,'o');
        xlabel('x');
        ylabel('y');
        legend('true solution','assimilated solution','std of error');
        legend('Location','northeastoutside'); 
        %pause
        drawnow
    end
end
toc

if graph == 2
    subplot(2,1,1);
    plot(x, v(1,:),'Linewidth',1),xlabel('x'),ylabel('solution');
    hold on;
    plot(x, u(1,:),'Linewidth',1),legend('initial guess','true initial condition');
    subplot(2,1,2);
    plot(x, v(nt,:),'Linewidth',1),xlabel('x'),ylabel('solution');
    hold on;
    plot(x, u(nt,:),'Linewidth',1),legend('final assimilation result','true final state');
    legend('Location','northeastoutside');
    figure
    pcolor(x,t,v),shading flat,colorbar;
end

end