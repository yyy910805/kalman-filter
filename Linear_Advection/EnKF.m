function v = EnKF(filename,graph)

%{
@param filename: name of file that contains parameters
@param graph: 0 for no graph, 1 to display only dynamic graph,
              2 to display only initial and final comparison graph
@return v: assimilated result matrix
%}

% read paramaters from file
[xmin,xmax,tmin,tmax,dx,dt,c,Ld,xsd,esd,freq] = readParams(filename);

% first get the true solution, for comparison purpose later
[u,g] = bc1linear_1d(xmin,xmax,tmin,tmax,dx,dt,c,0);

% set up the grid
x = xmin:dx:(xmax-dx);
t = tmin:dt:tmax;
nx = length(x);
nt = length(t);

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
z = H*u'; % true solution at observation locations
% set observation error covariance matrix
R = esd^2*eye(m);

% initial prior estimate error covariance matrix (diff between true state
% and a priori state estimate)
P = gaussian(nx,dx,Ld); % Gaussian correlation matrix
P_ = xsd^2*P; % covariance matrix

% model error covariance matrix
Q = 0.01*eye(nx);

% create matrix v to store assimilation results
v = zeros(nt,nx);

% ensemble Kalman filter
N = 80;  % number of samples for each ensemble
mu = zeros(1,nx); % set initial mean to 0
v(1,:) = mu;
E0 = zeros(nx,N);
E0(1,:) = ones(1,N);
X0 = mvnrnd(mu,P_,N)'; % initial ensemble, each column is a sample
X = X0;
tic
for i = 2:nt
    % evolve each ensemble member forward
    X = A*X + 2*E0*dt/dx*g;
    % reconstruct covariance matrix
    P_ = cov(X');
    % update 
    X = X + P_*H'*linsolve((H*P_*H' + R),(repmat(z(:,i),1,N) + ...
        mvnrnd(mu*H',R,N)' - H*X));
    v(i,:) = mean(X,2)';
    
    if graph == 1
        plot(x,u(i,:),x,v(i,:),x,sqrt(diag(P_)),x(obs),0*obs,'o');
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
    legend('Location','northeastoutside');
    subplot(2,1,2);
    plot(x, v(nt,:),'Linewidth',1),xlabel('x'),ylabel('solution');
    hold on;
    plot(x, u(nt,:),'Linewidth',1),legend('final assimilation result','true final state');
    legend('Location','northeastoutside');
    figure
    pcolor(x,t,v), shading flat, colorbar;
end

end

% Conjugate gradient code for computing lambda
% gamma*lambda = nu for each column of nu
%{
    gamma = H*P_*H' + R;
    nu = repmat(z(:,i),1,N) + mvnrnd(mu*H',R,N)' - H*X;
    L = zeros(m,N);
    eps = 10e-5;
    for j = 1:N
        lambda = zeros(m,1);
        r = nu(:,j);
        p = r;
        while norm(gamma*lambda - r,2) > eps
            alpha = r'*r/(p'*gamma*p);
            lambda = lambda + alpha*p;
            r_ = r - alpha*gamma*p;
            beta = r_'*r_/(r'*r);
            r = r_;
            p = r + beta*p;
        end
        L(:,j) = lambda;
    end
    X = X + P_*H'*L;
    %}
