% unstaggered grid
function [u,v] = kfmain1(filename, order, graph, filter)
%{
@param filename: name of file that contains parameters
@param graph: 0 for no graph, 1 to display only dynamic graph,
              2 to display only initial and final comparison graph
@param filter: 1 to use dynamic filter, 2 to set Kalman gain matrix fixed
@return v: assimilated result matrix
%}

%%
% read paramaters from file
[xmin, xmax, tmin, tmax, dx, ~, Ld, xsd, esd, freq] = readParams2(filename);

% first get the true solution, for comparison purpose later
[u, dt, T] = wave_solve1(filename, order, 0);

% set up the grid
x = xmin:dx:xmax;
t = tmin:dt:tmax;
nx = length(x);
nt = length(t);

% create matrix v to store assimilation results
v = zeros(2*nx,nt);
% set initial guess of wave height
v((nx+1:2*nx),1) = zeros(nx,1);

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

%%
% no update of Kalman gain matrix
if filter == 2
    K = P*H'*pinv(H*P*H' + R);
end

% Kalman Filter assimilation through time
for n = 2:nt
    % time update
    qh = T*v(:,n-1);
    h = qh(nx+1:2*nx);
    
    if filter == 1
        P = T(nx+1:2*nx,nx+1:2*nx)*P*T(nx+1:2*nx,nx+1:2*nx)';  % not sure about structure of covariance matrix for qh
        % measurement update
        K = P*H'*pinv(H*P*H' + R);
        P = (eye(nx) - K*H)*P;
    end
    h = h + K*(z(:,n) - H*h);
    
    % store assimilated solution in v
    qh(nx+1:2*nx) = h;
    v(:,n) = qh;
    
    if graph == 1
        plot(x,u((nx+1:2*nx),n),x,v((nx+1:2*nx),n),x,sqrt(diag(P)),x(obs),0*obs,'o');
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
    plot(x, v((nx+1:2*nx),1),'Linewidth',1),xlabel('x'),ylabel('solution');
    hold on;
    plot(x, u((nx+1:2*nx),1),'Linewidth',1),legend('initial guess','true initial condition');
    subplot(2,1,2);
    plot(x, v((nx+1:2*nx),nt),'Linewidth',1),xlabel('x'),ylabel('solution');
    hold on;
    plot(x, u((nx+1:2*nx),nt),'Linewidth',1),legend('final assimilation result','true final state');
    legend('Location','northeastoutside');
    figure
    pcolor(x,t,v),shading flat,colorbar;
end

for j = 1:nt
    plot(x,u((nx+1:2*nx),j),'Linewidth',1);
    hold on;
    plot(x,v((nx+1:2*nx),j),'Linewidth',1);
    hold on;
    plot(x,u((nx+1:2*nx),j)-v((nx+1:2*nx),j),'Linewidth',1);
    title(dt*j);
    legend('true solution', 'assimilated result','difference');
    hold off;
    drawnow
end

end