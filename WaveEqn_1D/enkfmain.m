% unstaggered grid
function out = enkfmain1(filename, order, graph)
%{
@param filename: name of file that contains parameters
@param graph: 0 for no graph, 1 to display only dynamic graph,
              2 to display only initial and final comparison graph
@return out: [u,v,x,t, h_save, v_save, std_save]
             u = true solution
             v = assimilated solution (apply assimilation at every time)
             x, t = spatial and temporal grid
             h_save = saved ensemble height at the coast for time ~ 1500s
             v_save = saved assimilated height at the coast for time ~ 1500s
             std_save = saved std at the coast for time ~ 1500s
%}

%%
% read paramaters from file
[xmin, xmax, tmin, tmax, dx, ~, Ld, xsd, esd, freq] = readParams(filename);

% first get the true solution, for comparison purpose later
[u, dt, T] = wave_solve_unstagg(filename, order, 0);

% set up the grid
x = xmin:dx:xmax;
t = tmin:dt:tmax;
nx = length(x);
nt = length(t);

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

% create vector v to store assimilation result for current time step
v = zeros(1,nx);

%%
% set initial ensemble
N = 100;
mu = zeros(1,nx);
h0 = mvnrnd(mu,P,N)';
q0 = zeros(nx,N);
X = [q0;h0];
v(1,:) = zeros(1,nx);
std = zeros(nt,nx);
std(1,:) = sqrt(diag(P));
h_save = zeros(1, N);

% Ensemble Kalman filter assimilation
for i = 2:nt
    % evolve ensemble forward
    X = T*X;
    h = X(1+nx:2*nx,:);
    % reconstruct covariance matrix
    P_ = cov(X');
    P = P_(1+nx:2*nx,1+nx:2*nx);
    % update 
    h = h + P*H'*linsolve((H*P*H' + R),(repmat(z(:,i),1,N) + ...
        mvnrnd(mu*H',R,N)' - H*h));
    X(1+nx:2*nx,:) = h;
    v(i,:) = mean(h,2)';
    std(i,:) = sqrt(diag(P));
    % save height at the coast for time ~ 1500s
    if i == ceil(nt/2)
        h_save = h(1,:);
        std_save = std(i,:);
        v_save = v(i,:);
    end
    
    if graph == 1
        plot(x,u((nx+1:2*nx),i),x,v(i,:),x,sqrt(diag(P)),x(obs),0*obs,'o');
        xlabel('distance from coast');
        ylabel('wave height');
        legend('true solution','assimilated solution','std');
        legend('Location','northeastoutside'); 
        %pause
        drawnow
    end
end

if graph == 2
    plot(x,u((nx+1:2*nx),i),x,v(i,:),x,sqrt(diag(P)),x(obs),0*obs,'o');
    xlabel('distance from coast (km)','FontSize',16);
    ylabel('wave height (m)','FontSize',16);
    title('Assimilated result vs. true solution','FontSize',16);
    figure
end

% space-time plots
u = u(nx+1:2*nx,:)';
p1 = pcolor(x,t,u);
set(p1,'LineStyle','none');
cmap();
xlabel('Distance from the Coast (km)');
ylabel('Time (s)');
title('Space-Time Plot of True Tsunami Wave Height (m)');
colorbar
caxis([-3 3])
figure

p2 = pcolor(x,t,v);
set(p2,'LineStyle','none');
cmap();
hold on
plot(x(obs),0*obs,'^','MarkerFaceColor','k','MarkerSize',5);
xlabel('Distance from the Coast (km)');
ylabel('Time (s)');
title('Space-Time Plot of Assimilated Tsunami Wave Height (m)');
colorbar
caxis([-3 3])
hold off
figure

% histogram plot for wave height at coast at t ~ 1500s
g = histfit(h_save,50,'kernel');
g(1).FaceColor = [.8 .8 1];
xlabel('Wave Height (m)')
ylabel('Frequency')
title('Probability Distribution of Wave Height at Coast at t = 1500s')
figure

% standard deviation at different times
plot(x,std(round(60/dt),:),'LineWidth',1);
hold on;
plot(x,std(round(300/dt),:),'LineWidth',1);
hold on;
plot(x,std(round(900/dt),:),'LineWidth',1);
hold on;
plot(x,std(round(1800/dt),:),'LineWidth',1);
hold on;
plot(x,std(round(2400/dt),:),'LineWidth',1);
hold on;
plot(x(obs),0*obs,'^','MarkerFaceColor','k','MarkerSize',5);
handle = legend('1min','5min','15min','35min','40min','observations');
set(handle,'FontSize',12);
%set(gca,'fontsize',16);
xlabel('Distance from Coast (km)')
ylabel('Wave Height Standard Deviation (m)')
title('Wave Height Standard Deviation at Various Times')
hold off

% save results in struct
out.u = u;
out.v = v;
out.x = x;
out.t = t;
out.h_save = h_save;
out.v_save = v_save;
out.std_save = std_save;

end