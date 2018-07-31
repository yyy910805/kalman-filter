%{
runs the enkf_cascadia results using optimal interpolation
analyzes the result in terms of:
 - first time the tsunami wave hits the coast
 - max height of tsunami and its arrival time at the coast
 - accuracy of prediction
%}

%% get outputs
addpath 'EnKF vs OptInterp'
order = 4;
assim_step = 1; % assimilate observations every x time steps
fcst_step = 20; % forecast every x time steps
out = optinterp_cascadia('params.txt',order,assim_step,fcst_step);

h = out.real;
all_res = out.all_res;
obs = out.obs;
x = out.x;
dx = out.dx;
t = out.t;
dt = out.dt;
assim_runs = out.assim_runs;
fcst_runs = out.fcst_runs;

nx = length(x);
nt = length(t);

% real prediction of max wave height and its arrival time
real_fcst = zeros(1,2);
[val,idx] = max(h(1,:));
real_fcst(1) = val;
real_fcst(2) = t(idx);

% array to store prediction results for the coast
% col 1 = max wave height
% col 2 = arrival time of max wave height
pred = zeros(fcst_runs,2);
for i = 1:fcst_runs
    [val,idx] = max(all_res(1,:,i));
    pred(i,1) = val;
    pred(i,2) = idx;
end

%% plot prediction over runs
% plot only until the actual max wave hits the coast
arrival_run = ceil(real_fcst(2)/dt/fcst_step) + 1;
ax = (1:arrival_run)*dt*fcst_step;
figure(1)
plot(ax, pred(1:arrival_run,1))
ylim([0 3.5])
hold on
plot([1,arrival_run]*dt*fcst_step,[real_fcst(1),real_fcst(1)],'--')
legend('Forecasted max wave height','True max wave height','Location','Best')
xlabel('Time (s)')
ylabel('Wave Height (m)')
title('Forecasted vs. True Max Wave Height')
hold off

figure(2)
plot(ax, pred(1:arrival_run,2))
ylim([0 1500])
hold on
plot([1,arrival_run]*dt*fcst_step,[real_fcst(2),real_fcst(2)],'--')
legend('Forecasted max wave height arrival time',...
       'True max wave height arrival time','Location','Best')
xlabel('Time (s)')
ylabel('Max Wave Height Arrival Time (s)')
title('Forecasted vs. True Arrival Time of Max Wave Height')
hold off

%% plot accuracy indicators over runs
% true vs assimilated waveform before each forecast
% over the entire grid, or at observation stations
K = zeros(1,fcst_runs);  % mean (0.8<K<1.2 is good)
kappa = zeros(1,fcst_runs);  % std (kappa<1.4 is good)
N = nx;  % over the entire grid
%N = length(obs);  % at observation stations

for i = 1:fcst_runs
    s = 0;
    std = 0;
    for j = 1:N
        s = s + log(abs(h(j,i*fcst_step)/all_res(j,i*fcst_step,i)));
        % s = s + log(abs(h(obs(j),i*fcst_step)/all_res(obs(j),i*fcst_step,i)));
    end
    K(i) = exp(s/N);
    for j = 1:N
        std = std + (log(abs(h(j,i*fcst_step)/all_res(j,i*fcst_step,i))))^2 - (s/N)^2;
        %std = std + (log(abs(h(obs(j),i*fcst_step)/all_res(obs(j),i*fcst_step,i))))^2 - (s/N)^2;
    end
    kappa(i) = exp(sqrt(std/N));
end
ax = (1:fcst_runs)*dt*fcst_step;
figure(3)
plot(ax,K)
hold on
plot(ax,kappa)
hold on
plot([1,fcst_runs]*dt*fcst_step,[1,1],'--')
legend('K','kappa')
xlabel('Time (s)')
title('Accuracy of Optimal Interpolation (across grid)')

% plotting accuracy at the coast for each forecast run
K = zeros(1,fcst_runs);  % mean (0.8<K<1.2 is good)
kappa = zeros(1,fcst_runs);  % std (kappa<1.4 is good)
for i = 1:fcst_runs
    s = 0;
    std = 0;
    N = nt - i*fcst_step + 1;
    for j = i*fcst_step:nt
        s = s + log(abs(h(1,j)/all_res(1,j,i)));
    end
    K(i) = exp(s/N);
    for j = i*fcst_step:nt
        std = std + (log(abs(h(1,j)/all_res(1,j,i))))^2 - (s/N)^2;
    end
    kappa(i) = exp(sqrt(std/N));
end
figure(4)
plot(ax,K)
hold on
plot(ax,kappa)
hold on
plot([1,fcst_runs]*dt*fcst_step,[1,1],'--')
legend('K','kappa')
xlabel('Time (s)')
title('Accuracy of Optimal Interpolation (at coast)')
hold off

%% Analyze assimilation at stations
%{
assimilation = zeros(nx,assim_runs);
initial = h(:,1:assim_runs);
figure(5)
for i = 1:assim_runs
    assimilation(:,i) = all_res(:,i*assim_step,i);
    plot(x,assimilation(:,i),'r',x,initial(:,i),'b')
    drawnow
end
diff = zeros(length(obs),asssim_runs);
figure(6)
for i = 1:runs
    diff(:,i) = assimilation(obs,i) - initial(obs,i);
    plot(x_obs,diff(:,i))
    ylim([-0.5,0.5])
    drawnow
end
%}