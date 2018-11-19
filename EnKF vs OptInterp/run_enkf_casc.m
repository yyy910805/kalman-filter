%{
runs the enkf_cascadia results using ensemble Kalman filter
analyzes the result in terms of:
 - first time the tsunami wave hits the coast
 - max height of tsunami and its arrival time at the coast
 - accuracy of prediction
%}

%% get outputs
order = 4;
assim_step = 2; % assimilate observations every x time steps
fcst_step = 20; % forecast every x time steps
out = enkf_cascadia('params.txt',order,assim_step,fcst_step);

h = out.real;
all_res = out.all_res;
obs = out.obs;
stdev = out.std;
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

%% plot prediction over forecasts
arrival_run = ceil(real_fcst(2)/dt/fcst_step) + 1;
ax = (1:arrival_run)*dt*fcst_step;
figure(1)
plot(ax, pred(1:arrival_run,1))
ylim([0 3])
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
%kappa = zeros(1,fcst_runs);  % std (kappa<1.4 is good)
N = nx;  % over the entire grid
%N = length(obs);  % at observation stations

for i = 1:fcst_runs
    s = 0;
    std = 0;
    for j = 1:N
        s = s + log(abs(h(j,i*fcst_step)/all_res(j,i*fcst_step,i)));
        % s = s + log(abs(h(obs(j),i*fcst_step)/all_res(obs(j),i*fcst_step,i)));
    end
    if exp(s/N) < 1
        K(i) = exp(s/N);
    else
        K(i) = 1/exp(s/N);
    end
    %{
    for j = 1:N
        std = std + (log(abs(h(j,i*fcst_step)/all_res(j,i*fcst_step,i))))^2 - (s/N)^2;
        %std = std + (log(abs(h(obs(j),i*fcst_step)/all_res(obs(j),i*fcst_step,i))))^2 - (s/N)^2;
    end
    kappa(i) = exp(sqrt(std/N));
    %}
end
ax = (1:fcst_runs)*dt*fcst_step;
figure(3)
plot(ax,K)
%{
hold on
plot(ax,kappa)
hold on
plot([1,fcst_runs]*dt*fcst_step,[1,1],'--')
legend('K','kappa')
%}
xlabel('Time (s)')
title('Accuracy of Ensemble Kalman Filter (across grid)')
% Convert y-axis values to percentage values
a = cellstr(num2str(get(gca,'ytick')'*100)); 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)

% plotting accuracy at the coast for each forecast run
K = zeros(1,fcst_runs);  % mean (0.8<K<1.2 is good)
%kappa = zeros(1,fcst_runs);  % std (kappa<1.4 is good)
for i = 1:fcst_runs
    s = 0;
    std = 0;
    N = nt - i*fcst_step + 1;
    for j = i*fcst_step:nt
        s = s + log(abs(h(1,j)/all_res(1,j,i)));
    end
    if exp(s/N) < 1
        K(i) = exp(s/N);
    else
        K(i) = 1/exp(s/N);
    end
    %{
    for j = i*fcst_step:nt
        std = std + (log(abs(h(1,j)/all_res(1,j,i))))^2 - (s/N)^2;
    end
    kappa(i) = exp(sqrt(std/N));
    %}
end
figure(4)
plot(ax,K)
%{
hold on
plot(ax,kappa)
hold on
plot([1,fcst_runs]*dt*fcst_step,[1,1],'--')
legend('K','kappa')
%}
xlabel('Time (s)')
title('Accuracy of Ensemble Kalman Filter (at coast)')
a = cellstr(num2str(get(gca,'ytick')'*100)); 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)

%% plot standard deviation evolution
figure(7)
t1 = round(60/dt/assim_step);
t2 = round(300/dt/assim_step);
t3 = round(900/dt/assim_step);
t4 = round(1800/dt/assim_step);
t5 = round(2700/dt/assim_step);
% convert runs to time, round to nearest 10
ts = roundn([t1 t2 t3 t4 t5]*assim_step*dt,1);
plot(x,stdev(t1,:))
hold on
plot(x,stdev(t2,:))
hold on
plot(x,stdev(t3,:))
hold on
plot(x,stdev(t4,:))
hold on
plot(x,stdev(t5,:))
legend([num2str(ts(1)) 's'],[num2str(ts(2)) 's'],[num2str(ts(3)) 's'],...
       [num2str(ts(4)) 's'],[num2str(ts(5)) 's'])
xlabel('Distance from Coast (km)')
ylabel('Standard Deviation')
title('EnKF Evolution of Standard Deviation across the Grid')
hold off

%% plot waveform at the fth forecast
figure(8)
f = 20;
[nx,nt,runs] = size(all_res);
step = nt/(runs + 1);
enkf = all_res(:,:,f);
%t_new = t(f*step:end);
%p = pcolor(x,t_new,enkf(:,f*step:end)');
p = pcolor(x,t,enkf');
set(p,'LineStyle','none');
cmap();
xlabel('Distance from the Coast (km)');
ylabel('Time (s)');
title('Space-Time Plot of EnKF Tsunami Wave Height (m)');
colorbar
caxis([-3 3])
hold on
plot([x(1),x(end)],[f*dt*fcst_step,f*dt*fcst_step],'k--')
hold off