%{
1. Runs the enkf_cascadia results using ensemble Kalman filter
and optimal interpolation, does comparison
2. Analyzes the result in terms of:
 - first time the tsunami wave hits the coast
 - max height of tsunami and its arrival time at the coast
 - accuracy of prediction
 - waveform at any forecast
%}

%% get outputs
addpath 'kalman-filter/EnKF vs OptInterp'
order = 4;
assim_step = 10; % assimilate observations every x time steps
fcst_step = 20; % forecast every x time steps
out1 = enkf_cascadia('params.txt',order,assim_step,fcst_step);
out2 = optinterp_cascadia('params.txt',order,assim_step,fcst_step);

all_res1 = out1.all_res;
all_res2 = out2.all_res;

h = out1.real;
obs = out1.obs;
stdev = out1.std;
x = out1.x;
dx = out1.dx;
t = out1.t;
dt = out1.dt;
assim_runs = out1.assim_runs;
fcst_runs = out1.fcst_runs;

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
pred1 = zeros(fcst_runs,2);
pred2 = zeros(fcst_runs,2);
for i = 1:fcst_runs
    [val,idx] = max(all_res1(1,:,i));
    pred1(i,1) = val;
    pred1(i,2) = idx;
    [val,idx] = max(all_res2(1,:,i));
    pred2(i,1) = val;
    pred2(i,2) = idx;
end

%% plot prediction over runs
% plot only until the actual max wave hits the coast
arrival_run = ceil(real_fcst(2)/dt/fcst_step) + 1;
ax = (1:arrival_run)*dt*fcst_step;
figure(1)
plot(ax, pred1(1:arrival_run,1), ax, pred2(1:arrival_run,1))
ylim([0 3.5])
hold on
plot([1,arrival_run]*dt*fcst_step,[real_fcst(1),real_fcst(1)],'--')
legend('EnKF','OptInterp','True Max Wave Height','Location','Best')
xlabel('Time (s)')
ylabel('Wave Height (m)')
title('EnKF vs. OptInterp Forecasted Max Wave Height')
hold off

figure(2)
plot(ax, pred1(1:arrival_run,2), ax, pred2(1:arrival_run,2))
ylim([0 1500])
hold on
plot([1,arrival_run]*dt*fcst_step,[real_fcst(2),real_fcst(2)],'--')
legend('EnKF','OptInterp','True Max Wave Height','Location','Best')
xlabel('Time (s)')
ylabel('Max Wave Height Arrival Time (s)')
title('EnKF vs. OptInterp Forecasted Max Wave Height Arrival Time')
hold off

%% plot accuracy indicators over runs
% true vs assimilated waveform before each forecast
% over the entire grid, or at observation stations
K1 = zeros(1,fcst_runs);
K2 = zeros(1,fcst_runs);
N = nx;  % over the entire grid
%N = length(obs);  % at observation stations

for i = 1:fcst_runs
    s1 = 0;
    s2 = 0;
    for j = 1:N
        s1 = s1 + log(abs(h(j,i*fcst_step)/all_res1(j,i*fcst_step,i)));
        s2 = s2 + log(abs(h(j,i*fcst_step)/all_res2(j,i*fcst_step,i)));
        %s1 = s1 + log(abs(h(obs(j),i*fcst_step)/all_res1(obs(j),i*fcst_step,i)));
        %s2 = s2 + log(abs(h(obs(j),i*fcst_step)/all_res2(obs(j),i*fcst_step,i)));
    end
    
    if exp(s1/N) < 1
        K1(i) = exp(s1/N);
    else
        K1(i) = 1/exp(s1/N);
    end
    
    if exp(s2/N) < 1
        K2(i) = exp(s2/N);
    else
        K2(i) = 1/exp(s2/N);
    end
end
ax = (1:fcst_runs)*dt*fcst_step;
figure(3)
plot(ax, K1, '-+','MarkerSize',3)
hold on
plot(ax, K2, '-*','MarkerSize',4)
xlabel('Time (s)')
title('Accuracy of EnKF vs. OptInterp across the Grid')
legend('EnKF','OptInterp')
ylim([0 1])
% Convert y-axis values to percentage values
a = cellstr(num2str(get(gca,'ytick')'*100)); 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)

% plotting accuracy at the coast for each forecast run
C1 = zeros(1,fcst_runs);
C2 = zeros(1,fcst_runs);
for i = 1:fcst_runs
    s1 = 0;
    s2 = 0;
    N = nt - i*fcst_step + 1;
    for j = i*fcst_step:nt
        s1 = s1 + log(abs(h(1,j)/all_res1(1,j,i)));
        s2 = s2 + log(abs(h(1,j)/all_res2(1,j,i)));
    end
    
    if exp(s1/N) < 1
        K1(i) = exp(s1/N);
    else
        K1(i) = 1/exp(s1/N);
    end
    
    if exp(s2/N) < 1
        K2(i) = exp(s2/N);
    else
        K2(i) = 1/exp(s2/N);
    end
end
figure(4)
plot(ax, K1, '-+','MarkerSize',3)
hold on
plot(ax, K2, '-*','MarkerSize',4)
xlabel('Time (s)')
title('Accuracy of EnKF vs. OptInterp at the Coast')
legend('EnKF','OptInterp')
ylim([0 1])
a = cellstr(num2str(get(gca,'ytick')'*100)); 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)

%% plot waveform at the fth forecast
f = 20;
[nx,nt,runs] = size(all_res1);
step = nt/(runs + 1);
enkf = all_res1(:,:,f);
opt = all_res2(:,:,f);

figure(5)
p1 = pcolor(x,t,enkf');
set(p1,'LineStyle','none');
cmap();
xlabel('Distance from the Coast (km)');
ylabel('Time (s)');
title('Space-Time Plot of EnKF Tsunami Wave Height Forecast at 600s (m)')
colorbar
caxis([-3 3])
hold on
plot([x(1),x(end)],[f*dt*fcst_step,f*dt*fcst_step],'k--')
hold off

figure(6)
p2 = pcolor(x,t,opt');
set(p2,'LineStyle','none');
cmap();
xlabel('Distance from the Coast (km)');
ylabel('Time (s)');
title('Space-Time Plot of OptInterp Tsunami Wave Height Forecast at 600s (m)')
colorbar
caxis([-3 3])
hold on
plot([x(1),x(end)],[f*dt*fcst_step,f*dt*fcst_step],'k--')
hold off

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