%{
runs the enkf_cascadia results using optimal interpolation
analyzes the result in terms of:
 - first time the tsunami wave hits the coast
 - max height of tsunami and its arrival time at the coast
 - accuracy of prediction
 - waveform at the fth forecast
%}

%% get outputs
order = 4;
assim_step = 10; % assimilate observations every x time steps
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
end
ax = (1:fcst_runs)*dt*fcst_step;
figure(3)
plot(ax,K)
xlabel('Time (s)')
title('Accuracy of Optimal Interpolation (across grid)')
ylim([0 1])
% Convert y-axis values to percentage values
a = cellstr(num2str(get(gca,'ytick')'*100)); 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)

% plotting accuracy at the coast for each forecast run
K = zeros(1,fcst_runs);  % mean (0.8<K<1.2 is good)
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
end
figure(4)
plot(ax,K)
xlabel('Time (s)')
title('Accuracy of Optimal Interpolation (at coast)')
ylim([0 1])
a = cellstr(num2str(get(gca,'ytick')'*100)); 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)

%% plot waveform at the fth forecast
figure(8)
f = 10;
[nx,nt,runs] = size(all_res);
opt = all_res(:,:,f);
p = pcolor(x,t,opt');
set(p,'LineStyle','none');
cmap();
xlabel('Distance from the Coast (km)');
ylabel('Time (s)');
title('Space-Time Plot of Optimal Interpolation Tsunami Wave Height (m)')
colorbar
caxis([-3 3])
hold on
plot([x(1),x(end)],[f*dt*fcst_step,f*dt*fcst_step],'k--')
hold off