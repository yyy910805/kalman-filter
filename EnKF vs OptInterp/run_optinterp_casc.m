%{
runs the enkf_cascadia results using optimal interpolation
analyzes the result in terms of:
 - first time the tsunami wave hits the coast
 - max height of tsunami and its arrival time at the coast
 - accuracy of prediction
%}

addpath 'EnKF vs OptInterp'
order = 4;
step = 10; %(assimilate observations every 10 time steps)
out = optinterp_cascadia('params.txt',order,step);

h = out.real;
x = out.x;
t = out.t;
obs = out.obs;
all_res = out.all_res;

% gives the dimension of all_res
[nx,nt,runs] = size(all_res);

% real prediction
real = zeros(1,3);
idx = find(h(1,:)>1,1);
real(1) = t(idx);
[val,idx] = max(h(1,:));
real(2) = val;
real(3) = t(idx);

% we only need the runs before the max wave height reaches the coast
% for plotting of max wave height
%n = floor(idx/10);
n = runs;

% array to store prediction results for the coast
% col 1 = first time wave height > 1.5m
% col 2 = max wave height
% col 3 = arrival time of max wave height

%
pred = zeros(n,3);
for i = 1:n
    if isempty(find(all_res(1,:,i)>1,1)) == 1
        idx = 0;
    else
       idx = find(all_res(1,:,i)>1,1);
    end
    pred(i,1) = idx;
    [val,idx] = max(all_res(1,:,i));
    pred(i,2) = val;
    pred(i,3) = idx;
end
%}

% plot prediction over runs
%{
ax = 1:n;
figure(1)
plot(ax, pred(:,1))
hold on
plot([1,n],[real(1),real(1)],'--')
title('First time wave height > 1m')
hold off

figure(2)
plot(ax, pred(:,2))
hold on
plot([1,n],[real(2),real(2)],'--')
title('Max wave height')
hold off

figure(3)
plot(ax, pred(:,3))
hold on
plot([1,n],[real(3),real(3)],'--')
title('Arrival time of max wave height (s)')
hold off
%}

% plot accuracy indicators over runs
% true vs assimilated waveform after each assimilation (as the new initial
% waveform)
% over the entire grid, or at observation stations
K = zeros(1,runs);  % mean (0.8<K<1.2 is good)
kappa = zeros(1,runs);  % std (kappa<1.4 is good)
N = length(obs);
for i = 1:runs
    s = 0;
    std = 0;
    for j = 1:N
        % s = s + log(abs(h(j,i*step)/all_res(j,i*step,i)));
        s = s + log(abs(h(obs(j),i*step)/all_res(obs(j),i*step,i)));
    end
    K(i) = exp(s/N);
    for j = 1:N
        % std = std + (log(abs(h(j,i*step)/all_res(j,i*step,i))))^2 - (s/N)^2;
        std = std + (log(abs(h(obs(j),i*step)/all_res(obs(j),i*step,i))))^2 - ...
              (s/N)^2;
    end
    kappa(i) = exp(sqrt(std/N));
end
ax = 1:runs;
figure
plot(ax,K)
hold on
plot(ax,kappa)
hold on
plot([1,runs],[1,1],'--')
legend('K','kappa')
title('Accuracy of Optimal Interpolation (at stations)')

% only plotting accuracy at the coast
% averaging over all time steps at each run (assimilation and propagation)
K = zeros(1,runs);  % mean (0.8<K<1.2 is good)
kappa = zeros(1,runs);  % std (kappa<1.4 is good)
for i = 1:runs
    s = 0;
    std = 0;
    N = nt - i*step + 1;
    for j = i*step:nt
        s = s + log(abs(h(1,j)/all_res(1,j,i)));
    end
    K(i) = exp(s/N);
    for j = i*step:nt
        std = std + (log(abs(h(1,j)/all_res(1,j,i))))^2;
    end
    kappa(i) = exp(sqrt(std/N - (s/N)^2));
end
ax = 1:runs;
figure
plot(ax,K)
hold on
plot(ax,kappa)
hold on
plot([1,runs],[1,1],'--')
legend('K','kappa')
title('Accuracy of Optimal Interpolation (at coast)')