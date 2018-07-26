% Waveform analysis using space-time plots

load('enkf40_result.mat')
load('optinterp40_result.mat')

% examine waveform at the 50th run
i = 50;
[nx,nt,runs] = size(all_res);
step = nt/(runs + 1);
opt_50 = all_res(:,:,50);
enkf_50 = all_res_enkf(:,:,50);
t_new = t(i*step:end);

figure
p1 = pcolor(x,t_new,h(:,i*step:end)');
set(p1,'LineStyle','none');
cmap();
xlabel('Distance from the Coast (km)');
ylabel('Time (s)');
title('Space-Time Plot of True Tsunami Wave Height (m)');
colorbar
caxis([-3 3])

figure
p2 = pcolor(x,t_new,opt_50(:,i*step:end)');
set(p2,'LineStyle','none');
cmap();
xlabel('Distance from the Coast (km)');
ylabel('Time (s)');
title('Space-Time Plot of Optimal Interpolation Tsunami Wave Height (m)');
colorbar
caxis([-3 3])

figure
p3 = pcolor(x,t_new,enkf_50(:,i*step:end)');
set(p3,'LineStyle','none');
cmap();
xlabel('Distance from the Coast (km)');
ylabel('Time (s)');
title('Space-Time Plot of EnKF Tsunami Wave Height (m)');
colorbar
caxis([-3 3])
