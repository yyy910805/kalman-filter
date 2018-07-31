% get PDF of wave height prediction for the nth forecast
% e.g. fcst_step = 20, n = 20 is after 400 time steps, i.e. 10 min after the
% tsunami generation
% assumes we've finished running 'run_enkf_casc.m'
n = 20;
elapsed_t = n*fcst_step;
fcst_pdf = zeros(nt-elapsed_t,nx,N);
v = assim;
for j = elapsed_t+1:nt
    v = T*v;
    fcst_pdf(j-elapsed_t,:,:) = v(1+nx:2*nx,:);
end

max_fcst = zeros(nt - elapsed_t,nx);
for i = 1:nt-elapsed_t
    for j = 1:nx
        max_fcst(i,j) = max(fcst_pdf(i,j,:));
    end
end

min_fcst = zeros(nt - elapsed_t,nx);
for i = 1:nt - elapsed_t
    for j = 1:nx
        min_fcst(i,j) = min(fcst_pdf(i,j,:));
    end
end

mean_fcst = zeros(nt - elapsed_t,nx);
for i = 1:nt - elapsed_t
    for j = 1:nx
        mean_fcst(i,j) = mean(fcst_pdf(i,j,:));
    end
end

% plot histogram of forecast when mean is max at coast
[val,idx] = max(mean_fcst(:,1));
[~, ~, N] = size(fcst_pdf(idx,1,:));
V = reshape(fcst_pdf(idx,1,:),N,1);
g = histfit(V,25,'kernel');
g(1).FaceColor = [0.3 0.75 0.93];
xlabel('Wave Height (m)')
ylabel('Frequency')
title('Distribution of Max Wave Height Forecast at Coast (10min)')

% plot mean, max, mean propagation of 10th forecast
for i = 1:nt - elapsed_t
    plot(x,mean_fcst(i,:),x,max_fcst(i,:),x,min_fcst(i,:))
    drawnow
end