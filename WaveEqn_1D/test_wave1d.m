% testing wave1d.m wave propagation for yuyun_cascadia.mat (which contains
% initial tsunami height, x grid and seafloor bathymetry)

[H,h,q,x] = load_cascadia('yuyun_cascadia.mat', 'params.txt');

out1 = wave1d_assimilate(q,h,H,x,4,1,'unstaggered');

%out2 = wave1d_assimilate(q,h,H,x,4,1,'staggered'); % does not work yet