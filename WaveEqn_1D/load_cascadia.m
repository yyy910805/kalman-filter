% for Cascadia
% loads the initial wave height, x grid and seafloor bathymetry from datafile,
% and interpolates data based on parameters file (params)

% outputs the bathymetry (H), initial wave height (h), initial flux (q),
% and x grid (x)

function [H,h,q,x] = load_cascadia(datafile, params)

% read paramaters from file
[xmin, xmax, ~, ~, dx, ~, ~, ~, ~, ~] = readParams(params);

% import water depth data
irreg = load(datafile);

% flip the axis first to make coast on the left
x_irreg = -flipud(irreg.x_SF);
x_irreg = x_irreg - x_irreg(1);
y_irreg = flipud(irreg.y_SF); % irregularly sampled y-position [km]

% initial wave height at the moment of rupture on irregular grid
h_irreg = flipud(irreg.uy_os);

% get rid of overlapping locations in the x-position
h_irreg(diff(x_irreg)==0) = [];
y_irreg(diff(x_irreg)==0) = [];
x_irreg(diff(x_irreg)==0) = [];

% interpolation on regular grid
x = xmin:dx:xmax;
nx = length(x);
H = -interp1(x_irreg,y_irreg,x'); % sea floor depth [km]
h_reg = interp1(x_irreg,h_irreg,x'); % wave height

% compute weighting to remove the body waves in Gabe's simulation
weight = zeros(nx,1);
for i = 1:nx
    easting = xmin + dx*(i-1);
    if (easting > 300)
        weight(i,1) = sin(3.14/2*(xmax - easting)/(xmax-300))^10; 
    else
        weight(i,1) = 1;
    end
end
h = h_reg.*weight;
q = zeros(length(h),1);

% plots initial water height and flux, and seafloor bathymetry
plot(x,h)
hold on
plot(x,-H)
hold on
plot(x,q)
hold off