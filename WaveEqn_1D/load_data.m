function [d,h] = load_data(depth_file, height_file, xmin, dx, xmax)

% loads the water depth and initial tsunami height data

% import water depth data
irreg = load(depth_file);
x_irreg = irreg.x_es; % irregularly sampled x-position [km]
y_irreg = irreg.y_es; % irregularly sampled y-position [km]

% load initial wave height at the moment of rupture on irregular grid
h_irreg = load(height_file);  
h_irreg_load = h_irreg.tsunami;

% get rid of overlapping locations in the x-position
h_irreg_load(diff(x_irreg)==0) = [];
y_irreg(diff(x_irreg)==0) = [];
x_irreg(diff(x_irreg)==0) = [];

% interpolation on regular grid
x1 = xmin:dx:xmax;
x2 = [xmin,(xmin+dx/2):dx:(xmax-dx/2), xmax];
nx = length(x1);
d = -interp1(x_irreg,y_irreg,x2'); % sea floor depth [km]
h_reg = interp1(x_irreg,h_irreg_load,x1'); % wave height

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