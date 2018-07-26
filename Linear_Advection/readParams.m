function [xmin, xmax, tmin, tmax, dx, dt, c, Ld, xsd, esd, freq] = readParams(filename)
%{
@return xmin: minimum value of x
@return xmax: maximum value of x in the period
@return tmin: minimum value of time
@return tmax: maximum time over which the assimilation is done
@return dx: step size of x
@return dt: step size of t
@return c: speed of propagation
@return Ld: spatial decorrelation length
@return xsd: initial estimate error standard deviation
@return esd: measurement error standard deviation
@return freq: frequency of observation/measurement
%}

fid = fopen(filename);
params = textscan(fid, '%[^= ]%*[= ]%f', 'CommentStyle', '%');
fclose(fid);
xmin = params{2}(strcmp(params{1}, 'xmin'));
xmax = params{2}(strcmp(params{1}, 'xmax'));
tmin = params{2}(strcmp(params{1}, 'tmin'));
tmax = params{2}(strcmp(params{1}, 'tmax'));
dx = params{2}(strcmp(params{1}, 'dx'));
dt = params{2}(strcmp(params{1}, 'dt'));
c = params{2}(strcmp(params{1}, 'c'));
Ld = params{2}(strcmp(params{1}, 'Ld'));
xsd = params{2}(strcmp(params{1}, 'xsd'));
esd = params{2}(strcmp(params{1}, 'esd'));
freq = params{2}(strcmp(params{1}, 'freq'));

end