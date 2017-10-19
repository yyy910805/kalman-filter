function [h,io] = obspat ( nobs, jdim )

%OBSPAT    Set obs pattern depending on value of nobs 
%
%   Usage:  [h,y,mobs] = obspat ( nobs, jdim )
%
%   where   jdim  - size of domain
%           nobs  - number of grid points to consider for obs pattern
%
% On ouput:  h  - observation matrix
%            io - array of points where obs occur 
%             

if nobs<0    , io = 1:abs(nobs); end
if nobs>0    , io = 1:nobs:jdim; end
if nobs==jdim, io = 1:jdim     ; end
p = length(io);
h = zeros(p,jdim);
h(:,io(:)) = eye(p,p);