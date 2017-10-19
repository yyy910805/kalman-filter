function [ Q ] = gauss ( Q, Ld, Lc ) 

%GAUSS  Gaussian function (cutoff allowed)
% 
%       Usage: [ Q ] = gauss ( Q, Ld, [Lc] )
%
%       Inputs:  Q  -  a grid over which to compute function
%                Ld -  decorrelation length of gaussian function
%               [Lc]-  optional cutoff value of correlations
% 
%  25feb97  R. Todling   Initial code.

if nargin < 3
    Q   = exp( -0.5 * (Q./Ld).^2 );                     %  Gaussian
elseif nargin == 3
    Q  = (abs(Q)<Lc).*exp( -0.5 * (Q./Ld).^2 );         %  Gaussian w/ cutoff
else
   help gcorr
end

