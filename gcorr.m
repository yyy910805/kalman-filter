function [ Q, D ] = gcorr ( func, Lx, Ld, J, opt, Lc )

%GCORR  Construct a correlation matrix on a periodic domain.
%
%       Usage:  [ Q, D ] = gcorr ( func, Lx, Ld, J, opt [, Lc] )
%
%       Example:  [ Q, D ] = gcorrnp ( 'gauss', 1, 0.2, 32, 1 );
%         or      [ Q, D ] = gcorrnp ( 'gauss', 1, 0.2, 32, 1, 0.6);
%
%       Inputs:   func-  representation function, e.g., gauss, foar
%                 Lx  -  size of domain 
%                 Ld  -  (de)correlation length 
%                 J   -  number of grid points 
%                 opt -  0=no plots produced
%                        1=contour plot of Q
%                 [Lc] - optional argument referring to
%                        cutoff length; in this case, Q
%                        will not be positive semi-definite
%
%       Outputs:  Q  - correlation matrix on the grid
%                 D  - eigenvalues of Q (should all be >=0, unless
%                      Lc has been specified) 
%
%    Some of the functions available can be used by 
%    commenting/uncommenting lines below.
%
%  25feb97  R. Todling   Initial code.

Q   =  zeros(J);
Ix = (-J/2+1:J/2);
Iy = (-J/2+1:J/2)';
[I1,I2] = meshgrid(Ix,Iy);
Idiff   = (I1-I2);
Dist    = mod(Idiff,J);
Q       = Lx*sin(pi*Dist/J);      
%Q       = (abs(Q)<Lc).*Q;      %   Cutoff by force
                                %   resulting Q won't be positive
%  
%  Pick representing correlation function
%    include cutoff if that's the case
%
if nargin < 6
    [ Q ] = feval ( func, Q, Ld );
elseif nargin == 6
    [ Q ] = feval ( func, Q, Ld, Lc );
else
  help gcorr
end
%
%  Eigenvalues and plotting follows here
%
[V,D] = eig(Q);
%
%  Check for positiveness
%
if(real(diag(D))<-1.e10)
   echo on,
   pause % Matrix is not positive definite ...
   min(D)
end   
   
if ( opt == 0 )
    return
end

if ( opt == 1 )
  figure(3), clf; figure(2), clf;
  figure(3);contour(Q)
  ind1 = J/2; ind2 = J;
  figure(2);subplot(2,1,1),plot(Q(:,ind1)),hold on, ...
     subplot(2,1,2),plot(Q(:,ind2)),hold off
  D = sort(real(diag(D)));
  D = D(J:-1:1);
  figure(1);plot(D),title('Eigenvalues of Q');
end

