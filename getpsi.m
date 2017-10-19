function psi = getpsi ( func, jdim )

%GETPSI  builds dynamics matrix by appying function func to columns
%        of the identity matrix.
%
%R. Todling

psi = eye(jdim);
for j = 1:jdim
  psi(:,j) = feval ( func, psi(:,j), 0.0 );
end
