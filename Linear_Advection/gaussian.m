function C = gaussian(nx, dx, Ld)
%
% @param nx: number of grid points in space
% @param dx: spacing of grid points
% @param Ld: spatial decorrelation length
% @return C: correlation matrix
%

C = zeros(nx,nx);
for i = 1:nx
    for j = 1:nx
        C(i,j) = exp(-((j-i)*dx)^2/(2*Ld^2));
    end
end

end