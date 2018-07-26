function X = CG2(A,B)
tol = 1e-8;
maxiter = 20;
m = size(B,1);
n = size(B,2);
X = zeros(m,n);
L = ichol(A);
for i = 1:n
    [x, flag] = pcg(A,B(:,i),tol,maxiter,L,L');
    X(:,i) = x;
end

end