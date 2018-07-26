function X = CG(A,B)
tol = 1e-8;
m = size(B,1);
n = size(B,2);
X = zeros(m,n);
for i = 1:n
    x = zeros(m,1);
    r = B(:,i);
    p = r;
    n0 = norm(r,2);
    while norm(r,2)/n0 > tol
        alpha = r'*r/(p'*A*p);
        x = x + alpha*p;
        r_ = r - alpha*A*p;
        beta = r_'*r_/(r'*r);
        r = r_;
        p = r + beta*p;
    end
    X(:,i) = x;
end

end