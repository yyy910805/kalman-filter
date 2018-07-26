function x = CG1(A,b)
tol = 10e-8;
x = zeros(length(b),1);
r = b;
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

end