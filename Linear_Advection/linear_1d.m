function u = linear_1d(xmin, xmax, tmin, tmax, dx, dt, c, opt)
% number of steps for x and t
nx = (xmax-xmin)/dx;
nt = (tmax-tmin)/dt + 1;
% create matrix u to store results
u = zeros(nt,nx);
x = xmin:dx:(xmax-dx);
t = tmin:dt:tmax;
% true initial condition
w = (xmax-xmin)/2;  % wavelength
u(1,:) = sin((2*pi/w)*x);
% set matrix A to compute iteratively the solution
A = diag((1-c*dt/dx)*ones(1,nx));
A(1,nx) = c*dt/dx;
for i = 2:nx
    A(i,i-1) = c*dt/dx;
end
% iterate
for n = 2:nt
    v = A*u(n-1,:)';
    u(n,:) = v;
end

if opt == 1
    subplot(2,1,1),plot(x,u(1,:)),title('t=0');
    hold on;
    subplot(2,1,2),plot(x,u(nt,:)),title('t=1');
    xlabel('x');
    ylabel('u');
end
end