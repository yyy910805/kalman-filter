function [u, g] = bc1linear_1d(xmin, xmax, tmin, tmax, dx, dt, c, opt)

% @param opt: 1 if want to see the graphs, 0 if not

% number of steps for x and t
nx = (xmax-xmin)/dx;
nt = (tmax-tmin)/dt + 1;
% create matrix u to store results
u = zeros(nt,nx);
x = xmin:dx:(xmax-dx);
t = tmin:dt:tmax;
% true initial condition
w = (xmax-xmin)/2;  % wavelength
u(1,:) = abs(sin((2*pi/w)*x));
% boundary condition
g = 1;
u(:,1) = g;

% set matrix A to compute iteratively the solution
e0 = zeros(nx,1);
e0(1) = 1;
M = eye(nx);
M(1,1) = -1;
for i = 2:nx
    M(i,i-1) = -1;
end
M = c*M/dx;

% 1st order Runge Kutta for time step
for j = 2:nt
    k1 = -M*u(j-1,:)' - 2/dx*e0*((e0')*u(j-1,:)' - g);
    u(j,:) = (u(j-1,:)' + k1*dt)';
    
    if opt == 1
        plot(x,u(j,:))
        %pause
        drawnow
    end
end

if opt == 1
    subplot(2,1,1),plot(x,u(1,:)),title('t=0');
    hold on;
    subplot(2,1,2),plot(x,u(nt,:)),title('t=1');
    xlabel('x');
    ylabel('u');
end

end