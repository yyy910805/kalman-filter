function u = bc4linear_1d(xmin, xmax, tmin, tmax, dx, dt, c, opt)

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
u(1,:) = sin((2*pi/w)*x);
% boundary condition
u(:,1) = exp(-t);

% set matrix A to compute iteratively the solution
e0 = zeros(nx);
e0(1) = 1;
A = eye(nx);
A(1,1) = -1;
for i = 2:nx
    A(i,i-1) = -1;
end
A = c*A/dx;

% 4th order Runge Kutta for time step
for j = 2:nt
    k1 = -A*u(j-1,:)' - 2/dx*e0*(e0'*u(j-1,:)' - exp(-(j-1)*dt));
    y1 = u(j-1,:)' + k1*dt/2;
    k2 = -A*y1 - 2/dx*e0*(e0'*y1 - exp(-(j-1)*dt));
    y2 = u(j-1,:)' + k2*dt/2;
    k3 = -A*y2 - 2/dx*e0*(e0'*y2 - exp(-(j-1)*dt));
    y3 = u(j-1,:)' + k3*dt;
    k4 = -A*y3 - 2/dx*e0*(e0'*y3 - exp(-(j-1)*dt));
    u(j,:) = (u(j-1,:)' + (k1 + 2*k2 + 2*k3 + k4)*dt/6)';
    
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