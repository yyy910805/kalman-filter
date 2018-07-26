function [Dq Dh]=RHS(t,q,h,M,Dp,Dm,SAT)

% contribution from PDE, no boundary conditions
% D*h is approximation to partial h / partial x
% Dp = Dm = D for unstaggered grid

Dq = -M.c.^2.*(Dm*h);
Dh = -Dp*q;

% forcing from seafloor motion
Dh = Dh + 0; % ADD FORCING HERE

% contribution from boundary conditions

% left boundary

% solve for qhat and hhat
qhat = 0; % zero flux (= zero velocity)
hhat = h(1)-(q(1)-qhat)/M.c(1);

% add SAT terms (this will be the same for any BC)
Dq(1) = Dq(1)-SAT(1)*(q(1)-qhat);
Dh(1) = Dh(1)-SAT(1)*(h(1)-hhat);

% right boundary 

% solve for qhat and hhat
% (absorbing BC)
qhat = 0.5*(q(end)+M.c(end)*h(end));
hhat = qhat/M.c(end);

% add SAT terms (this will be the same for any BC)
Dq(end) = Dq(end)-SAT(end)*(q(end)-qhat);
Dh(end) = Dh(end)-SAT(end)*(h(end)-hhat);