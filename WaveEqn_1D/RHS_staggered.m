function [Dq, Dh]=RHS_staggered(q,h,M,D1_primal,D1_dual,SAT_primal,SAT_dual)

% Let q live on primal grid and h on dual
% D1_primal*h is approximation to partial h / partial x on primal grid

Dq = -M.c_primal.^2.*(D1_primal*h);
Dh = -D1_dual*q;

% forcing from seafloor motion
Dh = Dh + 0; % ADD FORCING HERE

% contribution from boundary conditions
L = M.L;

% left boundary

% solve for qhat and hhat
qhat = 0; % zero flux (= zero velocity)
hhat = h(1)-(q(1)-qhat)/M.c(1);

% add SAT terms (this will be the same for any BC)
Dq(1) = Dq(1)-SAT_primal(1)*(q(1)-qhat);
Dh(1) = Dh(1)-SAT_dual(1)*(h(1)-hhat);

% right boundary 

% solve for qhat and hhat
% (absorbing BC)
qhat = 0.5*(q(end)+M.c(L)*h(end));
hhat = qhat/M.c(L);

% add SAT terms (this will be the same for any BC)
Dq(end) = Dq(end)-SAT_primal(end)*(q(end)-qhat);
Dh(end) = Dh(end)-SAT_dual(end)*(h(end)-hhat);