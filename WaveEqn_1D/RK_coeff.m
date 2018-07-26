function coeff = RK_coeff(A,B,dt)

    alpha = zeros(5,5);
    beta = zeros(5,5);
    coeff = zeros(1,6);
    
    % substep 1
        % rates
        beta(1,2)=1.0;
        
        % flux and height
        alpha(1,1) = 1.0;
        alpha(1,2) = dt*B(1);
    
    % substep 2
        % rates
        beta(2,2) = A(2)*beta(1,2)+alpha(1,1);
        beta(2,3) = alpha(1,2);

        % flux and height
        alpha(2,1) = alpha(1,1);
        alpha(2,2) = alpha(1,2)+dt*B(2)*beta(2,2);
        alpha(2,3) = dt*B(2)*beta(2,3);
        
    % substep 3
        % rates
        beta(3,2) = A(3)*beta(2,2)+alpha(2,1);
        beta(3,3) = A(3)*beta(2,3)+alpha(2,2);
        beta(3,4) = alpha(2,3);

        % flux and height
        alpha(3,1) = alpha(2,1);
        alpha(3,2) = alpha(2,2)+dt*B(3)*beta(3,2);
        alpha(3,3) = alpha(2,3)+dt*B(3)*beta(3,3);
        alpha(3,4) = dt*B(3)*beta(3,4);
        
    % substep 4
        % rates
        beta(4,2) = A(4)*beta(3,2)+alpha(3,1);
        beta(4,3) = A(4)*beta(3,3)+alpha(3,2);
        beta(4,4) = A(4)*beta(3,4)+alpha(3,3);
        beta(4,5) = alpha(3,4);
        
        % flux and height
        alpha(4,1) = alpha(3,1);
        alpha(4,2) = alpha(3,2)+dt*B(4)*beta(4,2);
        alpha(4,3) = alpha(3,3)+dt*B(4)*beta(4,3);
        alpha(4,4) = alpha(3,4)+dt*B(4)*beta(4,4);
        alpha(4,5) = dt*B(4)*beta(4,5);
        
    % substep 5
        % rates
        beta(5,2) = A(5)*beta(4,2)+alpha(4,1);
        beta(5,3) = A(5)*beta(4,3)+alpha(4,2);
        beta(5,4) = A(5)*beta(4,4)+alpha(4,3);
        beta(5,5) = A(5)*beta(4,5)+alpha(4,4);
        beta(5,6) = alpha(4,5);
        
        % flux and height
        alpha(5,1) = alpha(4,1);
        alpha(5,2) = alpha(4,2)+dt*B(5)*beta(5,2);
        alpha(5,3) = alpha(4,3)+dt*B(5)*beta(5,3);
        alpha(5,4) = alpha(4,4)+dt*B(5)*beta(5,4);
        alpha(5,5) = alpha(4,5)+dt*B(5)*beta(5,5);
        alpha(5,6) = dt*B(5)*beta(5,6);
        
    coeff(1,:) = alpha(5,:);
    
end