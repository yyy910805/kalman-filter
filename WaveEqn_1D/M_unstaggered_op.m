function [q_out,h_out] = M_unstaggered_op(adj,add,q_in,h_in,M,D,SAT,A,cdiss)

    n = size(q_in,1);

    % FORWARD
    if adj == 0
        if add == 0 
            q_out = 0.0;
            h_out = 0.0;
        end
                   
        q_out = -M.c.^2.*(D*h_in) + cdiss*A*q_in;
        h_out = -D*q_in + cdiss*A*h_in;
                
        % left boundary
        qhats = 0;
        hhats = h_in(1,1) - (q_in(1,1)-qhats)/M.c(1);      
        q_out(1,1) = q_out(1,1) - SAT(1)*(q_in(1,1)-qhats);
        h_out(1,1) = h_out(1,1) - SAT(1)*(h_in(1,1)-hhats);

        % right boundary
        qhate = 0.5*(q_in(n,1)+M.c(n)*h_in(n,1));
        hhate = qhate/M.c(n);
        q_out(n,1) = q_out(n,1) - SAT(2)*(q_in(n,1)-qhate);
        h_out(n,1) = h_out(n,1) - SAT(2)*(h_in(n,1)-hhate);        

    % ADJOINT
    else
        if add == 0 
            q_out = 0.0;
            h_out = 0.0;
        end
        q_out = -D'*h_in + cdiss*A'*q_in;
        h_out = -D'*(M.c.^2.*q_in) + cdiss*A'*h_in;

        % left boundary
        qhats = 0;
        hhats = h_in(1,1) - (q_in(1,1)-qhats)/M.c(1);          
        q_out(1,1) = q_out(1,1)-SAT(1)*q_in(1,1)-SAT(1)*h_in(1,1)/M.c(1);

        % right boundary
        qhate = 0.5*(q_in(n,1)+M.c(n)*h_in(n,1));
        hhate = qhate/M.c(n);
        q_out(n,1) = q_out(n,1) - SAT(2)*0.5*q_in(n,1) + SAT(2)/M.c(n)/2*h_in(n,1);
        h_out(n,1) = h_out(n,1) + SAT(2)*M.c(n)*0.5*q_in(n,1) - 0.5*SAT(2)*h_in(n,1);
    end
    
end