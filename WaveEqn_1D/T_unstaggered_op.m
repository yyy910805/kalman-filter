function [q_out,h_out] = T_unstaggered_op(adj,add,q_in,h_in,alpha,M,D,SAT,A,cdiss)

    % FORWARD
    if adj == 0
        if add == 0 
            q_out = 0.0;
            h_out = 0.0;
        end
        
        [q_temp1,h_temp1] = M_unstaggered_op(0,0,q_in,h_in,M,D,SAT,A,cdiss);
        [q_temp2,h_temp2] = M_unstaggered_op(0,0,q_temp1,h_temp1,M,D,SAT,A,cdiss);
        [q_temp3,h_temp3] = M_unstaggered_op(0,0,q_temp2,h_temp2,M,D,SAT,A,cdiss);        
        [q_temp4,h_temp4] = M_unstaggered_op(0,0,q_temp3,h_temp3,M,D,SAT,A,cdiss);
        [q_temp5,h_temp5] = M_unstaggered_op(0,0,q_temp4,h_temp4,M,D,SAT,A,cdiss);        
        
        q_out = alpha(1,1)*q_in + alpha(1,2)*q_temp1 + alpha(1,3)*q_temp2 + ...
        alpha(1,4)*q_temp3 + alpha(1,5)*q_temp4 + alpha(1,6)*q_temp5;

        h_out = alpha(1,1)*h_in + alpha(1,2)*h_temp1 + alpha(1,3)*h_temp2 + ...
        alpha(1,4)*h_temp3 + alpha(1,5)*h_temp4 + alpha(1,6)*h_temp5; 

    % ADJOINT
    else
        if add == 0
            q_out = 0.0;
            h_out = 0.0;
        end
        
        [q_temp1,h_temp1] = M_unstaggered_op(1,0,q_in,h_in,M,D,SAT,A,cdiss);
        [q_temp2,h_temp2] = M_unstaggered_op(1,0,q_temp1,h_temp1,M,D,SAT,A,cdiss);
        [q_temp3,h_temp3] = M_unstaggered_op(1,0,q_temp2,h_temp2,M,D,SAT,A,cdiss);        
        [q_temp4,h_temp4] = M_unstaggered_op(1,0,q_temp3,h_temp3,M,D,SAT,A,cdiss);
        [q_temp5,h_temp5] = M_unstaggered_op(1,0,q_temp4,h_temp4,M,D,SAT,A,cdiss);        
        
        q_out = alpha(1,1)*q_in + alpha(1,2)*q_temp1 + alpha(1,3)*q_temp2 + ...
        alpha(1,4)*q_temp3 + alpha(1,5)*q_temp4 + alpha(1,6)*q_temp5;

        h_out = alpha(1,1)*h_in + alpha(1,2)*h_temp1 + alpha(1,3)*h_temp2 + ...
        alpha(1,4)*h_temp3 + alpha(1,5)*h_temp4 + alpha(1,6)*h_temp5;        
    end
        
end 