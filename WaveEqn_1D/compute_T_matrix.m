function T_matrix = compute_T_matrix(N,alpha,M,D,SAT,A,cdiss)
    
    % At a given time, 
    % model vector is defined as:
    % [flux1,flux2,...,fluxN+1,height1,...,heightN+1]

    q_temp = zeros(N+1,1); % temporary q
    h_temp = zeros(N+1,1); % temporary h
    
    h = zeros(N+1,1); % set height to 0 - find the first N+1 columns of T

    for k=1:N+1 
        
       q(:,1) = zeros(N+1,1); % set flux to 0
       q(k,1) = 1; % set the k-th component to 1
       [q_temp,h_temp] = T_unstaggered_op(0,0,q(:,1),h(:,1),alpha,M,D,SAT,A,cdiss);
       T_matrix(1:N+1,k) = q_temp;
       T_matrix(N+2:2*N+2,k) = h_temp;       
 
    end

    q = zeros(N+1,1); % set flux to 0 - find the last N+1 columns of T   

    for k=1:N+1  
        
       h(:,1) = zeros(N+1,1); % set flux to 0
       h(k,1) = 1; % set the k-th component to 1  
       [q_temp,h_temp] = T_unstaggered_op(0,0,q(:,1),h(:,1),alpha,M,D,SAT,A,cdiss);        
       T_matrix(1:N+1,k+N+1) = q_temp;
       T_matrix(N+2:2*N+2,k+N+1) = h_temp; 
       
    end
    
end 