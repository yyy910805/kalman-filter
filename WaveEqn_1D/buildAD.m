function DD = buildAD(m,order)
% Constructs artificial dissipation based on undivided derivatives
% of order 'order'.
    
switch order
        
    case 1

        DD = zeros(m,m);
        d1 = [-1 1];
        DD(1,1:2) = d1;
        DD(2,1:2) = d1;
        for i=3:m-1
          DD(i,(i-1):i) = d1;
        end
        DD(m,m-1:m) = d1;
        
       
    case 2

        dd2=[1 -2 1];
        DD=(diag(ones(m-1,1),-1)-2*diag(ones(m,1),0)+ ...
              diag(ones(m-1,1),1));
        DD(1,1:3)=dd2; DD(m,m-2:m)=dd2;
        
    case 3

        d3=[-1 3 -3 1];
        DD=(-diag(ones(m-2,1),-2)+3*diag(ones(m-1,1),-1)-3*diag(ones(m,1),0)+ ...
              diag(ones(m-1,1),1));
        DD(1:2,1:4)=[d3;d3];
        DD(m,m-3:m)=d3;
        
    case 4

        d4=[1 -4 6 -4 1];
        DD=(diag(ones(m-2,1),2)-4*diag(ones(m-1,1),1)+6*diag(ones(m,1),0)-4*diag(ones(m-1,1),-1)+diag(ones(m-2,1),-2));
        DD(1:2,1:5)=[d4;d4]; DD(m-1:m,m-4:m)=[d4;d4];
        
    case 5

        d5=[-1 5 -10 10 -5 1];
        DD=(-diag(ones(m-3,1),-3)+5*diag(ones(m-2,1),-2)-10*diag(ones(m-1,1),-1)+10*diag(ones(m,1),0)-5*diag(ones(m-1,1),1)+diag(ones(m-2,1),2));
        DD(1:3,1:6)=[d5;d5;d5];
        DD(m-1:m,m-5:m)=[d5;d5];
        
    case 6

        d6=[1 -6 15 -20 15 -6 1];
        DD=(diag(ones(m-3,1),3)-6*diag(ones(m-2,1),2)+15*diag(ones(m-1,1),1)-20*diag(ones(m,1),0)+15*diag(ones(m-1,1),-1)-6*diag(ones(m-2,1),-2)+diag(ones(m-3,1),-3));
        DD(1:3,1:7)=[d6;d6;d6];DD(m-2:m,m-6:m)=[d6;d6;d6];
        
    case 7

        d7=[-1 7 -21 35 -35 21 -7 1]; 
        DD=(-diag(ones(m-4,1),-4)+7*diag(ones(m-3,1),-3)-21*diag(ones(m-2,1),-2)+35*diag(ones(m-1,1),-1)-35*diag(ones(m,1),0)+21*diag(ones(m-1,1),1)-7*diag(ones(m-2,1),2)+diag(ones(m-3,1),3));
        DD(1:4,1:8)=[d7;d7;d7;d7];
        DD(m-2:m,m-7:m)=[d7;d7;d7];
        
    case 9

        d9=[-1 9 -36 84 -126 126 -84 36 -9 1]; 
        DD=(-diag(ones(m-5,1),-5)+9*diag(ones(m-4,1),-4)-36*diag(ones(m-3,1),-3)+84*diag(ones(m-2,1),-2)-126*diag(ones(m-1,1),-1)+126*diag(ones(m,1),0)-84*diag(ones(m-1,1),1)+36*diag(ones(m-2,1),2)-9*diag(ones(m-3,1),3)+diag(ones(m-4,1),4));
        DD(1:5,1:10)=[d9;d9;d9;d9;d9];
        DD(m-3:m,m-9:m)=[d9;d9;d9;d9];
        
    case 40
    % Fourth derivative for the optimal 8th order operator
    % Nonequidistant grid.
        
        mmm=8;
        DD=(diag(ones(m-2,1),2)-4*diag(ones(m-1,1),1)+6*diag(ones(m,1),0)-4*diag(ones(m-1,1),-1)+diag(ones(m-2,1),-2));
        DD(1:6,1:8)=[0.70921010190504348684e1 -0.14196080536361841322e2 0.11072881931325435634e2 -0.50473576941871051066e1 0.10784552801730759259e1 0 0 0; 0.70921010190504348684e1 -0.14196080536361841322e2 0.11072881931325435634e2 -0.50473576941871051066e1 0.10784552801730759259e1 0 0 0; 0.70921010190504348684e1 -0.14196080536361841322e2 0.11072881931325435634e2 -0.50473576941871051066e1 0.10784552801730759259e1 0 0 0; 0 0.13740993382151221352e1 -0.42105600869792757010e1 0.54761010136211975317e1 -0.35797005751940657417e1 0.94006031033702177578e0 0 0; 0 0 0.82467928104463767301e0 -0.33274694995849432461e1 0.52587584638857303123e1 -0.37020511582893568152e1 0.94608291294393207601e0 0; 0 0 0 0.86436129166612654748e0 -0.37325441295306179390e1 0.57924699560798105338e1 -0.39066885960487908497e1 0.98240147783347170744e0;];
        DD(m-mmm+1:m,m-mmm+1:m)=rot90( DD(1:mmm,1:mmm) ,2 );
end
    
DD = sparse(DD);
        
end