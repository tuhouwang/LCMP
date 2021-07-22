function [kr, eigvector, z] = ChebEigenValueVector(Ns, Layers,...
    interface, dep, k, rho)

    Dsave = cell(Layers, 1);
    U     = zeros(sum(Ns+1), sum(Ns+1));
    z     = cell(Layers, 1);
    
    n = 1;
    for i = 1 : Layers
        [D, x] = ChebDifferenceMatrix(Ns(i));
        z(i) = {( ( dep{i}(1) + dep{i}(end) ) / ( dep{i}(end) - dep{i}(1) ) - x )...
                                               * ( dep{i}(end) - dep{i}(1) ) / 2.0};
        A = 4.0 / ( dep{i}(end) - dep{i}(1) ) ^ 2 * diag(rho{i}) * D ...
                               * diag(1 ./ rho{i}) * D + diag(k{i}.^ 2); 
        Dsave(i) = {D};       
        
        U(n:n+Ns(i)-2,     n:n+Ns(i)-2) = A(2:Ns(i), 2:Ns(i));
        U(n:n+Ns(i)-2, sum(Ns-1)+2*i-1) = A(2:Ns(i),       1);
        U(n:n+Ns(i)-2, sum(Ns-1)+2*i  ) = A(2:Ns(i), Ns(i)+1); 
        n = n + Ns(i) - 1;
    end
    
    % boundary condition
    n = 1;
    for i = 1 : Layers - 1
        
        U(sum(Ns-1)+2*i,   sum(Ns-1)+2*i  )   =  1;
        U(sum(Ns-1)+2*i,   sum(Ns-1)+2*i+1)   = -1;
        
        if (i==1)
            left = 1 / rho{i}(end) / interface(i) .* Dsave{i}(end, :);
        else
            left = 1 / rho{i}(end) / (interface(i+1)-interface(i)) .* Dsave{i}(end, :);       
        end
        
        right   = -1 / rho{i+1}(1) / (interface(i+1)-interface(i)) .* Dsave{i+1}(1, :);   
        
        U(sum(Ns-1)+2*i+1,     n:n+Ns(i)-2)  = left(2:Ns(i));
        U(sum(Ns-1)+2*i+1, sum(Ns-1)+2*i-1)  = left(1);
        U(sum(Ns-1)+2*i+1,   sum(Ns-1)+2*i)  = left(end);

        U(sum(Ns-1)+2*i+1, n+Ns(i)-1:n+Ns(i)+Ns(i+1)-3) = right(2:Ns(i+1));
        U(sum(Ns-1)+2*i+1, sum(Ns-1)+2*i+1) = right(1);
        U(sum(Ns-1)+2*i+1, sum(Ns-1)+2*i+2) = right(end);
       
        n = n + Ns(i) - 1;
    end  
        
    U(sum(Ns-1)+1, sum(Ns-1)+1) = 1.0; %surface
    U(sum(Ns+1),     sum(Ns+1)) = 1.0; %bottom

    %blocking
    L11 = U(1          :sum(Ns-1), 1          :sum(Ns-1));   
    L12 = U(1          :sum(Ns-1), sum(Ns-1)+1:sum(Ns+1));
    L21 = U(sum(Ns-1)+1:sum(Ns+1), 1          :sum(Ns-1));
    L22 = U(sum(Ns-1)+1:sum(Ns+1), sum(Ns-1)+1:sum(Ns+1));

    L = L11 - L12 * (L22 \ L21);
    [v, k2] = eig(L);

    v2 = -(L22 \ L21) * v;
    
    eigvector = cell(Layers, 1);
    
    n = 1;
    for i = 1 : Layers
        eigvector(i) = {[v2(2*i-1, :); v(n:n+Ns(i)-2, :); v2(2*i, :)]};
        n = n + Ns(i) - 1;
    end

    k2 = sqrt(diag(k2));
    [~, ind] = sort(real(k2), 'descend');
    kr = k2(ind);
    
    for i = 1 : Layers
        eigvector(i) = {eigvector{i}(:, ind)};
    end

end
