   
function [D, x]=DifferenceMatrix(N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % DifferenceMatrix.m
    %
    % Computes the Legendre differentiation matrix with collocation at the 
    % Legendre-Gauss-Lobatto nodes.
    %
    % Reference: 
    %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
    %   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
    %
    % Written by Houwang Tu - 01/22/2021
    % Contact: tuhouwang96@163.com
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [x, ~, ~] = LGLnodes(N);
    
    % The i-th row and j-th column of L and LP stores L_j(x_i)    
    [L, ~]    = legpoly(N, x);
    D  = zeros(N + 1, N + 1);

    for k = 1 : N + 1
        for j = 1 : N + 1
            if(k ~= j)
                D(k, j) = L(k, N + 1) / L(j, N + 1) / (x(k) - x(j));
            end
        end
    end
    
    D(1, 1)         =   N * (N + 1) / 4;
    D(N + 1, N + 1) = - N * (N + 1) / 4;
    
end
