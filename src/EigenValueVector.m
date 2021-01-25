function [kr, eigvectorw, eigvectorb, z1, z2] = EigenValueVector(Nw, Nb, ...
    interface, bottom, kw, kb, rhow, rhob, Lowerboundary)

    [D1, x1] = DifferenceMatrix(Nw);
    [D2, x2] = DifferenceMatrix(Nb);

    z1 = (1 - x1) * interface   / 2.0;
    z2 = (1 - x2) * (bottom - interface) / 2.0 + interface;   

    A  = 4.0 / interface ^ 2 * diag(rhow) * D1 * diag(1 ./ rhow) * D1 + diag(kw .^ 2);
    B  = 4.0 / (bottom - interface) ^2 * diag(rhob) * D2 * diag(1 ./ rhob) * D2 + diag(kb .^ 2);

    U = zeros(Nw+Nb+2, Nw+Nb+2);
    U(1:Nw-1,  1:Nw-1) = A(2:Nw, 2:Nw);
    U(1:Nw-1, Nw+Nb-1) = A(2:Nw,    1);
    U(1:Nw-1,   Nw+Nb) = A(2:Nw, Nw+1);

    U(Nw:Nw+Nb-2, Nw:Nw+Nb-2) = B(2:Nb, 2:Nb);
    U(Nw:Nw+Nb-2,    Nw+Nb+1) = B(2:Nb,    1);
    U(Nw:Nw+Nb-2,    Nw+Nb+2) = B(2:Nb, Nb+1);

    U(Nw+Nb-1, Nw+Nb-1) =  1.0;
    U(Nw+Nb,     Nw+Nb) =  1.0;
    U(Nw+Nb,   Nw+Nb+1) = -1.0;

    left  =  1 / rhow(Nw + 1) / interface .* D1(Nw + 1, :);
    right = -1 / rhob(1) / (bottom - interface) .* D2(1, :);

    U(Nw+Nb+1,  1:Nw-1)  = left(2:Nw);
    U(Nw+Nb+1, Nw+Nb-1)  = left(1);
    U(Nw+Nb+1,   Nw+Nb)  = left(Nw+1);

    U(Nw+Nb+1, Nw:Nw+Nb-2) = right(2:Nb);
    U(Nw+Nb+1,    Nw+Nb+1) = right(1);
    U(Nw+Nb+1,    Nw+Nb+2) = right(Nb+1);

    if(Lowerboundary==0)    %Dirichlet condition
        U(Nw+Nb+2, Nw+Nb+2) = 1.0;
    else                    %Neumann condition
        U(Nw+Nb+2, Nw:Nw+Nb-2) = D2(Nb+1, 2:Nb);
        U(Nw+Nb+2,    Nw+Nb+1) = D2(Nb+1,    1);
        U(Nw+Nb+2,    Nw+Nb+2) = D2(Nb+1, Nb+1);
    end

    %blocking
    L11 = U(1      :Nw+Nb-2, 1      :Nw+Nb-2);
    L12 = U(1      :Nw+Nb-2, Nw+Nb-1:Nw+Nb+2);
    L21 = U(Nw+Nb-1:Nw+Nb+2, 1      :Nw+Nb-2);
    L22 = U(Nw+Nb-1:Nw+Nb+2, Nw+Nb-1:Nw+Nb+2);

    L = L11 - L12 * (L22 \ L21);
    [v, k2] = eig(L);

    v2 = -(L22 \ L21) * v;

    eigvectorw = [v2(1, :); v(1:Nw-1, :)    ; v2(2, :)];
    eigvectorb = [v2(3, :); v(Nw:Nw+Nb-2, :); v2(4, :)];

    k2 = sqrt(diag(k2));
    [~, ind] = sort(real(k2), 'descend');
    kr = k2(ind);

    eigvectorw = eigvectorw(:, ind);
    eigvectorb = eigvectorb(:, ind);
end
