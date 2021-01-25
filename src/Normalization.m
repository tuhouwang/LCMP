function [psi, psizs, z] = Normalization(eigvectorw, eigvectorb, z1, z2, zs, rhow, rhob)

    z   = [z1(1 : length(z1) - 1); z2];
    psi = [eigvectorw(1 : size(eigvectorw, 1) - 1, :); eigvectorb];

    for k = 1 : size(eigvectorw, 2)
        f1 = eigvectorw(:, k) .^ 2;
        f2 = eigvectorb(:, k) .^ 2;

        f1 = diag(1.0 ./ rhow) * f1;
        f2 = diag(1.0 ./ rhob) * f2;

        norm1 = LGLQuadrature(f1) * max(z1) / 2;
        norm2 = LGLQuadrature(f2) * (max(z2) - max(z1)) / 2;

        psi(:, k) = psi(:, k) ./ sqrt(norm1 + norm2);
    end
    
    psizs = interp1(z, psi, zs, 'linear');

end
