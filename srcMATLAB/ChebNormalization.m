function [psi, psizs, zl] = ChebNormalization(Layers, nmodes, eigvector, z, zs, rho)

   for i = 1 : Layers            
       if (i == 1)
           zl  = z{i};
           psi = eigvector{i};
       else
           zl  = [zl; z{i}(2 : end)];
           psi = [psi; eigvector{i}(2 : end, :)];
       end
   end

    for k = 1 : nmodes
        norm = 0.0;
        for i = 1 : Layers 
            f    = eigvector{i}(:, k) .^ 2;
            f    = diag(1.0 ./ rho{i}) * f;
            norm = norm + ChebGaussQuadrature(f) * (max(z{i}) - min(z{i})) / 2;
        end
        psi(:, k) = psi(:, k) ./ sqrt(norm);
    end
    
    psizs = interp1(zl, psi, zs, 'linear');

end
