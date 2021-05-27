function integral = LGLQuadrature(f)

    integral = 0;
    n  = length(f);
    [~, w, ~] = LGLnodes(n - 1);
    for k = 1 : n
         integral = integral + w(k) * f(k);
    end

end
