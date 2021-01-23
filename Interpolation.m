function [c, rho, alpha] = Interpolation(dep, c, rho, alpha, N, s, t)

    [x, ~, ~] = LGLnodes(N);
    z = ( (t + s) / (t - s) - x ) * (t - s) / 2.0;

    c     = interp1(dep, c,     z, 'linear', 'extrap');
    rho   = interp1(dep, rho,   z, 'linear', 'extrap');
    alpha = interp1(dep, alpha, z, 'linear', 'extrap');

end
