function [tl, tl_zr, z] = SynthesizeSoundField(r, z, dz, kr, psizs, rhozs, psi, zr)

    bessel = besselh(0, 1, kr * r);
    p      = psi * diag( psizs ) * bessel * 1i * pi / rhozs;
    tl     = -20 * log10( abs(p) );
    tl_zr  = interp1(z, tl, zr, 'linear');
    
    % Interpolate TL to equidistant points
    zi     = 0 : dz : max(z);
    tl     = interp1(z, tl, zi, 'linear');
    z      = zi;
    
end
