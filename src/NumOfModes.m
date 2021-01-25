function [nmodes, kr, v1, v2] = NumOfModes(w, kr, v1, v2, cpmax)

    cp     = w ./ real(kr);

    mode = find( cp <= cpmax );
    kr = kr(   mode);
    v1 = v1(:, mode);
    v2 = v2(:, mode);
    
    mode = find( imag(kr) >= 0 & imag(kr) < real(kr) );
    kr = kr(   mode);
    v1 = v1(:, mode);
    v2 = v2(:, mode);
    
    nmodes = length( kr );

    if(nmodes == 0)
        error('Incorrect maximum phase speed input!');
    end
    
end
