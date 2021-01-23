function [nmodes, kr, v1, v2] = NumOfModes(w, kr, v1, v2, cpmax)

    cp     = w ./ real(kr);

    nmodes = length( find( cp <= cpmax ) );

    if(nmodes == 0)
        error('Incorrect maximum phase speed input!');
    end

    kr = kr(   1 : nmodes);
    v1 = v1(:, 1 : nmodes);
    v2 = v2(:, 1 : nmodes);
    
%     mode = find( imag(kr) >= 0 ) ;
%     nmodes = length( find( imag(kr) >= 0 ) );
%     kr = kr(   mode);
%     v1 = v1(:, mode);
%     v2 = v2(:, mode);
    
end
