function [nmodes, kr, eigvector] = NumOfModes(Layers, w, kr, eigvector, cpmax)

    cp     = w ./ real(kr);

    mode = find( cp <= cpmax );
    kr = kr(mode);
   
    for i = 1 : Layers
        eigvector(i) = {eigvector{i}(:, mode)};
    end
 
    
    mode = find( imag(kr) >= 0 & imag(kr) < real(kr) );
    kr = kr(mode);

    for i = 1 : Layers
        eigvector(i) = {eigvector{i}(:, mode)};
    end
    
    nmodes = length( kr );
    if(nmodes == 0)
        error('Incorrect maximum phase speed input!');
    end
    
end
