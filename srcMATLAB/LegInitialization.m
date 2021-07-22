function [nr, r, rhozs, k, w] = LegInitialization(Layers, Ns, freq, rmax, ...
                                       dr, zs, rho, dep, c, alpha, interface)

    w  = 2 * pi * freq;
    r  = dr : dr : rmax;
    nr = length(r);

    for i = 1 : Layers
        if( zs <= interface(i) )
            [x, ~, ~] = LGLnodes( Ns(i) );
            z = ( (dep{i}(1) + dep{i}(end) ) / ( dep{i}(end) - dep{i}(1)) ...
                                     - x ) * (dep{i}(end) - dep{i}(1)) / 2.0;
        
            rhozs = interp1(z, rho{i}, zs, 'linear');
            break
        end
    end
    
    k   = cell(Layers, 1);  
    for i = 1 : Layers
        k(i) = {w ./ c{i} .* (1.0 + 1i * alpha{i} / (40.0 * pi * log10(exp(1.0))))};
    end
end
