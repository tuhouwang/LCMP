function ShowSoundField(r, z, tl, tlmin, tlmax, casename, interface)

    disp('plot the transmission loss field!');
    r = r ./ 1000;
    figure; hold on;    
    pcolor( r , z, tl); 
    for i = 1 : length(interface) - 1
        plot([0, max(r)], [interface(i), interface(i)], 'k--', 'Linewidth', 1.5);
    end
    
    title(casename);
    caxis([tlmin tlmax]); 
    shading flat; view(0, -90);
    xlabel('Range (km)'); ylabel('Depth (m)');
    colormap(flipud(jet)); colorbar('YDir', 'Reverse');
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');

end
