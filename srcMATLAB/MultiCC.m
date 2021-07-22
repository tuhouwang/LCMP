% edit 'input.txt';
clear;
close all;
clc;
tic

[casename, Layers, Ns, cpmax, freq, zs, zr, rmax, dr, interface, dz, ...
      tlmin, tlmax, dep, c, rho, alpha] = ReadEnvParameter('input.txt');

[c, rho, alpha] = ChebInterpolation(dep, c, rho, alpha, Layers, Ns);

[nr, r, rhozs, k, w] = ChebInitialization(Layers, Ns, freq, ...
                  rmax, dr, zs, rho, dep, c, alpha, interface);

[kr, eigvector, z] = ChebEigenValueVector(Ns, Layers, interface, dep, k, rho);

[nmodes, kr, eigvector] = NumOfModes(Layers, w, kr, eigvector, cpmax);

[psi, psizs, z] = ChebNormalization(Layers, nmodes, eigvector, z, zs, rho);

% ShowMode(psi, z);

[tl, tl_zr, z]  = SynthesizeSoundField(r, z, dz, kr, psizs, rhozs, psi, zr);

% ShowWavenumbers(kr, casename);
% ShowTLcurve(r, zr, tl_zr);
ShowSoundField(r, z, tl, tlmin, tlmax, casename, interface);
% SaveSoundField('tl.bin', tlmin, tlmax, r, z, tl);
toc;
