% edit 'input.txt';
clear;
close all;
clc;
tic
[casename, Nw, Nb, cpmax, freq, zs, zr, rmax, dr, interface, ...
    bottom, dz, Lowerboundary, tlmin, tlmax, depw, cw, rhow, ...
    alphaw, depb, cb, rhob, alphab] = ReadEnvParameter('input.txt');

[cw, rhow, alphaw] = Interpolation(depw, cw, rhow, alphaw, Nw,   0,    interface);
[cb, rhob, alphab] = Interpolation(depb, cb, rhob, alphab, Nb, interface, bottom);

[nr, r, rhozs, kw, kb, w] = Initialization(Nw, Nb, freq, rmax, ...
    dr, zs, rhow, rhob, cw, cb, alphaw, alphab, interface, bottom);

[kr, eigvectorw, eigvectorb, z1, z2] = EigenValueVector(Nw, ...
    Nb, interface, bottom, kw, kb, rhow, rhob, Lowerboundary);

[nmodes, kr, eigvectorw, eigvectorb] = NumOfModes(w, kr, ...
    eigvectorw,eigvectorb, cpmax);

[psi, psizs, z] = Normalization(eigvectorw, eigvectorb, z1, z2, zs, rhow, rhob);

% ShowMode(psi, z);

[tl, tl_zr, z]  = SynthesizeSoundField(r, z, dz, kr, psizs, rhozs, psi, zr);

% ShowWavenumbers(kr, casename);
% ShowTLcurve(r, zr, tl_zr);
ShowSoundField(r, z, tl, tlmin, tlmax, casename, interface);
% SaveSoundField('tl.bin', tlmin, tlmax, r, z, tl);
toc;
