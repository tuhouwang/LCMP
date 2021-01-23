% Ocean acoustic normal modes.

% Copyright (C) 2020 Houwang Tu
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify it |
% under the terms of the GNU General Public License as published by the   |
% Free Software Foundation, either version 3 of the License, or (at your  |
% option) any later version.                                              |
%                                                                         |
% This code is distributed in the hope that it will be useful, but without|
% any warranty; without even the implied warranty of merchantability or   |
% fitness for a particular purpose. See the GNU General Public License for|
% more details.                                                           |
%                                                                         |
% You should have received a copy of the GNU General Public License along |
% with this program. If not, see <http://www.gnu.org/licenses/>.          |
%                                                                         |
% Originally developed as part of the author's article (H.Tu, Y.Wang, Q.  |
% Lan et al., Applying a Legendre collocation method based on domain      |
% decomposition for calculating underwater sound propagation in a         |
% horizontally stratified environment, arXiv:2011.02850) under the        |
% supervision of Prof. Yongxian Wang, National University of Defense      |
% Technology, China.                                                      |
%									                                      |
% This Matlab/Scilab style code computes the layered and range-independent|
% modal acoustic field using the Legendre collocation spectral method     |
% based on the normal modes.                                              |
% -------------------------------------------------------------------------
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
