% Ocean acoustic normal modes.

% Copyright (C) 2021 Houwang Tu
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
% Originally developed as the author's article (H.Tu, Y. Wang, Q. Lan et  |
% al., Applying a Legendre collocation method based on domain             |
% decomposition to calculate underwater sound propagation in a            |
% horizontally stratified environment, Journal of Sound and Vibration,    |
% https://doi.org/10.1016/j.jsv.2021.116364) under the supervision of     |
% Prof. Yongxian Wang, National University of Defense Technology, China.  |
%																		  |
% This Matlab/Scilab style code computes the horizontally stratified and  |
% range-independent modal acoustic field using the Legendre collocation   | 
% method based on the normal modes.                                                       |
% -------------------------------------------------------------------------
% edit 'input.txt';
clear;
close all;
clc;
tic

[casename, Layers, Ns, cpmax, freq, zs, zr, rmax, dr, interface, dz, ...
      tlmin, tlmax, dep, c, rho, alpha] = ReadEnvParameter('input.txt');

[c, rho, alpha] = LegInterpolation(dep, c, rho, alpha, Layers, Ns);

[nr, r, rhozs, k, w] = LegInitialization(Layers, Ns, freq, rmax, ...
                             dr, zs, rho, dep, c, alpha, interface);

[kr, eigvector, z] = LegEigenValueVector(Ns, Layers, interface, dep, k, rho);

[nmodes, kr, eigvector] = NumOfModes(Layers, w, kr, eigvector, cpmax);

[psi, psizs, z] = LegNormalization(Layers, nmodes, eigvector, z, zs, rho);

% ShowMode(psi, z);

[tl, tl_zr, z]  = SynthesizeSoundField(r, z, dz, kr, psizs, rhozs, psi, zr);

% ShowWavenumbers(kr, casename);
% ShowTLcurve(r, zr, tl_zr);
ShowSoundField(r, z, tl, tlmin, tlmax, casename, interface);
% SaveSoundField('tl.bin', tlmin, tlmax, r, z, tl);
toc;
