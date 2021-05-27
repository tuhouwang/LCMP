! Ocean acoustic normal modes.

! Copyright (C) 2021 Houwang Tu
! -------------------------------------------------------------------------
! This program is free software: you can redistribute it and/or modify it |
! under the terms of the GNU General Public License as published by the   |
! Free Software Foundation, either version 3 of the License, or (at your  |
! option) any later version.                                              |
!                                                                         |
! This code is distributed in the hope that it will be useful, but without|
! any warranty; without even the implied warranty of merchantability or   |
! fitness for a particular purpose. See the GNU General Public License for|
! more details.                                                           |
!                                                                         |
! You should have received a copy of the GNU General Public License along |
! with this program. If not, see <http://www.gnu.org/licenses/>.          |
!                                                                         |
! Originally developed as part of the author's article (H.Tu, Y.Wang, Q.  |
! Lan et al., Applying a Legendre collocation method based on domain     |
! decomposition for calculating underwater sound propagation in a         |
! horizontally stratified environment, arXiv:2011.02850) under the        |
! supervision  of Prof. Yongxian Wang, National University of Defense     |
! Technology, China.                                                      |
!																		  |
! This Matlab/Scilab style code computes the layered and range-independent|
! modal acoustic field using the Legendre collocation spectral method    |
! based on the normal modes.                                              |
! -------------------------------------------------------------------------
!*********************************************************************

program MultiLC
	use parameters
	use function_MultiLC
	implicit none
	external zgesv
	external zgeev
	!-------------------------------------------------------------
	! Declare the variable needed later.
	character(len=MAX_FILENAME_LEN)  :: data_file = "input.txt"
	character(len=MAX_FILENAME_LEN)  :: filename  = 'tl.bin'
	character(len=MAX_FILENAME_LEN)  :: casename
	integer(rkind)                   :: nr
	integer(rkind)                   :: nmodes
	integer(rkind)                   :: Layers
	integer(rkind),allocatable       :: Ns(:)
	real(rkind)                      :: cpmax
	real(rkind)                      :: freq
	real(rkind)                      :: zs,zr,dz
	real(rkind)                      :: rmax,dr
	real(rkind)                      :: tlmin,tlmax
	real(rkind)                      :: rhozs	
	real(rkind),   allocatable,dimension(:)   :: hi,r,z
	real(rkind),   allocatable,dimension(:,:) :: tl
	complex(rkind),allocatable,dimension(:)   :: kr,psizs
	complex(rkind),allocatable,dimension(:,:) :: eigvector,psi
	real(rkind),   dimension(LL,CC)           :: dep,c,rho,alpha
	complex(rkind),dimension(LL,CC)           :: ki
	  
	call ReadEnvParameter(casename,Layers,Ns,cpmax,freq,zs,zr,rmax,dr,dz,&
							tlmin,tlmax,hi,data_file,dep,c,rho,alpha)

	call LegInitialization(freq,rmax,dr,nr,r,zs,rhozs,hi,ki,Layers,Ns,dep,c,rho,alpha)

	call LegEigenValueVector(Ns,Layers,hi,dep,rho,ki,kr,eigvector)

	call NumOfModes(Layers,freq,cpmax,kr,eigvector,nmodes)

	call LegNormalization(Layers,nmodes,eigvector,Ns,rho,dep,z,psi,zs,psizs)

	call SynthesizeSoundField(nmodes,nr,r,kr,z,psizs,rhozs,psi,tl)

	call SaveSoundField(filename,tlmin,tlmax,r,z,tl)

	deallocate(hi,r,z,tl,kr,eigvector,psi,psizs)

end program
