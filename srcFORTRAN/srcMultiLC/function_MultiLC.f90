module parameters

	implicit none
	integer, 		parameter :: rkind	= 8
	real(rkind),   	parameter :: pi		= 4.0d0*atan(1.0d0)
	real(rkind),   	parameter :: eps    = 2.2204e-16
	complex(rkind),	parameter :: ci		= cmplx(0.0d0,1.0d0)
	integer, 		parameter :: LL	    = 5
	integer, 		parameter :: CC	    = 2000
    integer,        parameter :: MAX_FILENAME_LEN = 200	
	
end module

module function_MultiLC
use parameters
implicit none

    interface diag
        module procedure diagreal
        module procedure diagcomplex
    end interface
	
	interface Interpolation_zs
        module procedure Interpolation_zs_R
        module procedure Interpolation_zs_C
    end interface
	
contains

subroutine ReadEnvParameter(casename,Layers,Ns,cpmax,freq,zs,zr,rmax,dr,dz,&
						tlmin,tlmax,hi,data_file,dep,c,rho,alpha)
					
	implicit none
	character(len=MAX_FILENAME_LEN), intent(out) :: casename,data_file	
	integer(rkind),                  intent(out) :: Layers
	real   (rkind),                  intent(out) :: cpmax
	real   (rkind),                  intent(out) :: freq
	real   (rkind),                  intent(out) :: zs
 	real   (rkind),                  intent(out) :: zr
	real   (rkind),                  intent(out) :: rmax
	real   (rkind),                  intent(out) :: dr
 	real   (rkind),                  intent(out) :: dz
	real   (rkind),                  intent(out) :: tlmin
	real   (rkind),                  intent(out) :: tlmax     
	integer(rkind),     allocatable, intent(out) :: Ns(:)
	real   (rkind),     allocatable, intent(out) :: hi(:)
	real   (rkind),dimension(LL,CC), intent(out) :: dep
	real   (rkind),dimension(LL,CC), intent(out) :: c
	real   (rkind),dimension(LL,CC), intent(out) :: rho
	real   (rkind),dimension(LL,CC), intent(out) :: alpha    
	integer(rkind),     allocatable	 		     :: nprofile(:)
	integer(rkind)                               :: i,j
	dep   = 0.0_rkind
	c     = 0.0_rkind
	rho   = 0.0_rkind
	alpha = 0.0_rkind
	
    open(unit=1, status='unknown', file=data_file) 
  
	read(1,*)casename
	read(1,*)Layers
	
	if(Layers > LL) then 
		stop 'Error! Layers must less than or equal to 5!'
    end if
	
	allocate(Ns(Layers),hi(Layers),nprofile(Layers))

	do i=1, Layers
        read(1,*) Ns(i)
    enddo
	read(1,*)cpmax
	read(1,*)freq
	read(1,*)zs
	read(1,*)zr	
	read(1,*)rmax
	read(1,*)dr
	if( rmax/dr - int(rmax/dr) /=0 ) then
		   stop 'Error! The input dr unsuitable !' 
	end if
	
	do i=1, Layers
        read(1,*) hi(i)
    enddo
	do i=2, Layers
        if( hi(1) /=0 .and. hi(i)<hi(i-1) ) then
		   stop 'Error! h(i) must greater than hi(i-1) !' 
		end if
    enddo	
	
	read(1,*)dz
	do i=1, Layers
        if( hi(i)/dz - int(hi(i)/dz) /=0 ) then
		   stop 'Error! The input dz unsuitable !' 
		end if
    enddo

	if( zs >= hi(Layers) .or. zs <= 0 .or. zr >= hi(Layers) .or. zr <= 0 ) then
		stop 'zs and zr must be greater than 0 and less than H !' 
	end if	
		
	read(1,*)tlmin
	read(1,*)tlmax
	if (tlmin >= tlmax) then
        stop 'tlmin must less than tlmax !'
    end if
	
	do i=1, Layers
        read(1,*) nprofile(i)
		if (nprofile(i) > CC) then
		   stop 'Error! nprofiles must less than or equal to 2500!'	
        end if		   
    enddo
	
	do i=1, Layers
		do j=1, nprofile(i)
            read(1,*) dep(i,j),	c(i,j), rho(i,j), alpha(i,j)
		enddo		
    enddo
	
	do i=2, Layers
		if (hi(i) /= dep(i,nprofile(i)) .and. hi(i-1) /= dep(i,1)) then
		   stop 'Error! input sound profile is unsuitable !'	
        end if			
    enddo

	do i=1, Layers
		call Interpolation(nprofile(i),dep(i,:),c(i,:),rho(i,:),alpha(i,:),Ns(i))	
    enddo		
		
end subroutine

subroutine Interpolation(m,dep,c,rho,alpha,N)
    implicit none
    integer(rkind), intent(in)    :: N,m
    real   (rkind), intent(inout) :: dep(N+1)
    real   (rkind), intent(inout) :: c(N+1)
    real   (rkind), intent(inout) :: rho(N+1)
    real   (rkind), intent(inout) :: alpha(N+1)
    real   (rkind)                :: x(N+1)
    real   (rkind)                :: z(N+1)
    real   (rkind)                :: b2(N+1)
    real   (rkind)                :: c2(N+1)
    real   (rkind)                :: d2(N+1)    
    integer(rkind)                :: i,j
	
	x = LGLnodes(N)
    do i=1, N+1
        z(i)=((dep(m)+dep(1))/(dep(m)-dep(1))-x(i))*(dep(m)-dep(1))/2.0_rkind
    end do
	
    do i=1, N+1	
        do j=1, m-1
            if((z(i) >= dep(j)) .and. (z(i) <= dep(j+1))) then
                b2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * c(j+1)&
                        +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * c(j)
                
                c2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * rho(j+1)&
                        +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * rho(j)
                
                d2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * alpha(j+1)&
                        +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * alpha(j)				              
            end if				
        end do 
        
		if(z(i) > dep(m)) then
			b2(i) = c(m)
			c2(i) = rho(m)
			d2(i) = alpha(m)							
		end if
			
		if(z(i) < dep(1)) then
			b2(i) = c(1)
			c2(i) = rho(1)	
			d2(i) = alpha(1)					
		endif
    end do
	
	dep   = z
	c     = b2
	rho   = c2
	alpha = d2
	       
end subroutine Interpolation

function LGLnodes(N)
    implicit none
    integer(rkind),intent(in)   :: N
    real   (rkind)              :: LGLnodes(N+1)
	real   (rkind)              :: P(N+1,N+1)
	real   (rkind)              :: xold(N+1)
	integer(rkind)              :: j,k
	
	P = 0.0_rkind
	do j=1, N + 1
        LGLnodes(j) = cos((j - 1) * pi / N)
    end do
	xold = 2.0_rkind
	
	do while( maxval(abs(LGLnodes-xold))>eps )
        xold = LGLnodes

        P(:, 1) = 1.0_rkind    
		P(:, 2) = LGLnodes

        do k=2, N
		    do j=1, N+1
				P(j,k+1)=( (2*k-1)*LGLnodes(j)*P(j,k)-(k-1)*P(j,k-1) )/k
			end do
		enddo

        do j=1, N+1
            LGLnodes(j) = xold(j)-( LGLnodes(j) * P(j,N+1)-P(j,N) )/( (N+1)*P(j,N+1) )
		end do

    end do

end function

subroutine LegInitialization(freq,rmax,dr,nr,r,zs,rhozs,hi,ki,Layers,Ns,dep,c,rho,alpha)
    implicit none
    integer(rkind), intent(in)   :: Layers
    integer(rkind), intent(in)   :: Ns(:)
    integer(rkind), intent(out)  :: nr	
    real   (rkind), intent(in)   :: dep(LL,CC)
    real   (rkind), intent(in)   :: c(LL,CC)
    real   (rkind), intent(in)   :: rho(LL,CC)
    real   (rkind), intent(in)   :: alpha(LL,CC)
    real   (rkind), intent(in)   :: freq
    real   (rkind), intent(in)   :: rmax
    real   (rkind), intent(in)   :: dr
    real   (rkind), intent(in)   :: zs
    real   (rkind), intent(in)   :: hi(:)
    real   (rkind), allocatable  :: r(:)
    real   (rkind), intent(out)  :: rhozs
    complex(rkind), intent(out)  :: ki(LL,CC)
    integer(rkind)               :: i, j
	
	nr = int(rmax / dr)
	allocate(r(nr))
	do i=1, nr
	   r(i)=i*dr
	enddo
	
	do i=1, Layers
	   r(i)=i*dr
	enddo	

    do i=1, Layers
		do j=1, Ns(i)+1
			 ki(i,j) = 2.0_rkind*pi*freq/c(i,j)*(1.0+ci*alpha(i,j)/(40.0*pi*log10(exp(1.0))))
		end do
    end	do

    if( zs <= hi(1) ) then       
        rhozs = Interpolation_zs(dep(1,1:Ns(1)+1),zs,rho(1,1:Ns(1)+1))
    end if	
	do i=2, Layers
        if( zs <= hi(i) .and. zs >= hi(i-1)) then     
            rhozs = Interpolation_zs(dep(i,1:Ns(i)+1),zs,rho(i,1:Ns(i)+1))
        end if
    end do
	
end subroutine

function Interpolation_zs_R(z,zs,v)
	implicit none
	real   (rkind), intent(in) :: z(:)
    real   (rkind), intent(in) :: zs
	real   (rkind), intent(in) :: v(:)
	real   (rkind)             :: Interpolation_zs_R
	integer(rkind)             :: j
    Interpolation_zs_R = 0.0_rkind

	if(zs>z(size(z)) .or. zs < z(1)) then
		stop 'zs is not in the range of z, interpolation failed !'
	end if
	
	if(size(z) /= size(v)) then
		stop 'The dimensions of the interpolation array do not match !'
	end if	
		
    do j=1,size(z)-1
        if (zs>=z(j).and.zs<=z(j+1)) then
            Interpolation_zs_R=(zs-z(j))/(z(j+1)-z(j))*v(j+1)+(z(j+1)-zs)/(z(j+1)-z(j))*v(j)
        endif				
    end do 

end function

function Interpolation_zs_C(z,zs,v)
	implicit none
	real   (rkind), intent(in) :: z(:)
    real   (rkind), intent(in) :: zs
	complex(rkind), intent(in) :: v(:)
	complex(rkind)             :: Interpolation_zs_C
	integer(rkind)             :: j
    Interpolation_zs_C = 0.0_rkind

	if(zs>z(size(z)) .or. zs < z(1)) then
		stop 'zs is not in the range of z, interpolation failed !'
	end if
	
	if(size(z) /= size(v)) then
		stop 'The dimensions of the interpolation array do not match !'
	end if	
	
    do j=1,size(z)-1
        if (zs>=z(j).and.zs<=z(j+1)) then
            Interpolation_zs_C=(zs-z(j))/(z(j+1)-z(j))*v(j+1)+(z(j+1)-zs)/(z(j+1)-z(j))*v(j)
        endif				
    end do 

end function

subroutine LegEigenValueVector(Ns,Layers,hi,dep,rho,ki,kr,eigvector)
    implicit none
    integer(rkind), intent(in)  :: Layers
    integer(rkind), intent(in)  :: Ns(:)
    real   (rkind), intent(in)  :: hi(:)	
    real   (rkind), intent(in)  :: dep(LL,CC)
    real   (rkind), intent(in)  :: rho(LL,CC)    
    complex(rkind), intent(out) :: ki(LL,CC)
	complex(rkind), allocatable :: kr(:)
	complex(rkind), allocatable :: A(:,:)
	complex(rkind), allocatable :: U(:,:)
	complex(rkind), allocatable :: left(:)
	complex(rkind), allocatable :: right(:)
	complex(rkind), allocatable :: eigvector(:,:)
	complex(rkind), allocatable :: L11(:,:)
	complex(rkind), allocatable :: L12(:,:)
	complex(rkind), allocatable :: L21(:,:)
	complex(rkind), allocatable :: L22(:,:)
	complex(rkind), allocatable :: L(:,:)
	complex(rkind), allocatable :: v2(:,:)    
	real   (rkind), allocatable :: D(:,:)
	real   (rkind), allocatable :: D1(:,:)    
    integer(rkind)              :: i,n
	
	!------------for zgeev and zgesv--------------------------------
	integer        :: info
    integer        :: j(1)
    integer        :: IPIV(2*size(Ns))
	complex(rkind) :: VL(sum(Ns-1))
	complex(rkind) :: VR(sum(Ns-1),sum(Ns-1)),WORK(2*sum(Ns-1))
	real(rkind)    :: RWORK(2*sum(Ns-1))
	!-----------------------------------------------------	
	
	allocate(U(sum(Ns+1),sum(Ns+1)))
	
	n = 1
	do i = 1, Layers
        allocate(D(Ns(i)+1,Ns(i)+1),A(Ns(i)+1,Ns(i)+1))
		D = LegDifferenceMatrix(Ns(i))	
		A = matmul(diag(rho(i,1:Ns(i)+1)),D)
		A = matmul(A,diag(1.0_rkind/rho(i,1:Ns(i)+1)))
		A = matmul(A,D)
        A = 4.0_rkind/(dep(i,Ns(i)+1)-dep(i,1))**2*A + diag(ki(i,1:Ns(i)+1)**2) 
        
        U(n:n+Ns(i)-2,     n:n+Ns(i)-2) = A(2:Ns(i), 2:Ns(i))
        U(n:n+Ns(i)-2, sum(Ns-1)+2*i-1) = A(2:Ns(i),       1)
        U(n:n+Ns(i)-2, sum(Ns-1)+2*i  ) = A(2:Ns(i), Ns(i)+1) 
		
        n = n + Ns(i) - 1		 
		deallocate(D,A)
	end do
	
    ! boundary condition
    n = 1
    do i = 1, Layers-1
        
        U(sum(Ns-1)+2*i,   sum(Ns-1)+2*i  )   =  1.0_rkind
        U(sum(Ns-1)+2*i,   sum(Ns-1)+2*i+1)   = -1.0_rkind
		
		allocate(D(Ns(i)+1,Ns(i)+1),D1(Ns(i+1)+1,Ns(i+1)+1),left(Ns(i)+1),right(Ns(i+1)+1))
		D  = LegDifferenceMatrix(Ns(i))	
		D1 = LegDifferenceMatrix(Ns(i+1))
        
        if (i==1) then
            left = 1.0_rkind/rho(i,Ns(i)+1)/hi(i)*D(Ns(i)+1, :)
        else
            left = 1.0_rkind/rho(i,Ns(i)+1)/(hi(i+1)-hi(i))*D(Ns(i)+1, :)      
        end if
        
        right = -1.0_rkind/rho(i+1,1)/(hi(i+1)-hi(i))*D1(1, :) 
        
        U(sum(Ns-1)+2*i+1,     n:n+Ns(i)-2)  = left(2:Ns(i))
        U(sum(Ns-1)+2*i+1, sum(Ns-1)+2*i-1)  = left(1)
        U(sum(Ns-1)+2*i+1,   sum(Ns-1)+2*i)  = left(Ns(i)+1)

        U(sum(Ns-1)+2*i+1, n+Ns(i)-1:n+Ns(i)+Ns(i+1)-3) = right(2:Ns(i+1))
        U(sum(Ns-1)+2*i+1, sum(Ns-1)+2*i+1) = right(1)
        U(sum(Ns-1)+2*i+1, sum(Ns-1)+2*i+2) = right(Ns(i+1)+1)
       
        n = n + Ns(i) - 1
		deallocate(D,D1,left,right)
    end do	
	
	U(sum(Ns-1)+1, sum(Ns-1)+1) = 1.0_rkind !surface
    U(sum(Ns+1),     sum(Ns+1)) = 1.0_rkind !bottom
	
	allocate(L11(sum(Ns-1),sum(Ns-1)),L12(sum(Ns-1),2*size(Ns)),L21(2*size(Ns),sum(Ns-1)),&
	                      L22(2*size(Ns),2*size(Ns)),L(sum(Ns-1),sum(Ns-1)),kr(sum(Ns-1)),&
	        eigvector(sum(Ns+1),sum(Ns-1)),v2(2*size(Ns),sum(Ns-1)))
    !blocking
    L11 = U(1          :sum(Ns-1), 1          :sum(Ns-1))   
    L12 = U(1          :sum(Ns-1), sum(Ns-1)+1:sum(Ns+1))
    L21 = U(sum(Ns-1)+1:sum(Ns+1), 1          :sum(Ns-1))
    L22 = U(sum(Ns-1)+1:sum(Ns+1), sum(Ns-1)+1:sum(Ns+1))	
	
	call zgesv(2*size(Ns),sum(Ns-1),L22,2*size(Ns),IPIV,L21,2*size(Ns),info)	
	L = L11-matmul(L12,L21)
	call zgeev('N','V',sum(Ns-1),L,sum(Ns-1),kr,VL,1,VR,sum(Ns-1),WORK,2*sum(Ns-1),RWORK,info)  	
	v2= -matmul(L21,VR)
	
	deallocate(U,L11,L12,L21,L22,L)
	
	n = 0
    do i = 1, Layers
		eigvector(n+i              ,:) = v2(2*i-1, :)
		eigvector(n+i+1:n+Ns(i)+i-1,:) = VR(n-i+2:n+Ns(i)-i, :)
		eigvector(n+Ns(i)+i        ,:) = v2(2*i, :)
		n = n + Ns(i)			
    end do	
	
	kr = sqrt(kr)
	
	!L store the sorted eigenvectors respectively. VL stores the sorted eigenvalus
	allocate(L(sum(Ns+1),sum(Ns-1)))
	do i=1, sum(Ns-1)
	    j       =maxloc(real(kr))
		VL(i)   =kr(j(1))
		L(:,i)  =eigvector(:,j(1))
		kr(j(1))=-1.0_rkind
    enddo
	
	kr        = VL
	eigvector = L		
	deallocate(L,v2)	

end subroutine

function diagreal(rho)
    implicit none
    real   (rkind), intent(in) :: rho(:)
    real   (rkind)             :: diagreal(size(rho),size(rho))
	integer(rkind)             :: i
	diagreal = 0.0_rkind
	forall(i=1:size(rho)) diagreal(i,i)=rho(i)
	
end function

function diagcomplex(rho)
    implicit none
    complex(rkind), intent(in) :: rho(:)
    complex(rkind)             :: diagcomplex(size(rho),size(rho))
	integer(rkind)             :: i
	diagcomplex = 0.0_rkind
	forall(i=1:size(rho)) diagcomplex(i,i)=rho(i)
	
end function

function LegDifferenceMatrix(N)
    implicit none
    integer(rkind), intent(in) ::N
    real   (rkind)			   ::LegDifferenceMatrix(N+1,N+1)
    real   (rkind)             ::x(N+1),L(N+1,N+1)	
    integer(rkind)             ::j,k
	
	x = LGLnodes(N)
    L = legpoly(N, x)

	LegDifferenceMatrix = 0.0_rkind
	
    do k = 1 , N + 1
        do j = 1 , N + 1
            if(k .ne. j) then
                LegDifferenceMatrix(k, j) = L(k, N + 1) / L(j, N + 1) / (x(k) - x(j))
            end if
        end do
    end	do
	
    LegDifferenceMatrix(1, 1)         =   N * (N + 1) / 4
    LegDifferenceMatrix(N + 1, N + 1) = - N * (N + 1) / 4	
	
end function

function legpoly(N, x)
    implicit none
    integer(rkind), intent(in) :: N
    real   (rkind), intent(in) :: x(:)
	real   (rkind)             :: legpoly(size(x),N+1)	
    integer(rkind)             :: j,k
	
	legpoly = 0.0_rkind
	
	legpoly(:,1) = 1.0_rkind
    legpoly(:,2) = x
	
    do k = 2, N
		do j = 1,size(x)
			legpoly(j,k+1) = (2.0_rkind*k-1)/k*x(j)*legpoly(j,k) &
			                   - (k-1.0_rkind)/k*legpoly(j,k-1)
		end do
    end	do
	
end function

subroutine NumOfModes(Layers,freq,cpmax,kr,eigvector,nmodes)
    implicit none
    integer(rkind), intent(in)    :: Layers
    integer(rkind), intent(inout) :: nmodes	
    real   (rkind), intent(in)    :: freq
    real   (rkind), intent(in)    :: cpmax
	real   (rkind), allocatable   :: cp(:)
    complex(rkind), intent(inout) :: kr(:)
    complex(rkind), intent(inout) :: eigvector(:,:)
    integer(rkind)                :: i
	
	allocate(cp(size(kr)))
	cp = 2 * pi * freq / real(kr)
	
	where(cp>cpmax) kr = -1.0_rkind 
	where( imag(kr) < 0 .and. imag(kr) >= real(kr) ) kr = -1.0_rkind 
	
	nmodes = 0
	do i=1,size(kr)
	    if (kr(i) /= -1.0_rkind) then
		   nmodes              = nmodes + 1
		   kr       (nmodes)   =kr(i)
		   eigvector(:,nmodes) = eigvector(:,i)		   
		end if
    end do

end subroutine

subroutine LegNormalization(Layers,nmodes,eigvector,Ns,rho,dep,z,psi,zs,psizs)
    implicit none
    integer(rkind), intent(in)    :: Layers
    integer(rkind), intent(in)    :: nmodes
    integer(rkind), intent(in)    :: Ns(:)
	real   (rkind), intent(in)    :: rho(LL,CC)
	real   (rkind), intent(in)    :: dep(LL,CC)
 	real   (rkind), intent(in)    :: zs   
	real   (rkind), allocatable   :: z(:)
    complex(rkind), intent(inout) :: eigvector(:,:)
	complex(rkind), allocatable   :: psi(:,:)
	complex(rkind), allocatable   :: psizs(:)
	complex(rkind), allocatable   :: f(:)    
	real   (rkind)                :: norm
    integer(rkind)                :: i,j,k
	
	allocate(z(sum(Ns)+1),psi(sum(Ns)+1,nmodes),psizs(nmodes))
	
	j = 0
	do i=1, Layers            
        z  (j+1:j+Ns(i)  ) = dep(i,1:Ns(i))
        psi(j+1:j+Ns(i),:) = eigvector(j+1:j+Ns(i),1:nmodes)
		j = j + Ns(i)
    end do
    z  (sum(Ns)+1)   = dep(Layers,Ns(Layers)+1)
    psi(sum(Ns)+1,:) = eigvector(sum(Ns+1),1:nmodes)     
	
    do k = 1 , nmodes
        norm = 0.0_rkind
		j    = 0
        do i = 1 , Layers
           allocate(f(Ns(i)+1))		
           f = eigvector(j+1:j+Ns(i)+1, k) ** 2
           f = matmul(diag(1.0_rkind / rho(i,1:Ns(i)+1)), f)
           norm = norm + LGLQuadrature(f) * (maxval(dep(i,1:Ns(i)+1))-minval(dep(i,1:Ns(i)+1))) / 2			
		   j = j + Ns(i) + 1
		   deallocate(f)		   
        end do
        psi(:, k) = psi(:, k) / sqrt(norm)			
		psizs( k) = Interpolation_zs(z,zs,psi(:,k))
    end do

end subroutine

function LGLQuadrature(f)
    implicit none
    complex(rkind), intent(in)   :: f(:)
	complex(rkind)               :: LGLQuadrature
	real   (rkind), allocatable  :: w(:)
	integer(rkind)               :: j,k
	k = size(f)
	allocate(w(k))
	w = LGLweights(k-1)
	LGLQuadrature = 0.0_rkind
    do j = 1,k
         LGLQuadrature = LGLQuadrature + w(j) * f(j)
    end do
    
    deallocate(w)	
end function

function LGLweights(N)
    implicit none
    integer(rkind), intent(in) :: N
    real   (rkind)             :: x(N+1)
    real   (rkind)             :: LGLweights(N+1)
    real   (rkind)             :: P(N+1,N+1)
    real   (rkind)             :: xold(N+1)    
	integer                    :: j,k
	
	P = 0.0_rkind
	do j = 1, N + 1
        x(j) = cos((j - 1) * pi / N)
    end do
	xold = 2.0_rkind
	
	do while( maxval(abs(x-xold))>eps )
        xold = x

        P(:, 1) = 1.0_rkind    
		P(:, 2) = x

        do k=2, N
		    do j = 1, N + 1
				P(j,k+1)=( (2*k-1)*x(j)*P(j,k)-(k-1)*P(j,k-1) )/k
			end do
		enddo

        do j=1, N + 1
            x(j) = xold(j)-( x(j) * P(j,N+1)-P(j,N) )/( (N+1)*P(j,N+1) )
		end do

    end do
	
	LGLweights = 2.0_rkind/(N*(N+1)*P(:,N+1)**2)

end function

subroutine SynthesizeSoundField(nmodes,nr,r,kr,z,psizs,rhozs,psi,tl)
	implicit none
	integer(rkind), intent(in) :: nmodes
	integer(rkind), intent(in) :: nr    
	real   (rkind), intent(in) :: r(nr)
	real   (rkind), intent(in) :: z(:)
	real   (rkind), intent(in) :: rhozs    
	complex(rkind), intent(in) :: kr(:)
	complex(rkind), intent(in) :: psizs(nmodes)
	complex(rkind), intent(in) :: psi(:,:)    
	real   (rkind), intent(out),allocatable ::tl(:,:)
	complex(rkind)             :: p(size(z),nr)
	complex(rkind)             :: bessel(nmodes,nr)
	integer                    :: IERR1
	integer                    :: IERR2     
	real   (rkind)             :: CYR
	real   (rkind)             :: CYI    
	integer(rkind)             :: i,k
	
	
    allocate(tl(size(z),nr))
    do k=1, nmodes
	   do i=1, nr
	       bessel(k,i)=r(i)*kr(k)
		   call ZBESH(real(bessel(k,i)),aimag(bessel(k,i)),0.0d0,1,1,1,CYR,CYI,IERR1,IERR2)
		   bessel(k,i)=cmplx(CYR,CYI)
	    enddo
	enddo
	
	p  = matmul(matmul(psi,diag(psizs)),bessel)*ci*pi/rhozs
	tl = -20.0_rkind*log10(abs(p))

end subroutine

subroutine SaveSoundField(filename,tlmin,tlmax,r,z,tl)
	implicit none
	character(len=MAX_FILENAME_LEN), intent(in) :: filename
	real(rkind),                     intent(in) :: tlmin
	real(rkind),                     intent(in) :: tlmax
	real(rkind),                     intent(in) :: r(:)
	real(rkind),                     intent(in) :: z(:)
	real(rkind),                     intent(in) :: tl(:,:)    

    open(unit=2,status='unknown',file=filename,access='stream',form='unformatted')
       
	   write(2)  size(z),size(r),tlmin,tlmax,z,r,tl
     
	close(2)
	  
end subroutine

end module
