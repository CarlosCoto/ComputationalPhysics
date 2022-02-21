
!*****************************************************************
!*		Program that simulates a quantum particle				 *
!*				on a squared potential	well	 				 *
!*****************************************************************

	program schroedinger
	
	implicit none
	double complex phi(0:1000), beta(0:1000), alpha(0:1000)
	double complex i, auxi
	real*8 V(0:1000)
	real*8 lambda, snorm, ncycles, k0, sigma, norm
	integer N, t, j
	parameter(N=1000)
	common/parameters/ ncycles,auxi

	open(111,file="alpha.dat")
	open(222,file="beta0.dat")
	open(333,file="phi0.dat")
	open(444,file="data.dat")
	open(555,file="norm.dat")

! Complex number
	i=(0.0d+0,1.0d+0)

! Input N
!	write(*,*) 'Input N'
!	read(*,*) N

! Lambda, sigma, ...
	lambda=0.001d0
	sigma=N/16.0
	ncycles=10.0
	k0=2.0*dacos(-1.0d0)*ncycles/N
	snorm=1.0/(4.0*k0**2)
	auxi=2.0*i/snorm
		
	write(444,*) lambda, sigma, ncycles, k0, snorm

! Input the potential.
	do j=0,N-1
	   V(j)=0.0d0
	   if ((j.ge.(2.0*N/5.0)).and.(j.le.(3.0*N/5.0))) then
	      V(j)=lambda
	   endif
	enddo

! Initial contitions.
	phi(0)=(0.0d+0,0.0d+0)
	phi(N)=(0.0d+0,0.0d+0)
	alpha(N-1)=(0.0d+0,0.0d+0)
	beta(N-1)=(0.0d+0,0.0d+0)

! Input initial wave function.
	do j=1,N-1,1
	   phi(j)=exp(i*k0*j)*exp(-(j-N/4.0)**2/(2.0*sigma*sigma))
	   write(333,*) j, phi(j)
	enddo

! Calculate alpha for all space that remains during sim.
	do j=N-1,1,-1
	   alpha(j-1)=-1.0/(alpha(j)-2.0+auxi-V(j))
	   write(111,*) j-1, alpha(j-1)
	enddo

! Call the method to calculate the wave function
	do t=0,10000
		call fphi(beta,alpha,phi,snorm,V,i,norm)

! Check if norm conservates.
!		write(555,*) t, norm
!		write(*,*) t, norm
	enddo
	
	close(111)
	close(222)
	close(333)
	close(444)
	close(555)

	end program

!***********************************************************************	


! Subroutine that calculates the gi's & phi's
	subroutine fphi(beta,alpha,phi,snorm,V,i,norm)	
	
	implicit none
	double complex phi(0:1000), beta(0:1000), alpha(0:1000)
	double complex i, gi, auxi
	real*8 V(0:1000)
	real*8 snorm, k0, ncycles, aux, norm, dn
	integer N, j
	parameter(N=1000)
	common/parameters/ ncycles,auxi
	
	i=(0.0d+0,1.0d+0)
	dn=1000.0
	aux=(abs(phi(0)))**2

! Re-calculate beta with the new wave function in n+ncycles
	do j=N-1,1,-1
		beta(j-1)=(-beta(j)+phi(j)*2.0*auxi)/ &
     	(alpha(j)-2.0+auxi-V(j))
		write(222,*) j-1, beta(j-1)
	enddo
	gi=0
	do j=1,N-1
	   gi=alpha(j-1)*gi+beta(j-1)
	   phi(j)=gi-phi(j)
	   aux=aux+(abs(phi(j)))**2

! Print square modulus values of wave function and V
	   write(*,*) j/dn, (abs(phi(j)))**2
	   write(*,*) j/dn, V(j)
	enddo
	aux=aux+(abs(phi(N-1)))**2
	norm=aux/dn
	write(*,*)
	
	end subroutine

!************************************************************************
