

!*****************************************************************
!*     Program that simulates the Solar System with the         *
!*     motion equations of the planets and the Sun,               *
!*     using the Verlet algorithm for differential equations     *
!*****************************************************************

      PROGRAM SOLARSYSTEM

! Variables declaration
	IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 L, SUNMASS
      REAL*8 MASS(10), RADIUS(10), DIST0(10), SPEED0(10)
      REAL*8 DIST(10,2), SPEED(10,2), W(10,2), F(10,2)
      INTEGER STEP
	PARAMETER (N=10)
	COMMON MASS
	
! Constant values (C is the astronomical unit).
      SUNMASS=1.99D+33
      C=149.6D+0
      G=6.6674D-14
      GM=132712.0D+6

	OPEN(111,FILE='positions.dat')
	OPEN(222,FILE='momentum.dat')
	OPEN(666,FILE='w.dat')
	OPEN(555,FILE='MASS.dat')
	OPEN(444,FILE='force.dat')
	OPEN(333,FILE='SPEED.dat')

! Initialize variables.
	DO I=1,N
		MASS(I)=0.0D0
		RADIUS(I)=0.0D0
		DIST0(I)=0.0D0
		SPEED0(I)=0.0D0
		DO J=1,2
			DIST(I,J)=0.0D0
			SPEED(I,J)=0.0D0
			W(I,J)=0.0D0
			F(I,J)=0.0D0
		ENDDO
	ENDDO

! Input planets and Sun masses
! .Sun, 2.Mercury, 3.Venus,
! 4.Earth, 5.Mars, 6.Jupiter, 7.Saturn, 8.Uranu, 9.Neptunus,
! 10.Pluto.

      MASS(1) = SUNMASS
      MASS(2) = 0.3302D+24
      MASS(3) = 4.8685D+24
      MASS(4) = 5.9736D+24
      MASS(5) = 0.64185D+24
      MASS(6) = 1898.6D+24
      MASS(7) = 568.46D+24
      MASS(8) = 86.832D+24
      MASS(9) = 102.43D+24
      MASS(10) = 0.0125D+24      

! Let's input the initial speeds (in 10⁶ kms unit)
! and the initial speed.
      DIST0(1) = 0.00D0
      DIST0(2) = 69.82D0
      DIST0(3) = 108.94D0
      DIST0(4) = 152.10D0
      DIST0(5) = 249.23D0
      DIST0(6) = 816.62D0
      DIST0(7) = 1514.5D0
      DIST0(8) = 3003.62D0
      DIST0(9) = 4545.67D0
      DIST0(10) = 7375.93D0

      SPEED0(1) = 0.00D0
      SPEED0(2) = 38.86D0
      SPEED0(3) = 34.79D0
      SPEED0(4) = 29.29D0
      SPEED0(5) = 21.97D0
      SPEED0(6) = 12.44D0
      SPEED0(7) = 9.09D0
      SPEED0(8) = 6.49D0
      SPEED0(9) = 5.37D0
      SPEED0(10) = 3.71D0

! RADIUS of the planets ( 10⁶kms.)
	RADIUS(1) = 1392000D-6
	RADIUS(2) = 4879D-6
	RADIUS(3) = 12104D-6
	RADIUS(4) = 12756D-6
	RADIUS(5) = 6794D-6
	RADIUS(6) = 142984D-6
	RADIUS(7) = 120536D-6
	RADIUS(8) = 51118D-6
	RADIUS(9) = 49528D-6
	RADIUS(10) = 2390D-6

! Rescale time.
	C3=(C*1.0D+6)**3
	TIME=DSQRT(GM/C3)
!	WRITE(*,*) 'Scaled time= ', TIME
!	WRITE(*,*)
	
! Sending to file the initial positions, distance and speed
! 
	DO I=1,N
!		WRITE(*,*) DIST0(I), 0.0D0, I
!		WRITE(111,*) DIST0(I)/C, 0.00D0, I
		WRITE(*,*) DIST0(I)/C, 0.0D0, I
		DIST(I,1)=DIST0(I)/C
		SPEED(I,2)=SPEED0(I)/(C*1D+6*TIME)
		
		RADIUS(I)=RADIUS(I)/C
		MASS(I)=MASS(I)/SUNMASS
!		WRITE(555,*) I, MASS(I)
	ENDDO
	WRITE(111,*)
	WRITE(*,*)

!**********     Verlet algorithm     **********

! Ask for temporal step
!	WRITE(*,*) 'Introduzca el STEP temporal en días terrestres:'
!	READ(*,*) STEP
!	TH=STEP*24*3600*TIME
	TH=24*3600*TIME
	
!	WRITE(*,*) 'Please enter number of years:'
!	READ(*,*) STEP
!	TFIN=STEP*365*24*3600*TIME
	TFIN=50*365*24*3600*TIME

!	WRITE(*,*)
!	WRITE(*,*) 'TH=', TH, 'TFIN=', TFIN
!	WRITE(*,*) 'Percentage=', (100.0-(TFIN-TH)*100.0D0/TFIN),'%'
!	WRITE(*,*)

! Calling function force'.
	CALL FORCE(F,DIST)

! Calculating W
	DO I=1,N
		DO J=1,2
			W(I,J)=SPEED(I,J)+(TH/2)*F(I,J)
		ENDDO
		WRITE(666,*) W(I,1), W(I,2), I
	ENDDO
!	GOTO 100

! Here the big loop begins.
	T=0.0D0
	DO WHILE (T.LE.TFIN)

! Postion after STEP, r(t+h).
	DO I=1,N 
		DO J=1,2
			DIST(I,J)=DIST(I,J)+TH*W(I,J)
		ENDDO
		WRITE(*,*) DIST(I,1), DIST(I,2), I
!		WRITE(111,*) DIST(I,1), DIST(I,2), I
	ENDDO
!	WRITE(111,*)
	WRITE(*,*)
	
! Stop the program if planets crash!!!
	DO I=1,N
		DO J=1,N
			IF (I.EQ.J) GOTO 300
			AUXX=0D0
			AUXY=0D0
			AUXRADIUS=0.0D0
			AUXX=DIST(I,1)-DIST(J,1)
			AUXY=DIST(I,2)-DIST(J,2)
			AUXRADIUS=RADIUS(I)+RADIUS(J)
			IF (DISTANCE(AUXX,AUXY).LE.AUXRADIUS) THEN
				WRITE(*,*) '¡¡¡ATENCIÓN: Los planetas han chocado!!!'
				STOP
			END IF
300		CONTINUE
		ENDDO
	ENDDO

! Calling function FORCE in t+h
	CALL FORCE(F,DIST)

! New w en t+h
	DO I=1,N
		DO J=1,2
			W(I,J)=W(I,J)+TH*F(I,J)
		ENDDO
	ENDDO

! SPEED in t+h
	DO I=1,N
		DO J=1,2
			SPEED(I,J)=W(I,J)-TH*F(I,J)/2
		ENDDO
		WRITE(333,*) SPEED(I,1), SPEED(I,2), I
	ENDDO
	WRITE(333,*)

! Time increment
	T=T+TH

	CALL MOMANG(DIST,SPEED,L)

	ENDDO
! End of big loop.

	CLOSE(111)
	CLOSE(222)
	CLOSE(666)
	CLOSE(555)
	CLOSE(444)
	CLOSE(333)

100	STOP
	END PROGRAM

!***********************************************************************

! Subroutine for FORCE.
      SUBROUTINE FORCE(F,DIST)
      
      PARAMETER (N=10)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DIST(10,2),F(10,2)
      REAL*8 MASS(10)
      REAL*8 DISTANCE,DX,DY
	COMMON MASS

      DX=0.00
      DY=0.00
	DO I=1,N
		DO J=1,2
			F(I,J)=0.0D0
		ENDDO
	ENDDO

	DO I=1,N
		DO J=1,N
			IF (I.EQ.J) GOTO 100
			AUXX=DIST(I,1)-DIST(J,1)
			AUXY=DIST(I,2)-DIST(J,2)
			AUXDIST=DISTANCE(AUXX,AUXY)
			F(I,1)=F(I,1)-MASS(J)*AUXX/(AUXDIST**3)
			F(I,2)=F(I,2)-MASS(J)*AUXY/(AUXDIST**3)
!			WRITE(444,*) F(I,1), F(I,2), I, J
			WRITE(444,*) 'Modulo distance=', AUXDIST
100			CONTINUE
		ENDDO
		WRITE(444,*)
		WRITE(444,*) 'Planet',I
		WRITE(444,*) 'FORCE=', 'Fx=', F(I,1), 'Fy=', F(I,2)
		WRITE(444,*) 'Modulo FORCE=', DISTANCE(F(I,1),F(I,2))
		WRITE(444,*)
	ENDDO
	
	RETURN
	END SUBROUTINE


! function that calculate the distance betweeen planets.
	FUNCTION DISTANCE(DX,DY)

	REAL*8 DISTANCE, DX, DY

	DISTANCE=DSQRT(DX**2+DY**2)
	
	RETURN
	END FUNCTION


! SUBROUTINE for angular momentum.
	SUBROUTINE MOMANG(DIST,SPEED,L)

	IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 DIST(10,2), SPEED(10,2), MASS(10), L
	PARAMETER (N=10)
	COMMON MASS

	AUXL=0.0D0
	L=0D+0
	DO I=1,N
		AUXV=DISTANCE(SPEED(I,1),SPEED(I,2))
		AUXD=DISTANCE(DIST(I,1),DIST(I,2))
		L=L+MASS(I)*AUXV*AUXD
	ENDDO
	
	WRITE(222,*) L
	
	RETURN
	END SUBROUTINE
      
! ------------------------------------------------------
! 		This is the end
! ------------------------------------------------------
