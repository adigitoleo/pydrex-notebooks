
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine VELOCITYCALC, calculation of velocity at a given point      !!!
!!! by interpolation method given in Numerical Recipies, Press et al., p96 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE velocitycalc(X1s,X2s,X3s,U1s,U2s,U3s,i1s,i2s,i3s)

   USE comvar
   USE GENERAL_VARS

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s,X3s
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: U1s,U2s,U3s
   ! interpolated velocity at the point

   INTEGER :: i1s,i2s,i3s
   ! indices of the UP-LEFT grid point closest to the extrapolation point

   DOUBLE PRECISION :: y1,y2,y3,y4,y5,y6,y7,y8

   ! dummies for interpolation
   IF (dimensions .EQ. 3) THEN
	y1 = Ui(1,i1s,i2s,i3s) ; 	y2 = Ui(1,i1s+1,i2s,i3s)
	y3 = Ui(1,i1s+1,i2s+1,i3s) ;	y4 = Ui(1,i1s,i2s+1,i3s)
	y5 = Ui(1,i1s,i2s,i3s+1) ; 	y6 = Ui(1,i1s+1,i2s,i3s+1)
	y7 = Ui(1,i1s+1,i2s+1,i3s+1);y8 = Ui(1,i1s,i2s+1,i3s+1)
	
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,U1s)
	
	y1 = Ui(2,i1s,i2s,i3s) ; 	y2 = Ui(2,i1s+1,i2s,i3s)
	y3 = Ui(2,i1s+1,i2s+1,i3s) ; y4 = Ui(2,i1s,i2s+1,i3s) 
	y5 = Ui(2,i1s,i2s,i3s+1) ; 	y6 = Ui(2,i1s+1,i2s,i3s+1)
	y7 = Ui(2,i1s+1,i2s+1,i3s+1);y8 = Ui(2,i1s,i2s+1,i3s+1)
	
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,U2s)
	
	y1 = Ui(3,i1s,i2s,i3s) ; 	y2 = Ui(3,i1s+1,i2s,i3s)
	y3 = Ui(3,i1s+1,i2s+1,i3s) ; y4 = Ui(3,i1s,i2s+1,i3s) 
	y5 = Ui(3,i1s,i2s,i3s+1) ; 	y6 = Ui(3,i1s+1,i2s,i3s+1)
	y7 = Ui(3,i1s+1,i2s+1,i3s+1);y8 = Ui(3,i1s,i2s+1,i3s+1)
	
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,U3s)
   ELSE	!!!2D
   	y3=0d0; y4=0d0; y7=0d0; y8=0d0; 
	i2s = 1;
	
	y1 = Ui(1,i1s,i2s,i3s) ; 	y2 = Ui(1,i1s+1,i2s,i3s)
	y5 = Ui(1,i1s,i2s,i3s+1) ; 	y6 = Ui(1,i1s+1,i2s,i3s+1)
 	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,U1s)
	
	y1 = Ui(2,i1s,i2s,i3s) ; 	y2 = Ui(2,i1s+1,i2s,i3s)
	y5 = Ui(2,i1s,i2s,i3s+1) ; 	y6 = Ui(2,i1s+1,i2s,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,U2s)
 	
	y1 = Ui(3,i1s,i2s,i3s) ; 	y2 = Ui(3,i1s+1,i2s,i3s)
	y5 = Ui(3,i1s,i2s,i3s+1) ; 	y6 = Ui(3,i1s+1,i2s,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,U3s)
   ENDIF
 
   RETURN

   END SUBROUTINE velocitycalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine interp, bilinear interpolation of a function                !!!
!!! following the method given in Numerical Recipies, Press et al., p96    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,res)

   USE comvar
   USE GENERAL_VARS

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s,X3s
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: res
   ! result of the interpolation at the point

   INTEGER :: i1s,i2s,i3s
   ! indices of the UP-LEFT-FRONT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4,y5,y6,y7,y8,tt,uu,vv
   ! dummies for interpolation (numerical recipies, p96)

   tt = (X1s-X1(i1s))/(X1(i1s+1)-X1(i1s))
   IF (dimensions .EQ. 3) THEN
  	uu = (X2s-X2(i2s))/(X2(i2s+1)-X2(i2s))
   ELSE
    	uu = 0
   END IF
   
   vv = (X3s-X3(i3s))/(X3(i3s+1)-X3(i3s))  
   
   res = (1-vv)*((1d0-tt)*(1d0-uu)*y1+tt*(1d0-uu)*y2+tt*uu*y3+(1d0-tt)*uu*y4) + 		&
	 vv    *((1d0-tt)*(1d0-uu)*y5+tt*(1d0-uu)*y6+tt*uu*y7+(1d0-tt)*uu*y8)
	 
    RETURN

   END SUBROUTINE interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine GRADIENTCALC, interpolation of velocity gradient tensor     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE gradientcalc(X1s,X2s,X3s,i1s,i2s,i3s)

   USE comvar
   USE GENERAL_VARS

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s,X3s
   ! x and y coordinates of the point on the streamline

   INTEGER :: i1s,i2s,i3s
   ! indices of the UPPER-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4,y5,y6,y7,y8
   ! dummies for interpolation

   DOUBLE PRECISION, DIMENSION(3,3) :: ee
   ! dummy for reference strain rate calculation

   INTEGER :: nrot
   ! number of rotations for the Jacobi
   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects
   ! eigen values and vectors in jacobi

   l = 0d0 ; e = 0d0

   IF (dimensions .EQ. 3) THEN
	
	y1 = Dij(1,1,i1s,i2s,i3s) ; 	    y2 = Dij(1,1,i1s+1,i2s,i3s)
	y3 = Dij(1,1,i1s+1,i2s+1,i3s);   y4 = Dij(1,1,i1s,i2s+1,i3s)
	y5 = Dij(1,1,i1s,i2s,i3s+1) ;    y6 = Dij(1,1,i1s+1,i2s,i3s+1)
	y7 = Dij(1,1,i1s+1,i2s+1,i3s+1); y8 = Dij(1,1,i1s,i2s+1,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(1,1))
	
	y1 = Dij(1,2,i1s,i2s,i3s) ; y2 = Dij(1,2,i1s+1,i2s,i3s)
	y3 = Dij(1,2,i1s+1,i2s+1,i3s) ; y4 = Dij(1,2,i1s,i2s+1,i3s)
	y5 = Dij(1,2,i1s,i2s,i3s+1) ;    y6 = Dij(1,2,i1s+1,i2s,i3s+1)
	y7 = Dij(1,2,i1s+1,i2s+1,i3s+1); y8 = Dij(1,2,i1s,i2s+1,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(1,2))   
	
	y1 = Dij(1,3,i1s,i2s,i3s) ; y2 = Dij(1,3,i1s+1,i2s,i3s)
	y3 = Dij(1,3,i1s+1,i2s+1,i3s) ; y4 = Dij(1,3,i1s,i2s+1,i3s)
	y5 = Dij(1,3,i1s,i2s,i3s+1) ;    y6 = Dij(1,3,i1s+1,i2s,i3s+1)
	y7 = Dij(1,3,i1s+1,i2s+1,i3s+1); y8 = Dij(1,3,i1s,i2s+1,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(1,3))
	
	y1 = Dij(2,1,i1s,i2s,i3s) ; y2 = Dij(2,1,i1s+1,i2s,i3s)
	y3 = Dij(2,1,i1s+1,i2s+1,i3s) ; y4 = Dij(2,1,i1s,i2s+1,i3s)
	y5 = Dij(2,1,i1s,i2s,i3s+1) ;    y6 = Dij(2,1,i1s+1,i2s,i3s+1)
	y7 = Dij(2,1,i1s+1,i2s+1,i3s+1); y8 = Dij(2,1,i1s,i2s+1,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(2,1))
	
	y1 = Dij(2,2,i1s,i2s,i3s) ; y2 = Dij(2,2,i1s+1,i2s,i3s)
	y3 = Dij(2,2,i1s+1,i2s+1,i3s) ; y4 = Dij(2,2,i1s,i2s+1,i3s)
	y5 = Dij(2,2,i1s,i2s,i3s+1) ;    y6 = Dij(2,2,i1s+1,i2s,i3s+1)
	y7 = Dij(2,2,i1s+1,i2s+1,i3s+1); y8 = Dij(2,2,i1s,i2s+1,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(2,2))
	
	y1 = Dij(3,1,i1s,i2s,i3s) ; y2 = Dij(3,1,i1s+1,i2s,i3s)
	y3 = Dij(3,1,i1s+1,i2s+1,i3s) ; y4 = Dij(3,1,i1s,i2s+1,i3s)
	y5 = Dij(3,1,i1s,i2s,i3s+1) ;    y6 = Dij(3,1,i1s+1,i2s,i3s+1)
	y7 = Dij(3,1,i1s+1,i2s+1,i3s+1); y8 = Dij(3,1,i1s,i2s+1,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(3,1))
	
	y1 = Dij(2,3,i1s,i2s,i3s) ; y2 = Dij(2,3,i1s+1,i2s,i3s)
	y3 = Dij(2,3,i1s+1,i2s+1,i3s) ; y4 = Dij(2,3,i1s,i2s+1,i3s)
	y5 = Dij(2,3,i1s,i2s,i3s+1) ;    y6 = Dij(2,3,i1s+1,i2s,i3s+1)
	y7 = Dij(2,3,i1s+1,i2s+1,i3s+1); y8 = Dij(2,3,i1s,i2s+1,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(2,3))
	
	y1 = Dij(3,2,i1s,i2s,i3s) ; y2 = Dij(3,2,i1s+1,i2s,i3s)
	y3 = Dij(3,2,i1s+1,i2s+1,i3s) ; y4 = Dij(3,2,i1s,i2s+1,i3s)
	y5 = Dij(3,2,i1s,i2s,i3s+1) ;    y6 = Dij(3,2,i1s+1,i2s,i3s+1)
	y7 = Dij(3,2,i1s+1,i2s+1,i3s+1); y8 = Dij(3,2,i1s,i2s+1,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(3,2))   
	
	l(3,3) = -l(1,1)-l(2,2)
	
   ELSE !2D case
   	i2s = 1
   	y3 = 0d0; y4 = 0d0; y7 = 0d0 ; y8 = 0d0;
   	 
   	y1 = Dij(1,1,i1s,i2s,i3s) ; 	 y2 = Dij(1,1,i1s+1,i2s,i3s)
   	y5 = Dij(1,1,i1s,i2s,i3s+1) ;    y6 = Dij(1,1,i1s+1,i2s,i3s+1)
   	
   	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(1,1))

	y1 = Dij(1,2,i1s,i2s,i3s) ; y2 = Dij(1,2,i1s+1,i2s,i3s)
	y5 = Dij(1,2,i1s,i2s,i3s+1) ;y6 = Dij(1,2,i1s+1,i2s,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(1,2))   
	
	y1 = Dij(1,3,i1s,i2s,i3s) ; y2 = Dij(1,3,i1s+1,i2s,i3s)
	y5 = Dij(1,3,i1s,i2s,i3s+1) ;    y6 = Dij(1,3,i1s+1,i2s,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(1,3))
	
	y1 = Dij(2,1,i1s,i2s,i3s) ; y2 = Dij(2,1,i1s+1,i2s,i3s)
	y5 = Dij(2,1,i1s,i2s,i3s+1) ;    y6 = Dij(2,1,i1s+1,i2s,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(2,1))
	
	y1 = Dij(2,2,i1s,i2s,i3s) ; y2 = Dij(2,2,i1s+1,i2s,i3s)
	y5 = Dij(2,2,i1s,i2s,i3s+1) ;    y6 = Dij(2,2,i1s+1,i2s,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(2,2))
	
	y1 = Dij(3,1,i1s,i2s,i3s) ; y2 = Dij(3,1,i1s+1,i2s,i3s)
	y5 = Dij(3,1,i1s,i2s,i3s+1) ;    y6 = Dij(3,1,i1s+1,i2s,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(3,1))
	
	y1 = Dij(2,3,i1s,i2s,i3s) ; y2 = Dij(2,3,i1s+1,i2s,i3s)
	y5 = Dij(2,3,i1s,i2s,i3s+1) ;    y6 = Dij(2,3,i1s+1,i2s,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(2,3))
	
	y1 = Dij(3,2,i1s,i2s,i3s) ; y2 = Dij(3,2,i1s+1,i2s,i3s)
	y5 = Dij(3,2,i1s,i2s,i3s+1) ;    y6 = Dij(3,2,i1s+1,i2s,i3s+1)
	CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(3,2)) 
	
	l(3,3) = -l(1,1)-l(2,2)
  
   END IF

!!! strain rate tensorc- symmetic
   e(1,1) = l(1,1) ; e(3,3) = l(3,3) ; e(2,2) = l(2,2)
   e(3,1) = (l(1,3)+l(3,1))/2d0 ; e(1,3) = e(3,1)
   e(1,2) = (l(2,1)+l(1,2))/2d0 ; e(2,1) = e(1,2)
   e(2,3) = (l(3,2)+l(2,3))/2d0 ; e(3,2) = e(2,3)

!! reference strain rate
   ee = e
   CALL JACOBI(ee,3,3,evals,evects,nrot)
   epsnot = MAXVAL(ABS(evals))

   RETURN

   END SUBROUTINE gradientcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine check_veloc - matlab plot of velocity field                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE check_veloc

   USE comvar
   USE GENERAL_VARS

   IMPLICIT NONE
   
   INTEGER :: i1,i3 ! loop counters

   OPEN(10,file='veloc.m',status='replace')

   WRITE(10,*)
   !!WRITE(10,*) "clear"
   WRITE(10,*) "figure"

   WRITE(10,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(10,"(81(1PE12.4))") Ui(1,:,1,i3)
   END DO
   WRITE(10,*) "];"
   
   WRITE(10,*) "Uz=["
   DO i3 = 1,nx3
      WRITE(10,"(81(1PE12.4))") Ui(3,:,1,i3)
   END DO
   WRITE(10,*) "];"

   WRITE(10,*) "X=["
   WRITE(10,"(81(1PE12.4))") X1
   WRITE(10,*) "];"

   WRITE(10,*) "Z=["
   WRITE(10,"(81(1PE12.4))") X3
   WRITE(10,*) "];"

   WRITE(10,*) "quiver(X,Z,Ux,Uz)"
   WRITE(10,*) "axis([-4 4 0 3])"
   WRITE(10,*) "xlabel('X1, dimensionless distance from the ridge axis)')"
   WRITE(10,*) "ylabel('X3, dimensionless depth')"
   WRITE(10,*) "title('Velocity field')"

   WRITE(10,*) "axis square"
   WRITE(10,*) "axis ij"

   CLOSE(10)

   RETURN

   END SUBROUTINE check_veloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine pipar - Calculates GOL parameter at grid point              !!!
!!! Function name comes from PI Parameter, defining the grain lag.	   !!!
!!! This subroutine goes over all the velocity grid points and calculates  !!!
!!! the PI parameter at each grid point.				   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE pipar(i1,i2,i3,X10,X20, X30)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i1,i2, i3
   ! coordinates of the grid point

   DOUBLE PRECISION :: X10,X20, X30,X1i,X2i, X3i
   ! initial and intermediate position of the point

   DOUBLE PRECISION :: U1i,U2i, U3i
   ! local velocity vector

   DOUBLE PRECISION, DIMENSION(3) :: veloc,isa
   ! local velocity vector and orientation of infinite strain axis

   DOUBLE PRECISION :: dt
   ! time step used to calculate GOL

   DOUBLE PRECISION :: thetaISA
   ! angle between ISA and flow direction

!!! RMQ if GOL not defined then GOL = -1d0 - Max value of GOL set to 10

   GOL(i1,i2,i3) = 1d1


   veloc = 0d0 ; isa = 0d0


 
   CALL velocitycalc(X10,X20,X30,U1i,U2i,U3i,i1,i2,i3)

   
   dt = 1d-2/SQRT(U1i**2+U2i**2+U3i**2)

!!! previous point on the streamline

   X1i = X10-dt*U1i
   X2i = X20-dt*U2i
   X3i = X30-dt*U3i

   CALL velocitycalc(X1i,X2i,X3i,U1i,U2i,U3i,i1,i2,i3)
   veloc(1) = U1i ; veloc(2) = U2i ; veloc(3) = U3i ; 
   veloc = veloc/SQRT(SUM(veloc**2)) 

   CALL gradientcalc(X1i,X2i,X3i,i1,i2,i3)
   l = l/epsnot

!!! calculation of the ISA
   CALL isacalc(isa,i1,i2,i3)
   IF (SUM(isa) .EQ. -3d0) isa=veloc

!!! angle between ISA and flow direction
   thetaISA = ACOS(SUM(veloc*isa))

!!! next point on the streamline

   X1i = X10+dt*U1i
   X2i = X20+dt*U2i
   X3i = X30+dt*U3i

   CALL velocitycalc(X1i,X2i, X3i,U1i,U2i, U3i,i1,i2,i3)
   veloc(1) = U1i ; veloc(2) = U2i ; veloc(3) = U3i ; 
   veloc = veloc/SQRT(SUM(veloc**2))

  CALL gradientcalc(X1i,X2i,X3i,i1,i2,i3)
   l = l/epsnot

!!! calculation of the ISA
   CALL isacalc(isa,i1,i2,i3)
   IF (SUM(isa) .EQ. -3d0) THEN
   	isa=veloc
   END IF

!!! Pi parameter
   CALL gradientcalc(X10,X20,X30,i1,i2,i3)
   thetaISA = ABS(thetaISA-ACOS(SUM(veloc*isa)))
   IF (thetaISA .GT. ACOS(-1d0)) thetaISA = thetaISA-ACOS(-1d0)

   GOL(i1,i2,i3) = MIN(GOL(i1,i2,i3),thetaISA/2d0/dt/epsnot)

   RETURN

   END SUBROUTINE pipar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine isacalc - Calculates ISA orienation at grid point           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE isacalc(isa,i1,i2,i3)

   USE comvar
   USE GENERAL_VARS

   IMPLICIT NONE

   INTEGER :: i1,i2,i3
   ! indices of the point of calculation

   DOUBLE PRECISION, DIMENSION(3) :: isa
   ! orientation local of infinite strain axis

   DOUBLE PRECISION, DIMENSION(3,3) :: F,U
   ! deformation gradient tensor and left strech tensor
   
   INTEGER :: nrot
   ! number of rotations for the Jacobi

   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects
   ! eigen values and vectors in jacobi

!!! Rmq: if GOL is not defined, then isa=0 and GOL(i1,i3)=-1d0

!!! calculation of the ISE orientation using Sylvester's formula
   
!!! Limit deformation gradient tensor for infinite time
   CALL eigen(F)
  
   IF (SUM(F) .EQ. -9d0) THEN
      isa = -1d0 
!!! 2. formation of the left-stretch tensor U = FFt
   ELSE IF (SUM(ABS(F)) .EQ. 0) THEN
      isa = 0d0 ; !!GOL(i1,i3)=-1d0
   ELSE
      U = MATMUL(F,TRANSPOSE(F)) 

!!! 3. eigen-values and eigen-vectors of U

      CALL JACOBI(U,3,3,evals,evects,nrot)

      IF (evals(1) .EQ. MAXVAL(evals)) THEN
         isa = evects(:,1)
         ELSE IF(evals(2) .EQ. MAXVAL(evals)) THEN
            isa = evects(:,2)
            ELSE
            isa = evects(:,3)
      END IF

   END IF

   RETURN

   END SUBROUTINE isacalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
