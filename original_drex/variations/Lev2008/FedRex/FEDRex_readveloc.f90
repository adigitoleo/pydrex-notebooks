!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine READVELOC, read a velocity input file for a given time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE readveloc (fname_base, step)

   USE comvar
   USE GENERAL_VARS

   IMPLICIT NONE

   CHARACTER (80) fname_base
   CHARACTER (100) :: fname_long
   	! baseline name of the imput file containing X1, X3, Ui and Dij; 
	! Format of name: baseline.step.dat, where 'step' is always 5 characters)
	! Format of file: X1, X3, U1, U2, U3, D11 D13 D23 D31 (assume D33=-D11 in 2D)
   CHARACTER(10) step_string
   INTEGER :: step, i,i1,i2,i3 ! counters

!!! Initialization of grid, velocity field and velocity gradient
   X1 = 0d0 ; X2 = 0d0; X3 = 0d0 ; Ui = 0d0 ; Dij = 0d0

   write(step_string,'(i5.5)') step

   !!fname_long = fname_base // '.' // trim(step_string) // '.dat'
    do i = len(fname_base), 1, -1
        if (fname_base(i:i).ne.' ') go to 500
    end do
500  continue

   fname_long(1:i)=fname_base(1:i)
   fname_long(i+1:i+1)='.'
   fname_long(i+2:i+6)=step_string(1:5)
   fname_long(i+7:i+10)='.dat'
   print *, fname_long

   OPEN(8,file=fname_long)

   IF (dimensions .EQ. 2) THEN
   	i2=1;
	   DO i1 = 1, nx1
	      DO i3 = 1, nx3
	         READ(8,*) X1(i1),X3(i3),Ui(1,i1,i2,i3),Ui(2,i1,i2,i3),Ui(3,i1,i2,i3), &
        	     Dij(1,1,i1,i2,i3),Dij(1,3,i1,i2,i3),Dij(2,3,i1,i2,i3),Dij(3,1,i1,i2,i3)
        	 Dij(3,3,i1,i2,i3) = -Dij(1,1,i1,i2,i3)-Dij(2,2,i1,i2,i3)
	      END DO
	   END DO
	   CLOSE(8)
   ELSE
   	print *, 'Trying to work in a 3D mode, using file ', fname_long
	 DO i1 = 1, nx1
	      DO i2 = 1, nx2
	     	DO i3 = 1, nx3
	       	  READ(8,*) X1(i1),X2(i2),X3(i3),Ui(1,i1,i2,i3),Ui(2,i1,i2,i3),Ui(3,i1,i2,i3), &
     	  	     Dij(1,1,i1,i2,i3),Dij(1,2,i1,i2,i3),Dij(1,3,i1,i2,i3), &
     		     Dij(2,1,i1,i2,i3),Dij(2,2,i1,i2,i3),Dij(2,3,i1,i2,i3), &
     		     Dij(3,1,i1,i2,i3),Dij(3,2,i1,i2,i3)
        	 Dij(3,3,i1,i2,i3) = -Dij(1,1,i1,i2,i3)-Dij(2,2,i1,i2,i3)
	      END DO
	   END DO
         END DO
	CLOSE (8)	   
   END IF

  RETURN

  7000 FORMAT(16(1pe14.6))
  6000 FORMAT(a200)

END SUBROUTINE readveloc

