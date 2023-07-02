!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine PATHLINE - Calculation of tracers pathlines                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Changed: Einat Lev, July 2007 
!! This subroutine only needs now to calculate the velocity gradient and timestep
!! the current location of the particle, and to calculate the location to which
!! it will be advected by the end of the timestep.
!!
!! Input arguments: 
!! X1s, X3s: current location of the particle
!! step_counter: the current timestep
!! 
!! Output parameters:
!! X1next, X3next: the location at the end of the timestep
!!
!! Global variables set by this subroutine:
!! strem_Lij, stream_dt, stream_e
!!
!!

   SUBROUTINE pathline(X1s,X2s,X3s, step_counter, X1next, X2next,X3next,box_exit_flag) 

   USE comvar 
   USE GENERAL_VARS

   IMPLICIT NONE

   INTEGER :: step_counter,box_exit_flag
   ! number of steps to construct the streamline

   INTEGER :: i1s,i2s,i3s,i1,i2,i3,i1r,i2r,i3r
   ! initial, intermediate and final indices of the closest upper-left grid point

   DOUBLE PRECISION :: ds,dt
   ! elementary displacement along the streamline and corresponding time step

   DOUBLE PRECISION :: X1s,X2s,X3s,X1i,X2i,X3i, X1next, X2next,X3next
   ! initial and intermediate position of the tracer

   DOUBLE PRECISION :: U1i,U2i,U3i
   ! intermediate velocity components on the streamline

   DOUBLE PRECISION :: kx1,kx2,kx3,kx4
   DOUBLE PRECISION :: ky1,ky2,ky3,ky4
   DOUBLE PRECISION :: kz1,kz2,kz3,kz4
   ! Runge-Kutta intermediate for streamlines calculation 

   DOUBLE PRECISION :: max_strain
   ! value of the dimensionless strain for which any former
   ! LPO has been erased. From Kaminski and Ribe 2002, 10
   ! is a safe value. Should be tested for a given flow though.

!! initial indices and position of the tracers

   max_strain = 0d0 

      X1i = X1s ; X2i = X2s ; X3i = X3s
      i1 = nx1-1 ; i3 = nx3-1 ; i2 = nx2-1
	
!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X2(i2) .GT. X2i .AND. i2 .GT. 1) ; i2 = i2-1 ; ENDDO
      DO WHILE (X2(i2+1) .LT. X2i .AND. i2 .LT. nx2-1) ; i2 = i2+1 ; ENDDO
      
      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO
    
      CALL velocitycalc(X1i,X2i,X3i,U1i,U2i,U3i,i1,i2,i3) !! X1i, X3i - inputs, U1i,U3i, i1, i3 outputs
      CALL gradientcalc(X1i,X2i,X3i,i1,i2,i3)

      IF (dimensions .EQ. 3) THEN
      	ds = 0.25d0*MIN(ABS(X1(i1+1)-X1(i1)),MIN(ABS(X2(i2+1)-X2(i2)),ABS(X3(i3+1)-X3(i3))))
      ELSE
      	ds = 0.25d0*MIN(ABS(X1(i1+1)-X1(i1)),ABS(X3(i3+1)-X3(i3)))
      END IF
      dt = ds/SQRT(U1i**2+U2i**2+U3i**2)

! record of the local velocity gradient tensor and time spent at that point
! "l" is calculated in "gradientcalc" called above

      stream_Lij(:,:,step_counter) = l(:,:) ; stream_dt(step_counter) = dt 
      stream_e(step_counter) = epsnot
      max_strain = max_strain + dt*epsnot

      kx1 = U1i*dt ; ky1 = U2i*dt ; kz1 = U3i*dt

      X1i = X1s + 0.5d0*kx1
      X2i = X2s + 0.5d0*ky1
      X3i = X3s + 0.5d0*kz1
            
!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO
      
      DO WHILE (X2(i2) .GT. X2i .AND. i2 .GT. 1) ; i2 = i2-1 ; ENDDO
      DO WHILE (X2(i2+1) .LT. X2i .AND. i2 .LT. nx2-1) ; i2 = i2+1 ; ENDDO
      
      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X2i,X3i,U1i,U2i,U3i,i1,i2,i3) !! X1i, X3i - inputs, U1i,U3i, i1, i3 outputs
      
      kx2 = U1i*dt ; ky2 = U2i*dt ; kz2 = U3i*dt
      X1i = X1s + 0.5d0*kx2
      X2i = X2s + 0.5d0*ky2
      X3i = X3s + 0.5d0*kz2

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X2i,X3i,U1i,U2i,U3i,i1,i2,i3) !! X1i, X3i - inputs, U1i,U3i, i1, i3 outputs

      kx3 = U1i*dt ; ky3 = U2i*dt ; kz3 = U3i*dt

      X1i = X1s + kx3
      X2i = X2s + ky3
      X3i = X3s + kz3

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X2(i2) .GT. X2i .AND. i2 .GT. 1) ; i2 = i2-1 ; ENDDO
      DO WHILE (X2(i2+1) .LT. X2i .AND. i2 .LT. nx2-1) ; i2 = i2+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X2i,X3i,U1i,U2i,U3i,i1,i2,i3) !! X1i, X3i - inputs, U1i,U3i, i1, i3 outputs

      kx4 = U1i*dt ; ky4 = U2i*dt ; kz4 = U3i*dt

      X1s = X1s + (kx1/2d0+kx2+kx3+kx4/2d0)/3d0
      X2s = X2s + (ky1/2d0+ky2+ky3+ky4/2d0)/3d0
      X3s = X3s + (kz1/2d0+kz2+kz3+kz4/2d0)/3d0
      
      box_exit_flag = 0;
      IF ( (X1s .GT. X1(nx1)) .OR. (X2s .GT. X2(nx2)) .OR. (X3s .GT. X3(nx3)) &
        .OR. (X1s .LT. X1(1)) .OR. (X2s .LT. X2(1)) .OR. (X3s .LT. X3(1))) THEN
      	 	box_exit_flag = 1;
      	return;
      END IF
        
!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1s .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1s .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X2(i2) .GT. X2i .AND. i2 .GT. 1) ; i2 = i2-1 ; ENDDO
      DO WHILE (X2(i2+1) .LT. X2i .AND. i2 .LT. nx2-1) ; i2 = i2+1 ; ENDDO
      
      DO WHILE (X3(i3) .GT. X3s .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3s .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      i1r = i1 ; i2r = i2 ; i3r = i3	

	X1next = X1s
	X2next = X2s
	X3next = X3s
	
   RETURN

   END SUBROUTINE pathline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
