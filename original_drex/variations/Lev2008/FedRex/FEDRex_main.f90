!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FedRex - Version 1							   !!!
!!! MIT April 2008							   !!!
!!! This is a modification of the program DRex. The changes include:       !!!
!!! 1) time-dependent input velocity fields				   !!!
!!! 2) flag for non-random initial orientation of the grains	   	   !!!
!!! 3) each particle advected in the flow can have a specific composition  !!!
!!!	(percentage of olivine) and initial position (in X,Y coords)	   !!!
!!! 4) Output files now include also the time evolution of: A-axis, 	   !!!
!!!  	FSE, particle position. Also included is the final Euler angles of !!!
!!!     the grains.							   !!!
!!! 5) Analysis of computation time of various steps of the code.	   !!!
!!! 6) Fully three-dimentional flow fields now supported		   !!!
!!!									   !!!
!!! (future development may include prescribed initial grain orientation,  !!!
!!! to allow for existing LPO, e.g. LPO measured in field samples).	   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DRex  - version 2                                                      !!!
!!! Paris 07/21/2004 - Calculates the evolutions of texture an olivine     !!!
!!! aggregate in a ridge+shear flow (flow from Neil Ribe).                 !!!
!!! LPO evolves by dynamic recrystallization - Calculated using D-Rex      !!!
!!! Equivalent transverse isotropic media calculated using the theory of   !!!
!!! Jules Browaeys. A couple of bugs fixed thanks to Teresa Lassak         !!!
!!! Pathline subroutine modified to allow calculation in flows presenting  !!!
!!! stagnation points (additional parameter max_strain added). Relevant    !!!
!!! information on the streamline (Lij,strain rate, time) are recorded     !!!
!!! once for all during backward calculation of the streamline.            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Module of common variables                                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   MODULE comvar

   CHARACTER(80) :: fname
   ! name of the imput file containing X1, X3, Ui and Dij

   CHARACTER(80) :: EulerAnglesFileName
   ! name of the imput file containing initial Euler angles

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xe1,xe2,xe3
   ! matrixes of initial random eulerian angles



!!! Grid points - Velocity field - Velocity gradient tensor
   INTEGER :: nx1, nx2, nx3 !! grid size in X1 and X3 direction
   INTEGER :: n_particles !! number of particles in the model

   INTEGER :: timestep_first, timestep_last, timestep_jump
   !! which timesteps (and acompanying velocity files) should be used

   !! LPO calculates only every stepx1 and stepx3 point
   INTEGER :: dimensions, i1first,i3first,i1last,i3last
   !! indices of the first and last points at which LPO is calculated
   

   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Fij
   DOUBLE PRECISION, DIMENSION(1000,3,3) :: Fij_steps
   ! finite strain tensor - F11(particle)=Fij(1,1,particle)


   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: grout
   ! dummy used to activate graphic output only for calculated points

!!! LPO calculation

   DOUBLE PRECISION, DIMENSION(3,3) :: l,e
   ! velocity gradient tensor and strain rate tensor

   DOUBLE PRECISION :: epsnot
   ! reference strain rate

   DOUBLE PRECISION, DIMENSION(1000) :: pathline_x1
   DOUBLE PRECISION, DIMENSION(1000) :: pathline_x2
   DOUBLE PRECISION, DIMENSION(1000) :: pathline_x3
   DOUBLE PRECISION, DIMENSION(3,3,1000) :: stream_Lij
   DOUBLE PRECISION, DIMENSION(1000) :: stream_dt
   DOUBLE PRECISION, DIMENSION(1000) :: stream_e
   ! values of Lij, time spent and strain rate at each location on the streamline

!!! Dynamic recrystallization
   
   INTEGER :: size3, size ! size = size3^3
   ! number of points in the (metric) Eulerian space

   DOUBLE PRECISION :: lambda, Mob, chi, random_flag
   ! lambda = nucleation parameter
   ! Mob = grain mobility
   ! threshold volume fraction for activation of grain boundary sliding

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Xol
   ! fraction of olivine in the aggregate

 
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: fse_particles
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: odf,dotodf,odfi
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: odf_particles
   ! volume fraction of the olivine grains
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: odf_ens,dotodf_ens,odfi_ens
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: odf_steps
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: odf_ens_particles
   ! volume fraction of the enstatite grains

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rt,rt_ens
   ! dislocations density for olivine and enstatite

   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: acs,dotacs,acsi,acs0
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: acs_ens,dotacs_ens,acsi_ens
   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: acs_step
   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: acs_particles,acs_ens_particles
   !! matrix of direction cosine


!!! Output

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_fse
   ! orientation of the long axis of the FSE on the (X1,X3) plan
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fse_axis
   ! orientation of the long axis of the FSE on the (X1,X3) plan
   

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ln_fse
   ! finite strain = ln(longaxis/shortaxis)

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_a
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: hex_axis
   ! average orientation of a-axis on (X1,X3) plan

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: perc_a
   ! percentage of S wave anisotropy

!!EINAT MAY 11
   DOUBLE PRECISION, DIMENSION(1000) :: perc_a_steps
   ! percentage of S wave anisotropy

   DOUBLE PRECISION, DIMENSION(1000) :: ln_fse_steps
   ! percentage of S wave anisotropy

   DOUBLE PRECISION, DIMENSION(1000) :: phi_fse_steps
   ! percentage of S wave anisotropy

   DOUBLE PRECISION, DIMENSION(1000) :: phi_a_steps
   ! percentage of S wave anisotropy

!!! Elastic tensor 
   DOUBLE PRECISION, DIMENSION(6,6) :: Sav
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Sav_particles

   ! stiffness matrix
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: Cav
   !Cijkl tensor at point of calculation

!   END MODULE comvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   DOUBLE PRECISION, DIMENSION(3,3,3) :: alt ! \epsilon_{ijk} tensor
   DOUBLE PRECISION, DIMENSION(3,3) :: del ! \delta_{ij} tensor
   INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
   INTEGER, DIMENSION(6) :: l1,l2 ! tensor of indices to form Sij from Cijkl

   DOUBLE PRECISION, DIMENSION(4) :: tau
   ! RSS for the 4 slip systems of olivine (only 3 are activated)
   DOUBLE PRECISION :: tau_ens ! RSS of enstatite slip system
   DOUBLE PRECISION :: stressexp ! stress exponent for olivine and enstatite
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1
   ! grid points coordinates - X
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X2
   ! grid points coordinates - Y
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X3
   ! grid points coordinates - depth

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1p, X2p, X3p, box_drop_status
   !! particle locations	

   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Ui
   ! velocity vector - U1(x1,x2,x3)=Ui(1,x1,x2,x3)

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1gr, X2gr, X3gr
   !! actual coordinates of calculation point when not a grid point

   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: GOL
   ! Grain Orientation Lag i.e. Pi parameter

   DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE :: Dij
   ! velocity gradient tensor - D11(x1,x2,x3)=Dij(1,1,x1,x2,x3)
   DOUBLE PRECISION, DIMENSION(6,6) :: S0,S0_ens
   ! stiffness matrix

END MODULE comvar

MODULE GENERAL_VARS
END MODULE GENERAL_VARS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Main program                                                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   PROGRAM DR_CF 

   USE comvar
   USE GENERAL_VARS
   USE DECMOD

!-------------------------------------------------------------
   IMPLICIT NONE

   DOUBLE PRECISION :: X10,X20,X30,loc_x1,loc_x2,loc_x3
   !! coordinates of the grid point where calculation is started
   DOUBLE PRECISION :: temp_perc, temp_phi
   DOUBLE PRECISION, DIMENSION (3) :: temp_Aaxis

   INTEGER :: i,i1,i2,i3, step_counter, step, steady_state_flag, iter !! loop counters
   INTEGER :: i1a,i2a,i3a !! indices of the grid point at the end of streamline

   INTEGER :: step1,step2!! number of points on a streamline
   INTEGER :: particle !! particle loop counter
   INTEGER :: box_exit_flag

   DOUBLE PRECISION, DIMENSION(3,3) :: LSij
   ! left-strech tensor for FSE calculation

   INTEGER :: nrot
   ! number of rotations for the Jacobi

   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects, scc_rot


   !! Variables for TWB code
   double precision, dimension (3) :: tiaxis
   ! symmetry axis of approximated hexagonal tensor

   double precision :: k,g,favg
   double precision, dimension(9) :: vel ! various velocities from tensor
   double precision, dimension(6,6) :: dummyt
   double precision,dimension(6) :: frac_symm
!-------------------------------------------------------------

!!! initialization of computation cost timers
REAL time_begin, time_end, total_time_pathline, total_time_strain, total_time_voigt, total_time_decsym
total_time_pathline = 0
total_time_strain = 0 
total_time_voigt = 0
total_time_decsym = 0


  CALL init0

   X1gr = 0d0 ; X2gr = 0d0 ; X3gr = 0d0
   phi_fse = 0d0 ; ln_fse = 0d0 ; phi_a = 0d0 ; perc_a = 0d0

   OPEN (17, file='Locations.dat',status='replace')
   OPEN (18, file='FSE.dat',status='replace')
   OPEN (19, file='Aaxis.dat',status='replace')
   OPEN (21, file='EulerAngles.dat',status='replace')
   OPEN (22, file='GOL.dat',status='replace')
   WRITE (17, *), n_particles, n_particles, n_particles
   WRITE (18, *), n_particles, n_particles, n_particles, n_particles
!   WRITE (19, *), n_particles, n_particles, n_particles, n_particles, n_particles, n_particles, n_particles, n_particles, n_particles, n_particles
   WRITE (22, *), nx1, nx2, nx3, 0

!------------------------------------------------------
!!! FSE and LPO calculation
!------------------------------------------------------

   IF (timestep_jump .EQ. 0) THEN
	timestep_jump=1
	steady_state_flag = 1
   END IF
   step_counter = 0;

!!! Loop over timesteps - some may be repeated, according to input file

   DO step=timestep_first,timestep_last,timestep_jump
	DO iter=1,timestep_jump
		step_counter = step_counter + 1
	
		!!! We may allow for a constant velocity field to be re-read by raising the steadt_state_flag. 
		!!! Otherwise, a new velocity field is read
		IF (steady_state_flag .NE. 1) THEN 
			CALL readveloc (fname, step)
		END IF
		

		!!! Calculation of GOL parameter
		IF (dimensions .EQ. 3) THEN
  		  DO i1 = 1, nx1-1
		    DO i2 = 1, nx2-1
		      DO i3 = 1 , nx3-1
		
			!!! Location of the calculation point 
		         X10 = (X1(i1)+X1(i1+1))/2d0 ; 
		         X20 = (X2(i2)+X2(i2+1))/2d0 ; 
		         X30 = (X3(i3)+X3(i3+1))/2d0
		         X1gr(i1) = X10 ; X2gr(i2) = X20 ; X3gr(i3) = X30
		         CALL pipar(i1,i2,i3,X1gr(i1),X2gr(i2),X3gr(i3))
		      END DO
		    END DO
		  END DO
		ELSE
		  DO i1 = 1, nx1-1
		      DO i3 = 1 , nx3-1
			!!! Location of the calculation point 
			i2=1;
		         X10 = (X1(i1)+X1(i1+1))/2d0 ; 
		         X30 = (X3(i3)+X3(i3+1))/2d0 ;
			 X20 = 0;
		         X1gr(i1) = X10 ; X2gr(i2) = X20 ; X3gr(i3) = X30
		         CALL pipar(i1,i2,i3,X1gr(i1),X2gr(i2),X3gr(i3))
		      END DO
		   END DO
		END IF

 
  		!!! Calculation of GOL parameter
		IF (dimensions .EQ. 3) THEN
  		  DO i1 = 1, nx1-1
		    DO i2 = 1, nx2-1
		      DO i3 = 1 , nx3-1
			WRITE (22, "(4(1PE12.4))"), X1(i1), X2(i2), X3(i3), GOL(i1,i2,i3)
		      END DO
		    END DO
		  END DO
		ELSE
		  DO i1 = 1, nx1-1
		      DO i3 = 1 , nx3-1
			WRITE (22, "(4(1PE12.4))"), X1(i1),  0d0, X3(i3), GOL(i1,1,i3)
		      END DO
		   END DO
		END IF



				
		DO particle = 1,n_particles
		       IF (box_drop_status(particle) .EQ. 1) THEN
			       continue;
		       END IF
		       !! Initial location of the particles - using coordinates read from particles input file
	        	X10 = X1p(particle) ; X30 = X3p(particle) ; X20 = X2p(particle)

		       !!! Forward calculation of the pathline for each tracer
	               CALL CPU_TIME ( time_begin)
         		CALL pathline(X10,X20,X30,step_counter,X1p(particle), X2p(particle), X3p(particle),box_exit_flag) !! only advects one 	step along the path
	               CALL CPU_TIME ( time_end)
		       total_time_pathline = total_time_pathline + (time_end-time_begin)

			!!! we would like to ignore particles that exited the model box so we mark them as such for the LPO calculation
         	       IF (box_exit_flag .EQ. 1) THEN
         		       box_drop_status (particle) = 1;
         		       continue;
         	       END IF

	               CALL CPU_TIME ( time_begin)
	        	      CALL strain(step_counter,particle,Fij(:,:,particle))
	               CALL CPU_TIME ( time_end)
		       total_time_strain = total_time_strain + (time_end-time_begin)

		       !!! Left-stretch tensor for FSE calculation
        	       LSij = MATMUL(Fij(:,:,particle),TRANSPOSE(Fij(:,:,particle)))
        	       CALL JACOBI(LSij,3,3,evals,evects,nrot)

 		       !!!Pick up the orientation of the long axis of the FSE
       		       IF (evals(1) .EQ. MAXVAL(evals)) THEN
       			       phi_fse(particle) = ATAN2(evects(3,1),evects(1,1))
       			       fse_axis(particle,:) = evects(:,1)
		       ELSE IF(evals(2) .EQ. MAXVAL(evals)) THEN
			       phi_fse(particle) = ATAN2(evects(3,2),evects(1,2))
			       fse_axis(particle,:) = evects(:,2)
		       ELSE
			       phi_fse(particle) = ATAN2(evects(3,3),evects(1,3))
			       fse_axis(particle,:) = evects(:,3)
		       END IF


		       !!! natural strain = ln(a/c) where a is the long axis = maxval(evals)**0.5
         	       ln_fse(particle) = 0.5*LOG(MAXVAL(evals)/MINVAL(evals))


		       !!! Cijkl tensor (using Voigt average)
			CALL CPU_TIME ( time_begin )
        		       CALL voigt(particle)
	        	CALL CPU_TIME ( time_end)
			total_time_voigt = total_time_voigt + (time_end-time_begin)

 			!!! Percentage of anisotropy and orientation of axis of hexagonal symmetry
 	        	CALL CPU_TIME ( time_begin)
    		    		CALL DECSYM(Sav_particles(particle,:,:),temp_perc, temp_phi, temp_Aaxis, frac_symm)
	        	CALL CPU_TIME ( time_end)
			total_time_decsym = total_time_decsym + (time_end-time_begin)

         	 
			!!! If we encounter a problem with DECSYM still projecting onto the slow axis 
			!!! instead of the fast axis, we should replace the above call with the following
			!!  call to Thorsten Becker's modified version
		 
!		        CALL drex_decsym_ftrn(Sav_particles(particle,:,:),k,g,frac_symm,tiaxis,&
!         	        		     dummyt,dummyt,dummyt,dummyt,dummyt,dummyt,vel,dummyt,scc_rot)
!         	        temp_perc = frac_symm(2) * 100d0 ! VTI percent inclination angle with third axis
!   vel(1) = sqrt(xhex(3))                  ! vp1*sqrt(dens)
!   vel(2) = sqrt(xhex(1))                  ! vp2*sqrt(dens)

!         	        if (vel(2) .LT. vel(1)) then
!				temp_phi = asin(tiaxis(3))
!			else
!				temp_phi = asin(tiaxis(6))
!			end if
		       !!! ----------end of alternative code section----------------
		
		       perc_a(particle) = temp_perc
		       phi_a(particle)  = temp_phi
		       hex_axis(particle,:) = temp_Aaxis

		       WRITE (17, "(3(1PE12.4))"), X1p(particle), X2p(particle), X3p(particle)
		       WRITE (18, "(4(1PE12.4))"), fse_axis(particle,:), ln_fse(particle)
		       WRITE (19, "(8(1PE12.4))"), phi_a(particle), perc_a(particle) , frac_symm(1:6)
		     
		END DO !!! end of loop over timestep iterations
  	END DO	!!! end of loop over particles
   END DO		!! end of loop over timesteps

   print *, 'Total times: '
   print *, 'pathline: ', total_time_pathline
   print *, 'strain: ', total_time_strain
   print *, 'Voigt: ', total_time_voigt
   print *, 'DECSYM:', total_time_decsym

!   CALL RecalcEulerAnglesFromACS
  OPEN (23, file='ODF.dat', status='replace')

   DO i=1,size
   	WRITE (21, "(9(1PE12.4))"), acs(i,1:3,1:3)
	WRITE (23, "(1(1PE12.4))"), odf(i)
   END DO
   
   DEALLOCATE(xe1,xe2,xe3)

 
   STOP

   CLOSE (17)
   CLOSE (18)
   CLOSE (19)
   CLOSE (21)
   CLOSE (22)
   CLOSE (23)
   
  END PROGRAM DR_CF
  
  
SUBROUTINE RecalcEulerAnglesFromACS
   USE comvar
   USE GENERAL_VARS

   IMPLICIT NONE
   INTEGER :: i ! loop counters

 	DO i=1,size
		xe2(i) = ACOS(acs(i,3,3));			!theta
		xe1(i) = ACOS(acs(i,3,2)/(-SIN(xe2(i))));	!phi1
	 	xe3(i) = ACOS(acs(i,2,3)/(SIN(xe2(i))));	!phi2
	END DO

   RETURN

   END SUBROUTINE RecalcEulerAnglesFromACS
