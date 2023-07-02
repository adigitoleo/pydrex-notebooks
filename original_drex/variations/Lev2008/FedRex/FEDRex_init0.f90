!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine INIT0, initialization of variables, arrays and parameters   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE init0

   USE comvar
   USE GENERAL_VARS

   IMPLICIT NONE

   INTEGER :: i1,i3,i,j1,j2,j3, particle ! loop counters


   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ran0
   ! matrix of random numbers used to generate initial random LPO

   DOUBLE PRECISION :: phi1,theta,phi2
   ! eulerian angles

 
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: temp

!!! Name of input files from file input.dat

   OPEN(15,file="input.dat")

   READ(15,9000) fname
   READ(15,*) nx1,nx2,nx3
   READ(15,*) n_particles
   READ(15,*) timestep_first, timestep_last, timestep_jump
   READ(15,*) size3
   READ(15,*) tau(1),tau(2),tau(3),tau(4),tau_ens
   READ(15,*) Mob
   READ(15,*) chi
   READ(15,*) dimensions
   READ(15,*) random_flag


!!! Sanity checks
 IF ((dimensions .EQ. 2) .AND. (nx2 .NE. 1)) THEN
 	print *, 'Asked to run in two dimensions, but number of nodes in Y (into screen) direction not equal 1.' 
	print *, 'FedRex is now forcing nx2 to be 1, but please fix the input file'
	nx2 = 1;
END IF

!!! Nucleation parameter set to 5 as in Kaminski and Ribe, 2001
   lambda = 5d0

!!! initial size
   size = size3**3

!!! stress exponent for olivine and enstatite
   stressexp = 3.5d0

!!! initialization of streamline data
   stream_Lij=0d0 ; stream_dt=0d0 ; stream_e=0d0
   IF (dimensions .EQ. 2) THEN
	nx2 = 1;
   END IF
	
	
   ALLOCATE(X1(nx1),X2(nx2),X3(nx3),Ui(3,nx1,nx2,nx3),Dij(3,3,nx1,nx2,nx3))
   ALLOCATE(X1gr(nx1),X2gr(nx2),X3gr(nx3))
   ALLOCATE(GOL(nx1,nx2,nx3))
   
   ALLOCATE(temp(n_particles,4), X1p(n_particles), X2p(n_particles),X3p(n_particles), Xol(n_particles),box_drop_status(n_particles))
   ALLOCATE(Fij(3,3,n_particles))
   ALLOCATE(phi_fse(n_particles),ln_fse(n_particles),phi_a(n_particles),perc_a(n_particles))
   ALLOCATE(fse_axis(n_particles,3))
   ALLOCATE(hex_axis(n_particles,3))
   
!!! Initialization of grid, velocity field and velocity gradient
   X1 = 0d0 ; X3 = 0d0 ; Ui = 0d0 ; Dij = 0d0
   CALL readveloc (fname,  timestep_first, dimensions)

!! Read initial particle locations
   OPEN (9,file="particles.dat")
   DO particle = 1, n_particles
	read(9,8000) temp(particle,1), temp(particle,2), temp(particle,3), temp(particle,4)
	X1p(particle) = temp(particle,1)
	X2p(particle) = temp(particle,2)
	X3p(particle) = temp(particle,3)
	Xol(particle) = temp(particle,4)/1d2
   END DO
   CLOSE(9)
   box_drop_status = 0


!!! Initial deformation gradient tensor
   Fij = 0d0 ; Fij(1,1,:) = 1d0 ; Fij(2,2,:) = 1d0 ; Fij(3,3,:) = 1d0

!!! Initial Grain lag PI parameter

   GOL = 1d1 ;

!!! tensor \epsilon_{ijk}

   alt=0d0
   alt(1,2,3) = 1d0 ; alt(2,3,1) = 1d0 ; alt(3,1,2) = 1d0
   alt(1,3,2) = -1d0 ; alt(2,1,3) = -1d0 ; alt(3,2,1) = -1d0

!!! tensor \delta_{ij}

   del=0d0
   del(1,1) = 1d0 ; del(2,2) = 1d0 ; del(3,3) = 1d0

!!! tensors of indices

   ijkl(1,1) = 1 ; ijkl(1,2) = 6 ; ijkl(1,3) = 5
   ijkl(2,1) = 6 ; ijkl(2,2) = 2 ; ijkl(2,3) = 4
   ijkl(3,1) = 5 ; ijkl(3,2) = 4 ; ijkl(3,3) = 3

   l1(1) = 1 ; l1(2) = 2 ; l1(3) = 3
   l1(4) = 2 ; l1(5) = 3 ; l1(6) = 1
   l2(1) = 1 ; l2(2) = 2 ; l2(3) = 3
   l2(4) = 3 ; l2(5) = 1 ; l2(6) = 2

!!! Stiffness matrix for Olivine (GigaPascals)

   S0(1,1) = 320.71d0 ; S0(1,2) = 69.84d0  ; S0(1,3) = 71.22d0
   S0(2,1) = S0(1,2)  ; S0(2,2) = 197.25d0 ; S0(2,3) = 74.8d0
   S0(3,1) = S0(1,3)  ; S0(3,2) = S0(2,3)  ; S0(3,3) = 234.32d0
   S0(4,4) = 63.77d0  ; S0(5,5) = 77.67d0  ; S0(6,6) = 78.36d0

!!! Stiffness matrix for Enstatite (GPa)

   S0_ens(1,1) = 236.9d0 ; S0_ens(1,2) = 79.6d0  ; S0_ens(1,3) = 63.2d0
   S0_ens(2,1) = S0_ens(1,2)  ; S0_ens(2,2) = 180.5d0 ; S0_ens(2,3) = 56.8d0
   S0_ens(3,1) = S0_ens(1,3)  ; S0_ens(3,2) = S0_ens(2,3)  ; S0_ens(3,3) = 230.4d0
   S0_ens(4,4) = 84.3d0  ; S0_ens(5,5) = 79.4d0  ; S0_ens(6,6) = 80.1d0

!!! allocation of the dimensions of the arrays

   ALLOCATE(xe1(size),xe2(size),xe3(size))

   ALLOCATE(odf(size),odfi(size),dotodf(size))
   ALLOCATE(odf_ens(size),odfi_ens(size),dotodf_ens(size))

   ALLOCATE(rt(size),rt_ens(size))

   ALLOCATE(ran0(3*size))

   ALLOCATE(acs(size,3,3),dotacs(size,3,3),acsi(size,3,3),acs0(size,3,3))
   ALLOCATE(acs_ens(size,3,3),dotacs_ens(size,3,3),acsi_ens(size,3,3))
   ALLOCATE(acs_step(1000,size,3,3))
   ALLOCATE(odf_steps(1000,size))

   ALLOCATE(acs_particles(n_particles,size,3,3))
   ALLOCATE(odf_particles(n_particles,size))

   ALLOCATE(acs_ens_particles(n_particles,size,3,3))
   ALLOCATE(odf_ens_particles(n_particles,size))

  ALLOCATE(Sav_particles(n_particles,6,6))
 
!!! initialization of orientations - uniformally random distribution
!!! Rmq cos(theta) used to sample the metric Eulerian space

   CALL RANDOM_NUMBER(ran0)

   i = 1
   IF (random_flag .EQ. 0) THEN
       print *,'Random orientation used'
       DO j1 =1, size3 ; DO j2 =1, size3 ; DO j3 =1, size3          
            xe1(i) = (REAL(j1)-ran0(i))/REAL(size3)*ACOS(-1d0)
            xe2(i) = ACOS(-1d0 + (REAL(j2)-ran0(size+i))/REAL(size3)*2d0)
            xe3(i) = (REAL(j3)-ran0(i+2*size))/REAL(size3)*ACOS(-1d0)
            i = i + 1
        END DO ; END DO ; END DO
    ELSE
	IF (random_flag .EQ. 1) THEN
		print *,'45 degrees uniform orientation used'
		DO j1 =1, size3 ; DO j2 =1, size3 ; DO j3 =1, size3
        	  xe1(i) = 0d0
        	  xe2(i) = ACOS(-1d0)/4
        	  xe3(i) = ACOS(-1d0)/2
        	  i = i + 1
        	END DO ; END DO ; END DO
	ELSE
		READ(15, 9000), EulerAnglesFileName
		print *, 'Using Euler angles from file: ', EulerAnglesFileNAME
		OPEN (7, file=EulerAnglesFileName)
		OPEN (23, file='ODF.dat')
		acs0=0d0;
		DO particle=1,n_particles
			rt = 0d0 ; 
			rt_ens = rt ; 
			DO i=1,size
				READ (7, *), acs0(i,:,:)
				acs_particles(particle,i,:,:) = TRANSPOSE(acs0(i,:,:));
			END DO	

		   	READ (23, *), odf_particles(particle,:)
			odf_ens_particles(particle,:) = 1d0/REAL(size3**3)
			acs_ens_particles(particle,:,:,:) = acs0

		 END DO
	
		CLOSE (7)
	END IF
    END IF





IF (random_flag .NE. 2) THEN
   !!! Random initial LPO
   rt = 0d0 ; 
   rt_ens = rt ; 

   DO i = 1 , size
      phi1 = xe1(i) ; theta = xe2(i) ; phi2 = xe3(i)
!!! Direction cosine matrix
      acs0(i,1,1)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
      acs0(i,1,2)=COS(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*SIN(phi2)
      acs0(i,1,3)=SIN(phi2)*SIN(theta)

      acs0(i,2,1)=-SIN(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*COS(phi2)
      acs0(i,2,2)=-SIN(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*COS(phi2)
      acs0(i,2,3)=COS(phi2)*SIN(theta)

      acs0(i,3,1)=SIN(theta)*SIN(phi1)
      acs0(i,3,2)=-SIN(theta)*COS(phi1)
      acs0(i,3,3)=COS(theta)

   END DO
   DO particle=1,n_particles
	odf_particles(particle,:) = 1d0/REAL(size3**3)
	acs_particles(particle,:,:,:) = acs0
	odf_ens_particles(particle,:) = 1d0/REAL(size3**3)
	acs_ens_particles(particle,:,:,:) = acs0
   END DO
END IF


   DEALLOCATE(ran0)

  CLOSE(15) !that's input.dat
  close(23) !that's ODF.dat

   RETURN

 7000 FORMAT(16(1pe14.6))
 8000 FORMAT(13(1pe14.6))
 9000 FORMAT(A80)
 END SUBROUTINE init0
   
