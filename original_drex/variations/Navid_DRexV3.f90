
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DRex - Version used in http://dx.doi.org/10.1016/j.epsl.2016.12.004    !!!
!!!                                                                        !!!
!!! Implements some additional finite strain and stiffness calculations    !!!
!!!                                                                        !!!
!!! Calculates the evolutions of texture an olivine                        !!!
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

!!! Grid points - Velocity field - Velocity gradient tensor

   INTEGER :: nx1, nx3 !! grid size in X1 and X3 direction
   INTEGER :: stepx1, stepx3 
   !! interval of calculation of the LPO on the grid
   !! LPO calculates only every stepx1 and stepx3 point
   INTEGER :: i1first,i3first,i1last,i3last
   !! indices of the first and last points at which LPO is calculated

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1
   ! grid points coordinates - distance from the ridge
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X3
   ! grid points coordinates - depth

   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Ui
   ! velocity vector - U1(x1,x3)=Ui(1,x1,x3)

   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Dij
   ! velocity gradient tensor - D11(x1,x3)=Dij(1,1,x1,x3)

   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Fij
   ! finite strain tensor - F11(x1,x3)=Fij(1,1,x1,x3)

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Temperature
   ! grid points of temperature values

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DeformationMechanism
   ! grid points of deformation regime (100 = dislocation, 0 = diffusion)

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: radani
   ! grid points of radial anisotropy parameter

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: percani
   ! grid points of percentage azimuthal anisotropy parameter

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: grout
   ! dummy used to activate graphic output only for calculated points

!!! LPO calculation

   DOUBLE PRECISION, DIMENSION(3,3) :: l,e
   ! velocity gradient tensor and strain rate tensor

   DOUBLE PRECISION :: epsnot
   ! reference strain rate

   DOUBLE PRECISION, DIMENSION(3,3,10000) :: stream_Lij
   DOUBLE PRECISION, DIMENSION(10000) :: stream_dt
   DOUBLE PRECISION, DIMENSION(10000) :: stream_e
   ! values of Lij, time spent and strain rate at each location on the streamline

!!! Dynamic recrystallization
   
   INTEGER :: size3, size ! size = size3^3
   ! number of points in the (metric) Eulerian space

   DOUBLE PRECISION :: lambda, Mob, chi
   ! lambda = nucleation parameter
   ! Mob = grain mobility
   ! threshold volume fraction for activation of grain boundary sliding

   DOUBLE PRECISION :: Xol
   ! fraction of olivine in the aggregate

   DOUBLE PRECISION, DIMENSION(3,3,3) :: alt ! \epsilon_{ijk} tensor
   DOUBLE PRECISION, DIMENSION(3,3) :: del ! \delta_{ij} tensor
   INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
   INTEGER, DIMENSION(6) :: l1,l2 ! tensot of indices to form Sij from Cijkl

   DOUBLE PRECISION, DIMENSION(4) :: tau
   ! RSS for the 4 slip systems of olivine (only 3 are activated)
   DOUBLE PRECISION :: tau_ens ! RSS of enstatite slip system
   DOUBLE PRECISION :: stressexp ! stress exponent for olivine and enstatite

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: odf,dotodf,odfi
   ! volume fraction of the olivine grains
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: odf_ens,dotodf_ens,odfi_ens
   ! volume fraction of the enstatite grains

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rt,rt_ens
   ! dislocations density for olivine and enstatite

   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: acs,dotacs,acsi,acs0
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: acs_ens,dotacs_ens,acsi_ens
   !! matrix of direction cosine

!!! Output

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_fse
   ! orientation of the long axis of the FSE on the (X1,X3) plan

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ln_fse
   ! finite strain = ln(longaxis/shortaxis)

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_a
   ! average orientation of a-axis on (X1,X3) plan

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: perc_a
   ! percentage of S wave anisotropy

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GOL
   ! Grain Orientation Lag i.e. Pi parameter

!!! Elastic tensor 

   DOUBLE PRECISION, DIMENSION(6,6) :: S0,Sav,S0_ens
   ! stiffness matrix
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: Cav
   !Cijkl tensor at point of calculation

   END MODULE comvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Main program                                                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   PROGRAM DR_CF 

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: X10,X30 
   !! coordinates of the grid point where calculation is started

   INTEGER :: i,i1,i3 !! loop counters
   INTEGER :: i1a,i3a !! indices of the grid point at the end of streamline

   INTEGER :: alpha
   ! 1 or 0 depending on the type of deformation regime

   LOGICAL :: setDiffusionCreep
   ! true if under diffusion creep regime

   INTEGER :: step1 !! number of points on a streamline

   DOUBLE PRECISION, DIMENSION(3,3) :: LSij
   ! left-strech tensor for FSE calculation

   INTEGER :: nrot
   ! number of rotations for the Jacobi

   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects
   ! eigen values and vectors in jacobi

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1gr,X3gr
   !! actual coordinates of calculation point when not a grid point

!!! initialization

   CALL init0

   ALLOCATE(X1gr(nx1),X3gr(nx3))

   OPEN(18,FILE='Cijkl',form='unformatted')

   X1gr = 0d0 ; X3gr = 0d0
   phi_fse = 0d0 ; ln_fse = 0d0 ; phi_a = 0d0 ; perc_a = 0d0
   grout = 0d0 ; GOL = 0d0

!!! Calculation of GOL parameter

   DO i1 = 1, nx1-1
      DO i3 = 1 , nx3-1

!!! Location of the calculation point 
         X10 = (X1(i1)+X1(i1+1))/2d0 ; X30 = (X3(i3)+X3(i3+1))/2d0
         X1gr(i1) = X10 ; X3gr(i3) = X30

         CALL pipar(i1,i3,X1gr(i1),X3gr(i3))

      END DO
   END DO

!!! FSE and LPO calculation

   DO i1 = i1first , i1last-1 , stepx1
      DO i3 = i3first , i3last-1 , stepx3

!!! Initial location of the grid point 
         X10 = (X1(i1)+X1(i1+1))/2d0 ; X30 = (X3(i3)+X3(i3+1))/2d0
         X1gr(i1) = X10 ; X3gr(i3) = X30

!!! Backward calculation of the pathline for each tracer

         step1 = 0

         CALL pathline(X10,X30,i1,i3,step1,i1a,i3a)

!!! Inward calculation of the LPO

!!! Random initial LPO

         rt = 0d0 ; odf = 1d0/REAL(size3**3) ; acs = acs0
         rt_ens = rt ; odf_ens = odf ; acs_ens = acs0

         if (DeformationMechanism(i1,i3)>0.d0)  then
            alpha = 1
            setDiffusionCreep = .false.
         else
            alpha = 0
            setDiffusionCreep = .true.
         end if

         CALL strain(step1,Fij(:,:,i1,i3), alpha, setDiffusionCreep)

!!! Left-stretch tensor for FSE calculation

         LSij = MATMUL(Fij(:,:,i1,i3),TRANSPOSE(Fij(:,:,i1,i3)))

         CALL JACOBI(LSij,3,3,evals,evects,nrot)

!!! Pick up the orientation of the long axis of the FSE

         IF (evals(1) .EQ. MAXVAL(evals)) THEN
            phi_fse(i1,i3) = ATAN2(evects(3,1),evects(1,1))
            ELSE IF(evals(2) .EQ. MAXVAL(evals)) THEN
               phi_fse(i1,i3) = ATAN2(evects(3,2),evects(1,2))
               ELSE
               phi_fse(i1,i3) = ATAN2(evects(3,3),evects(1,3))
         END IF

!!! natural strain = ln(a/c) where a is the long axis = maxval(evals)**0.5

         ln_fse(i1,i3) = 0.5*LOG(MAXVAL(evals)/MINVAL(evals))

!!! Activate graphic output
         grout(i1,i3) = 1d0

!!! Cijkl tensor (using Voigt average)
         CALL voigt
         WRITE(18) X1gr(i1),X3gr(i3),Sav

!!! Percentage of anisotropy and orientation of axis of hexagonal symmetry
         CALL DECSYM(Sav,perc_a(i1,i3),phi_a(i1,i3),radani(i1,i3),percani(i1,i3))

      END DO
   END DO

   CLOSE(18)

!!! Graphic output
   
   CALL check_veloc

   OPEN(10,file='fse.m',status='replace')
   OPEN(11,file='a_axis.m',status='replace')
   OPEN(12,file='GOL.m',status='replace')
   OPEN(13,file='ani.m',status='replace')
   OPEN(14,file='sot.m',status='replace')
   

   WRITE(10,*)
   WRITE(10,*) "clear"
   WRITE(10,*) "figure"

   WRITE(11,*)
   WRITE(11,*) "clear"
   WRITE(11,*) "figure"

   WRITE(12,*)
   WRITE(12,*) "clear"
   WRITE(12,*) "figure"

   WRITE(13,*)
   WRITE(13,*) "clear"
   WRITE(13,*) "figure"

   WRITE(14,*)
   WRITE(14,*) "clear"
   WRITE(14,*) "figure"

   WRITE(10,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(10,"(121(1PE12.4))") grout(:,i3)*ln_fse(:,i3)*COS(phi_fse(:,i3))
   END DO
   WRITE(10,*) "];"

   WRITE(10,*) "Uz=["
   DO i3 = 1,nx3
      WRITE(10,"(121(1PE12.4))") grout(:,i3)*ln_fse(:,i3)*SIN(phi_fse(:,i3))
   END DO
   WRITE(10,*) "];"

   WRITE(11,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(11,"(121(1PE12.4))") grout(:,i3)*COS(phi_a(:,i3))*perc_a(:,i3)
   END DO
   WRITE(11,*) "];"

   WRITE(11,*) "Uz=["
   DO i3 = 1,nx3
      WRITE(11,"(121(1PE12.4))") grout(:,i3)*SIN(phi_a(:,i3))*perc_a(:,i3)
   END DO
   WRITE(11,*) "];"

   WRITE(12,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(12,"(121(1PE12.4))") GOL(:,i3)
   END DO
   WRITE(12,*) "];"

   WRITE(13,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(13,"(121(1PE12.4))") percani(:,i3)
   END DO
   WRITE(13,*) "];"

   WRITE(14,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(14,"(121(1PE12.4))") radani(:,i3)
   END DO
   WRITE(14,*) "];"

   WRITE(10,*) "X=["
   WRITE(10,"(121(1PE12.4))") X1gr
   WRITE(10,*) "];"

   WRITE(10,*) "Z=["
   WRITE(10,"(121(1PE12.4))") X3gr
   WRITE(10,*) "];"

   WRITE(11,*) "X=["
   WRITE(11,"(121(1PE12.4))") X1gr
   WRITE(11,*) "];"

   WRITE(11,*) "Z=["
   WRITE(11,"(121(1PE12.4))") X3gr
   WRITE(11,*) "];"

   WRITE(12,*) "X=["
   WRITE(12,"(121(1PE12.4))") X1gr
   WRITE(12,*) "];"

   WRITE(12,*) "Z=["
   WRITE(12,"(121(1PE12.4))") X3gr
   WRITE(12,*) "];"

   WRITE(13,*) "X=["
   WRITE(13,"(121(1PE12.4))") X1gr
   WRITE(13,*) "];"

   WRITE(13,*) "Z=["
   WRITE(13,"(121(1PE12.4))") X3gr
   WRITE(13,*) "];"

   WRITE(14,*) "X=["
   WRITE(14,"(121(1PE12.4))") X1gr
   WRITE(14,*) "];"

   WRITE(14,*) "Z=["
   WRITE(14,"(121(1PE12.4))") X3gr
   WRITE(14,*) "];"

   WRITE(10,*) "quiver(X,Z,Ux,Uz,'.')"
   WRITE(10,*) "axis([0 1200000 0 400000])"
   WRITE(10,*) "xlabel('X1, dimensionless distance from the ridge axis')"
   WRITE(10,*) "ylabel('X3, dimensionless depth')"
   WRITE(10,*) "title('Orientation of the long axis of the FSE')"
   
   WRITE(10,*) "axis image"
   WRITE(10,*) "axis ij"
   
   WRITE(11,*) "quiver(X,Z,Ux,Uz,'.')"
   WRITE(11,*) "axis([0 1200000 0 400000])"
   WRITE(11,*) "xlabel('X1, dimensionless distance from the ridge axis')"
   WRITE(11,*) "ylabel('X3, dimensionless depth')"
   WRITE(11,*) "title('Orientation of the hexagonal symmetry axis')"
   
   WRITE(11,*) "axis image"
   WRITE(11,*) "axis ij"

   WRITE(12,*) "surf(X,Z,log(abs(Ux)))"
   WRITE(12,*) "axis([0 1200000 0 400000])"
   WRITE(12,*) "caxis([-6 3])"
   WRITE(12,*) "shading 'interp'"
   WRITE(12,*) "colorbar 'horiz'"
   WRITE(12,*) "xlabel('X1, dimensionless distance from the ridge axis')"
   WRITE(12,*) "ylabel('X3, dimensionless depth')"
   WRITE(12,*) "title('Values of GOL parameter (in log units)')"
   
   WRITE(12,*) "axis image"
   WRITE(12,*) "axis ij"

   WRITE(13,*) "surf(X,Z,Ux)"
   WRITE(13,*) "axis([0 1200000 0 400000])"
   WRITE(13,*) "shading 'interp'"
   WRITE(13,*) "colorbar 'horiz'"
   WRITE(13,*) "xlabel('Distance from the ridge axis (m)')"
   WRITE(13,*) "ylabel('Depth (m)')"
   WRITE(13,*) "title('Percentage azimuthal anisotropy')"
   
   WRITE(13,*) "axis image"
   WRITE(13,*) "axis ij"

   WRITE(14,*) "surf(X,Z,Ux)"
   WRITE(14,*) "axis([0 1200000 0 400000])"
   WRITE(14,*) "shading 'interp'"
   WRITE(14,*) "colorbar 'horiz'"
   WRITE(14,*) "xlabel('Distance from the ridge axis (m)')"
   WRITE(14,*) "ylabel('Depth (m)')"
   WRITE(14,*) "title('Radial anisotropy')"
   
   WRITE(14,*) "axis image"
   WRITE(14,*) "axis ij"

   CLOSE(10)
   CLOSE(11)
   CLOSE(12)
   CLOSE(13)
   CLOSE(14)

   STOP

   END PROGRAM DR_CF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine INIT0, initialization of variables, arrays and parameters   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE init0

   USE comvar

   IMPLICIT NONE

   INTEGER :: i1,i3,i,j1,j2,j3 ! loop counters

   CHARACTER(30) :: fname
   ! name of the imput file containing X1, X3, Ui and Dij

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ran0
   ! matrix of random numbers used to generate initial random LPO

   DOUBLE PRECISION :: phi1,theta,phi2
   ! eulerian angles

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xe1,xe2,xe3
   ! matrixes of initial random eulerian angles

!!! Name of input files from file input.dat

   OPEN(15,file="input2.dat")

   READ(15,*) fname
   READ(15,*) nx1,nx3,stepx1,stepx3
   READ(15,*) i1first,i3first,i1last,i3last
   READ(15,*) size3,Xol
   READ(15,*) tau(1),tau(2),tau(3),tau(4),tau_ens
   READ(15,*) Mob
   READ(15,*) chi

   CLOSE(15)

!!! Nucleation parameter set to 5 as in Kaminski and Ribe, 2001
   lambda = 5d0

!!! initial size
   size = size3**3

!!! Xol given as a percentage
   Xol = Xol/1d2

!!! stress exponent for olivine and enstatite
   stressexp = 3.5d0

!!! initialization of streamline data
   stream_Lij=0d0 ; stream_dt=0d0 ; stream_e=0d0

   ALLOCATE(X1(nx1),X3(nx3),Ui(3,nx1,nx3),Dij(3,3,nx1,nx3))
   ALLOCATE(Fij(3,3,nx1,nx3),grout(nx1,nx3),GOL(nx1,nx3))
   ALLOCATE(phi_fse(nx1,nx3),ln_fse(nx1,nx3),phi_a(nx1,nx3),perc_a(nx1,nx3))
   ALLOCATE(Temperature(nx1,nx3))
   ALLOCATE(DeformationMechanism(nx1,nx3))
   ALLOCATE(radani(nx1,nx3))
   ALLOCATE(percani(nx1,nx3))

!!! Initialization of grid, velocity field and velocity gradient

   X1 = 0d0 ; X3 = 0d0 ; Ui = 0d0 ; Dij = 0d0
   radani = 0d0 ; percani = 0d0

   OPEN(8,file=fname)

   do i1 = 1, nx1
      do i3 = 1, nx3
         read(8,*) X1(i1),X3(i3),Ui(1,i1,i3),Ui(3,i1,i3), &
              Dij(1,1,i1,i3),Dij(1,3,i1,i3),Dij(3,1,i1,i3),Temperature(i1,i3), &
              DeformationMechanism(i1,i3)
         Dij(3,3,i1,i3) = -Dij(1,1,i1,i3)-Dij(2,2,i1,i3)
      END DO
   END DO

   CLOSE(8)

   !make variable alpha based on the deformation mechanism.

!!! Initial deformation gradient tensor
   Fij = 0d0 ; Fij(1,1,:,:) = 1d0 ; Fij(2,2,:,:) = 1d0 ; Fij(3,3,:,:) = 1d0

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

!!! initialization of orientations - uniformally random distribution
!!! Rmq cos(theta) used to sample the metric Eulerian space

   CALL RANDOM_NUMBER(ran0)

   i = 1

   DO j1 =1, size3 ; DO j2 =1, size3 ; DO j3 =1, size3
      xe1(i) = (REAL(j1)-ran0(i))/REAL(size3)*ACOS(-1d0)
      xe2(i) = ACOS(-1d0 + (REAL(j2)-ran0(size+i))/REAL(size3)*2d0)
      xe3(i) = (REAL(j3)-ran0(i+2*size))/REAL(size3)*ACOS(-1d0)
      i = i + 1
   END DO ; END DO ; END DO

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
   
   DEALLOCATE(xe1,xe2,xe3,ran0)

   RETURN

 1000 FORMAT(16(1pe14.6))

   END SUBROUTINE init0
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine PATHLINE - Calculation of tracers pathlines                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE pathline(X1s,X3s,i1s,i3s,step2,i1r,i3r)

   USE comvar 

   IMPLICIT NONE

   INTEGER :: step2
   ! number of steps to construct the streamline

   INTEGER :: i1s,i3s,i1,i3,i1r,i3r
   ! initial, intermediate and final indices of the closest upper-left grid point

   DOUBLE PRECISION :: ds,dt
   ! elementary displacement along the streamline and corresponding time step

   DOUBLE PRECISION :: X1s,X3s,X1i,X3i
   ! initial and intermediate position of the tracer

   DOUBLE PRECISION :: U1i,U3i
   ! intermediate velocity components on the streamline

   DOUBLE PRECISION :: kx1,kx2,kx3,kx4
   DOUBLE PRECISION :: kz1,kz2,kz3,kz4
   ! RGK intermediate for streamlines calculation 

   DOUBLE PRECISION :: max_strain
   ! value of the dimensionless strain for which any former
   ! LPO has been erased. From Kaminski and Ribe 2002, 10
   ! is a safe value. Should be tested for a given flow though.

!! initial indices and position of the tracers

   i1 = i1s ; i3 = i3s ; max_strain = 0d0 

!!! construction of the streamline

   DO WHILE (X3s .LT. X3(nx3) .AND. X1s .GT. X1(1) .AND. max_strain .LT. 1d1 .AND. step2 .LT. 10000)

      X1i = X1s ; X3i = X3s
      step2 = step2 + 1

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
      CALL gradientcalc(X1i,X3i,i1,i3)

      ds = 0.25d0*MIN(ABS(X1(i1+1)-X1(i1)),ABS(X3(i3+1)-X3(i3)))
      dt = ds/SQRT(U1i**2+U3i**2)

!! record of the local velocity gradient tensor and time spent at that point
      stream_Lij(:,:,step2) = l(:,:) ; stream_dt(step2) = dt 
      stream_e(step2) = epsnot

      max_strain = max_strain + dt*epsnot

      kx1 = -U1i*dt ; kz1 = -U3i*dt

      X1i = X1s + 0.5d0*kx1
      X3i = X3s + 0.5d0*kz1

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)

      kx2 = -U1i*dt ; kz2 = -U3i*dt

      X1i = X1s + 0.5d0*kx2
      X3i = X3s + 0.5d0*kz2

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)

      kx3 = -U1i*dt ; kz3 = -U3i*dt

      X1i = X1s + kx3
      X3i = X3s + kz3

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)

      kx4 = -U1i*dt ; kz4 = -U3i*dt

      X1s = X1s + (kx1/2d0+kx2+kx3+kx4/2d0)/3d0
      X3s = X3s + (kz1/2d0+kz2+kz3+kz4/2d0)/3d0

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1s .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1s .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3s .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3s .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      i1r = i1 ; i3r = i3

   END DO

   RETURN

   END SUBROUTINE pathline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine VELOCITYCALC, calculation of velocity at a given point      !!!
!!! by interpolation method given in Numerical Recipies, Press et al., p96 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE velocitycalc(X1s,X3s,U1s,U3s,i1s,i3s)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X3s
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: U1s,U3s
   ! interpolated velocity at the point

   INTEGER :: i1s,i3s
   ! indices of the UP-LEFT grid point closest to the extrapolation point

   DOUBLE PRECISION :: y1,y2,y3,y4
   ! dummies for interpolation

   y1 = Ui(1,i1s,i3s) ; y2 = Ui(1,i1s+1,i3s)
   y3 = Ui(1,i1s+1,i3s+1) ; y4 = Ui(1,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,U1s)

   y1 = Ui(3,i1s,i3s) ; y2 = Ui(3,i1s+1,i3s)
   y3 = Ui(3,i1s+1,i3s+1) ; y4 = Ui(3,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,U3s)

   RETURN

   END SUBROUTINE velocitycalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine interp, bilinear interpolation of a function                !!!
!!! following the method given in Numerical Recipies, Press et al., p96    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,res)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X3s
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: res
   ! result of the interpolation at the point

   INTEGER :: i1s,i3s
   ! indices of the UP-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4,tt,uu
   ! dummies for interpolation (numerical recipies, p96)

   tt = (X1s-X1(i1s))/(X1(i1s+1)-X1(i1s))
   uu = (X3s-X3(i3s))/(X3(i3s+1)-X3(i3s))

   res = (1d0-tt)*(1d0-uu)*y1+tt*(1d0-uu)*y2+tt*uu*y3+(1d0-tt)*uu*y4

   RETURN

   END SUBROUTINE interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine GRADIENTCALC, interpolation of velocity gradient tensor     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE gradientcalc(X1s,X3s,i1s,i3s)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X3s
   ! x and y coordinates of the point on the streamline

   INTEGER :: i1s,i3s
   ! indices of the UPPER-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4
   ! dummies for interpolation

   DOUBLE PRECISION, DIMENSION(3,3) :: ee
   ! dummy for reference strain rate calculation

   INTEGER :: nrot
   ! number of rotations for the Jacobi
   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects
   ! eigen values and vectors in jacobi

   l = 0d0 ; e = 0d0

   y1 = Dij(1,1,i1s,i3s) ; y2 = Dij(1,1,i1s+1,i3s)
   y3 = Dij(1,1,i1s+1,i3s+1) ; y4 = Dij(1,1,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,l(1,1))

   y1 = Dij(1,3,i1s,i3s) ; y2 = Dij(1,3,i1s+1,i3s)
   y3 = Dij(1,3,i1s+1,i3s+1) ; y4 = Dij(1,3,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,l(1,3))

   y1 = Dij(3,1,i1s,i3s) ; y2 = Dij(3,1,i1s+1,i3s)
   y3 = Dij(3,1,i1s+1,i3s+1) ; y4 = Dij(3,1,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,l(3,1))

   y1 = Dij(2,3,i1s,i3s) ; y2 = Dij(2,3,i1s+1,i3s)
   y3 = Dij(2,3,i1s+1,i3s+1) ; y4 = Dij(2,3,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,l(2,3))

   l(3,3) = -l(1,1)-l(2,2)

!!! strain rate tensor
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
!!! subroutine STRAIN - Calculation of strain along pathlines              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE strain(step2,fse, alpha, setDiffusionCreep)

   USE comvar 

   IMPLICIT NONE

   INTEGER :: j,j1,j2,nn,N_strain 
   ! loop counters

   INTEGER :: step2
   ! number of steps used to construct the streamline

   INTEGER :: alpha
   ! 1 or 0 depending on the type of deformation regime

   LOGICAL :: setDiffusionCreep
   ! true if in diffusion creep regime, false if otherwise
   
   DOUBLE PRECISION :: t_loc,dt
   ! time spent on a given point of the streamline and DRex time step

   DOUBLE PRECISION, DIMENSION(3,3) :: fse,fsei
   ! local finite deformation tensor

   DOUBLE PRECISION, DIMENSION(3,3) :: kfse1,kfse2,kfse3,kfse4
   ! RGK intermediates for Finite Strain Ellipsoid

   DOUBLE PRECISION, DIMENSION(size) :: kodf1,kodf2,kodf3,kodf4
   DOUBLE PRECISION, DIMENSION(size) :: kodf1_ens,kodf2_ens,kodf3_ens,kodf4_ens
   ! intermediates for odf

   DOUBLE PRECISION, DIMENSION(size,3,3) :: kac1,kac2,kac3,kac4
   DOUBLE PRECISION, DIMENSION(size,3,3) :: kac1_ens,kac2_ens,kac3_ens,kac4_ens
   ! intermediates for matrix of direction cosine

!! LPO calculation along the streamline

   DO WHILE (step2 .GT. 0)

!! Lij and time from record
      l = stream_Lij(:,:,step2)
      t_loc = stream_dt(step2)
      epsnot = stream_e(step2)
!! strain rate tensor
      e(1,1) = l(1,1) ; e(3,3) = l(3,3) ; e(2,2) = l(2,2)
      e(3,1) = (l(1,3)+l(3,1))/2d0 ; e(1,3) = e(3,1)
      e(1,2) = (l(2,1)+l(1,2))/2d0 ; e(2,1) = e(1,2)
      e(2,3) = (l(3,2)+l(2,3))/2d0 ; e(3,2) = e(2,3)

!! time stepping for LPO calculation
      dt = MIN(t_loc,1d-2/epsnot)
!! number of iterations in the LPO loop
      N_strain = NINT(t_loc/dt)
      dt=t_loc/dfloat(N_strain)

!! LPO loop on the point on the pathline

      DO nn = 1 , N_strain

      fsei = fse
      odfi = odf ; acsi = acs
      odfi_ens = odf_ens ; acsi_ens = acs_ens

      CALL deriv(alpha, setDiffusionCreep)

      kfse1 = MATMUL(l,fsei)*dt
      kodf1 = dotodf*dt*epsnot
      kodf1_ens = dotodf_ens*dt*epsnot
      kac1 = dotacs*dt*epsnot
      kac1_ens = dotacs_ens*dt*epsnot

      fsei = fse + 0.5d0*kfse1
      odfi = odf + 0.5d0*kodf1
      acsi = acs + 0.5d0*kac1
      odfi_ens = odf_ens + 0.5d0*kodf1_ens
      acsi_ens = acs_ens + 0.5d0*kac1_ens

      DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
         IF (acsi(j,j1,j2) .GT. 1d0) acsi(j,j1,j2) = 1d0
         IF (acsi_ens(j,j1,j2) .GT. 1d0) acsi_ens(j,j1,j2) = 1d0
         IF (acsi(j,j1,j2) .LT. -1d0) acsi(j,j1,j2) = -1d0
         IF (acsi_ens(j,j1,j2) .LT. -1d0) acsi_ens(j,j1,j2) = -1d0
      END DO ; END DO ; END DO
      DO j = 1 , size
         IF (odfi(j) .LE. 0 ) odfi(j) = 0d0
         IF (odfi_ens(j) .LE. 0 ) odfi_ens(j) = 0d0
      END DO
      odfi = odfi/SUM(odfi)
      odfi_ens = odfi_ens/SUM(odfi_ens)

      CALL deriv(alpha, setDiffusionCreep)

      kfse2 = MATMUL(l,fsei)*dt
      kodf2 = dotodf*dt*epsnot
      kodf2_ens = dotodf_ens*dt*epsnot
      kac2 = dotacs*dt*epsnot
      kac2_ens = dotacs_ens*dt*epsnot

      fsei = fse + 0.5d0*kfse2
      odfi = odf + 0.5d0*kodf2
      odfi_ens = odf_ens + 0.5d0*kodf2_ens
      acsi = acs + 0.5d0*kac2
      acsi_ens = acs_ens + 0.5d0*kac2_ens

      DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
         IF (acsi(j,j1,j2) .GT. 1d0) acsi(j,j1,j2) = 1d0
         IF (acsi_ens(j,j1,j2) .GT. 1d0) acsi_ens(j,j1,j2) = 1d0
         IF (acsi(j,j1,j2) .LT. -1d0) acsi(j,j1,j2) = -1d0
         IF (acsi_ens(j,j1,j2) .LT. -1d0) acsi_ens(j,j1,j2) = -1d0
      END DO ; END DO ; END DO
      DO j = 1 , size
         IF (odfi(j) .LE. 0 ) odfi(j) = 0d0
         IF (odfi_ens(j) .LE. 0 ) odfi_ens(j) = 0d0
      END DO
      odfi = odfi/SUM(odfi)
      odfi_ens = odfi_ens/SUM(odfi_ens) 

      CALL deriv(alpha, setDiffusionCreep)

      kfse3 = MATMUL(l,fsei)*dt 
      kodf3 = dotodf*dt*epsnot
      kodf3_ens = dotodf_ens*dt*epsnot
      kac3 = dotacs*dt*epsnot
      kac3_ens = dotacs_ens*dt*epsnot

      fsei = fse + kfse3
      odfi = odf + kodf3
      odfi_ens = odf_ens + kodf3_ens
      acsi = acs + kac3
      acsi_ens = acs_ens + kac3_ens

      DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
         IF (acsi(j,j1,j2) .GT. 1d0) acsi(j,j1,j2) = 1d0
         IF (acsi_ens(j,j1,j2) .GT. 1d0) acsi_ens(j,j1,j2) = 1d0
         IF (acsi(j,j1,j2) .LT. -1d0) acsi(j,j1,j2) = -1d0
         IF (acsi_ens(j,j1,j2) .LT. -1d0) acsi_ens(j,j1,j2) = -1d0
      END DO ; END DO ; END DO
      DO j = 1 , size
         IF (odfi(j) .LE. 0 ) odfi(j) = 0d0
         IF (odfi_ens(j) .LE. 0 ) odfi_ens(j) = 0d0
      END DO
      odfi = odfi/SUM(odfi)
      odfi_ens = odfi_ens/SUM(odfi_ens)
   
      CALL deriv(alpha, setDiffusionCreep)

      kfse4 = MATMUL(l,fsei)*dt
      kodf4 = dotodf*dt*epsnot
      kodf4_ens = dotodf_ens*dt*epsnot
      kac4 = dotacs*dt*epsnot
      kac4_ens = dotacs_ens*dt*epsnot

      fse = fse + (kfse1/2d0+kfse2+kfse3+kfse4/2d0)/3d0
      acs = acs + (kac1/2d0+kac2+kac3+kac4/2d0)/3d0
      acs_ens = acs_ens + (kac1_ens/2d0+kac2_ens+kac3_ens+kac4_ens/2d0)/3d0
      odf = odf + (kodf1/2d0+kodf2+kodf3+kodf4/2d0)/3d0
      odf_ens = odf_ens + (kodf1_ens/2d0+kodf2_ens+kodf3_ens+kodf4_ens/2d0)/3d0

      DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
         IF (acs(j,j1,j2) .GT. 1d0) acs(j,j1,j2) = 1d0
         IF (acs_ens(j,j1,j2) .GT. 1d0) acs_ens(j,j1,j2) = 1d0
         IF (acs(j,j1,j2) .LT. -1d0) acs(j,j1,j2) = -1d0
         IF (acs_ens(j,j1,j2) .LT. -1d0) acs_ens(j,j1,j2) = -1d0
      END DO; END DO ; END DO

      !!! Grain-boundary sliding

      odf = odf/SUM(odf)
      odf_ens = odf_ens/SUM(odf_ens)

      do j = 1 , size
         if(odf(j) .le. chi/dble(size)) then
            odf(j) = chi/dble(size)
            acs(j,:,:) = acs0(j,:,:)
         end IF
         if(odf_ens(j) .le. chi/dble(size)) then
            odf_ens(j) = chi/dble(size)
            acs_ens(j,:,:) = acs0(j,:,:)
         end IF
      end do

      odf = odf/sum(odf)
      odf_ens = odf_ens/sum(odf_ens)

      END DO

   step2 = step2-1

   END DO

   RETURN

   END SUBROUTINE strain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine check_veloc - matlab plot of velocity field                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE check_veloc

   USE comvar

   IMPLICIT NONE
   
   INTEGER :: i1,i3 ! loop counters

   OPEN(10,file='veloc.m',status='replace')

   WRITE(10,*)
   WRITE(10,*) "clear"
   WRITE(10,*) "figure"

   WRITE(10,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(10,"(121(1PE12.4))") Ui(1,:,i3)
   END DO
   WRITE(10,*) "];"
   
   WRITE(10,*) "Uz=["
   DO i3 = 1,nx3
      WRITE(10,"(121(1PE12.4))") Ui(3,:,i3)
   END DO
   WRITE(10,*) "];"

   WRITE(10,*) "X=["
   WRITE(10,"(121(1PE12.4))") X1
   WRITE(10,*) "];"

   WRITE(10,*) "Z=["
   WRITE(10,"(121(1PE12.4))") X3
   WRITE(10,*) "];"

   WRITE(10,*) "quiver(X,Z,Ux,Uz)"
   WRITE(10,*) "axis([0 1200000 0 400000])"
   WRITE(10,*) "xlabel('X1, dimensionless distance from the ridge axis)')"
   WRITE(10,*) "ylabel('X3, dimensionless depth')"
   WRITE(10,*) "title('Velocity field')"

   WRITE(10,*) "axis image"
   WRITE(10,*) "axis ij"

   CLOSE(10)

   RETURN

   END SUBROUTINE check_veloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine pipar - Calculates GOL parameter at grid point              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE pipar(i1,i3,X10,X30)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i1,i3
   ! coordinates of the grid point

   DOUBLE PRECISION :: X10,X30,X1i,X3i
   ! initial and intermediate position of the point

   DOUBLE PRECISION :: U1i,U3i
   ! local velocity vector

   DOUBLE PRECISION, DIMENSION(3) :: veloc,isa
   ! local velocity vector and orientation of infinite strain axis

   DOUBLE PRECISION :: dt
   ! time step used to calculate GOL

   DOUBLE PRECISION :: thetaISA
   ! angle between ISA and flow direction

!!! RMQ if GOL not defined then GOL = -1d0 - Max value of GOL set to 10
   GOL(i1,i3) = 1d1

   veloc = 0d0 ; isa = 0d0

   CALL velocitycalc(X10,X30,U1i,U3i,i1,i3)
   dt = 1d-2/SQRT(U1i**2+U3i**2)

!!! previous point on the streamline

   X1i = X10-dt*U1i
   X3i = X30-dt*U3i

   CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
   veloc(1) = U1i ; veloc(3) = U3i ; veloc = veloc/SQRT(SUM(veloc**2))

   CALL gradientcalc(X1i,X3i,i1,i3)
   l = l/epsnot

!!! calculation of the ISA
   CALL isacalc(isa,i1,i3)
   IF (SUM(isa) .EQ. -3d0) isa=veloc

!!! angle between ISA and flow direction
   thetaISA = ACOS(SUM(veloc*isa))

!!! next point on the streamline

   X1i = X10+dt*U1i
   X3i = X30+dt*U3i

   CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
   veloc(1) = U1i ; veloc(3) = U3i ; veloc = veloc/SQRT(SUM(veloc**2))

   CALL gradientcalc(X1i,X3i,i1,i3)
   l = l/epsnot

!!! calculation of the ISA
   CALL isacalc(isa,i1,i3)
   IF (SUM(isa) .EQ. -3d0) isa=veloc

!!! Pi parameter
   CALL gradientcalc(X10,X30,i1,i3)
   thetaISA = ABS(thetaISA-ACOS(SUM(veloc*isa)))
   IF (thetaISA .GT. ACOS(-1d0)) thetaISA = thetaISA-ACOS(-1d0)

   GOL(i1,i3) = MIN(GOL(i1,i3),thetaISA/2d0/dt/epsnot)

   RETURN

   END SUBROUTINE pipar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine isacalc - Calculates ISA orienation at grid point           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE isacalc(isa,i1,i3)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i1,i3
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
      isa = 0d0 ; GOL(i1,i3)=-1d0
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine  eigen - find 3 eigen values of velocity gradient tensor    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE eigen(F)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: a2,a3,Q,R,theta,xx,lambda1,lambda2,lambda3

   DOUBLE PRECISION, DIMENSION(3,3) :: F,Id

   Id = 0d0 ; Id(1,1) = 1d0 ; Id(2,2) = 1d0 ; Id(3,3) = 1d0

!!! looking for the eigen-values of L (using tr(l)=0)

   a2 = l(1,1)*l(2,2)+l(2,2)*l(3,3)+l(3,3)*l(1,1) &
       -l(1,2)*l(2,1)-l(2,3)*l(3,2)-l(3,1)*l(1,3)
   a3 = l(1,1)*l(2,3)*l(3,2)+l(1,2)*l(2,1)*l(3,3)+l(1,3)*l(2,2)*l(3,1) &
       -l(1,1)*l(2,2)*l(3,3)-l(1,2)*l(2,3)*l(3,1)-l(1,3)*l(2,1)*l(3,2)

   Q = -a2/3d0
   R = a3/2d0

   IF (ABS(Q) .LT. 1d-9) THEN
      !! simple shear, isa=veloc
      F = -1d0
   ELSE IF (Q**3-R**2 .GE. 0) THEN
         theta = ACOS(R/Q**1.5)
         lambda1 = -2d0*SQRT(Q)*COS(theta/3d0)
         lambda2 = -2d0*SQRT(Q)*COS((theta+2d0*ACOS(-1d0))/3d0)
         lambda3 = -2d0*SQRT(Q)*COS((theta+4d0*ACOS(-1d0))/3d0)

         IF (ABS(lambda1-lambda2) .LT. 1e-13) lambda1=lambda2
         IF (ABS(lambda2-lambda3) .LT. 1e-13) lambda2=lambda3
         IF (ABS(lambda3-lambda1) .LT. 1e-13) lambda3=lambda1

         IF (lambda1 .GT. lambda2 .AND. lambda1 .GT. lambda3) F=matmul(l-lambda2*Id,l-lambda3*Id)
         IF (lambda2 .GT. lambda3 .AND. lambda2 .GT. lambda1) F=matmul(l-lambda3*Id,l-lambda1*Id)
         IF (lambda3 .GT. lambda1 .AND. lambda3 .GT. lambda2) F=matmul(l-lambda1*Id,l-lambda2*Id)

         IF (lambda1 .EQ. lambda2 .AND. lambda3 .GT. lambda1) F=matmul(l-lambda1*Id,l-lambda2*Id)
         IF (lambda2 .EQ. lambda3 .AND. lambda1 .GT. lambda2) F=matmul(l-lambda2*Id,l-lambda3*Id)
         IF (lambda3 .EQ. lambda1 .AND. lambda2 .GT. lambda3) F=matmul(l-lambda3*Id,l-lambda1*Id)

         IF (lambda1 .EQ. lambda2 .AND. lambda3 .LT. lambda1) F=0d0
         IF (lambda2 .EQ. lambda3 .AND. lambda1 .LT. lambda2) F=0d0
         IF (lambda3 .EQ. lambda1 .AND. lambda2 .LT. lambda3) F=0d0
 
      ELSE

         xx = (SQRT(R**2-Q**3)+ABS(R))**(1d0/3d0)
         lambda1 = -SIGN(1d0,R)*(xx+Q/xx)
         lambda2 = -lambda1/2d0
         lambda3 = -lambda1/2d0
         IF (lambda1 .GT. 1d-9) THEN
            F=matmul(l-lambda2*Id,l-lambda3*Id)
         ELSE
            F = 0d0
         END IF

   END IF
   
   RETURN

   END SUBROUTINE eigen
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine DERIV, calculation of the rotation vector and slip rate     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE deriv(alpha, setDiffusionCreep)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,i1,i2,i3,i4,j,k ! counters
   INTEGER :: imax,iint,imin,iinac ! dummies
   INTEGER , DIMENSION(1) :: ti ! reordering array
   INTEGER :: alpha !variable set depending on deformation regime

   LOGICAL :: setDiffusionCreep !true if under diffusion creep regime

   DOUBLE PRECISION :: Emean,rt1,rt2,rt3,Emean_ens,rt0_ens
   !! surface averaged aggregate NRJ
   !! dislocation density for each slip system

   DOUBLE PRECISION :: gam0
   ! slip rate on the softest slip system

   DOUBLE PRECISION :: R1,R2
   DOUBLE PRECISION :: qint,qmin,sn1,rat
   !!! dummies

   DOUBLE PRECISION, DIMENSION(4) :: bigi,q,qab ! intermediates for G calc

   DOUBLE PRECISION, DIMENSION(4) :: gam
   ! ratios of strain between softest slip system and slip system s for Olivine

   DOUBLE PRECISION, DIMENSION(3) :: rot
   !! rotation rate vector

   DOUBLE PRECISION, DIMENSION(3,3) :: g
   ! slip tensor

   DOUBLE PRECISION, DIMENSION(3,3) :: lx,ex
   ! dimensionless velocity gradient and strain rate tensors

!!! Dimensionless strain rate and velocity gradient tensors
   lx = l/epsnot ; ex = e/epsnot

!!! Plastic deformation + dynamic recrystallization

   DO i=1,size

!!! Calculate invariants e_{pr} T_{pr} for the four slip systems of olivine

   bigi=0d0 ; gam = 0d0 ; g = 0d0

   DO i1 = 1,3 ; DO i2 = 1,3
      bigi(1) = bigi(1)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,2,i2)
      bigi(2) = bigi(2)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,3,i2)
      bigi(3) = bigi(3)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,2,i2)
      bigi(4) = bigi(4)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,1,i2)
   ENDDO ; ENDDO

!!! Quotients I/tau
   q = bigi/tau

!!! Reorder quotients I/tau according to absolute magnitude
   qab = ABS(q)
   ti = MAXLOC(qab) ; imax = ti(1) ; qab(imax)=-1d0
   ti = MAXLOC(qab) ; iint = ti(1) ; qab(iint)=-1d0
   ti = MAXLOC(qab) ; imin = ti(1) ; qab(imin)=-1d0
   ti = MAXLOC(qab) ; iinac = ti(1)

!!! Calculate weighting factors gam_s relative to value gam_i for which
!!! I/tau is largest

   gam(imax)=1d0

   rat = tau(imax)/bigi(imax)
   qint = rat*bigi(iint)/tau(iint)
   qmin = rat*bigi(imin)/tau(imin)
   sn1 = stressexp-1d0

   gam(iint)=qint*(abs(qint))**sn1
   gam(imin)=qmin*(abs(qmin))**sn1
   gam(iinac)=0d0

!!! calculation of G tensor

   DO i1 = 1,3 ; DO i2 = 1,3
      g(i1,i2)=2d0*(gam(1)*acsi(i,1,i1)*acsi(i,2,i2) + &
                    gam(2)*acsi(i,1,i1)*acsi(i,3,i2) + &
                    gam(3)*acsi(i,3,i1)*acsi(i,2,i2) + &
                    gam(4)*acsi(i,3,i1)*acsi(i,1,i2))
   END DO ; END DO

!!! calculation of strain rate on the softest slip system

   R1 = 0d0 ; R2 = 0d0

   DO j= 1 , 3
      i2 = j + 2
      IF (i2 .GT. 3) i2 = i2 - 3

      R1 = R1 - (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j))
      R2 = R2 - (g(j,i2)-g(i2,j))*(lx(j,i2)-lx(i2,j))

      DO k = 1 , 3

         R1 = R1 + 2d0*g(j,k)*g(j,k)
         R2 = R2 + 2d0*lx(j,k)*g(j,k)

      END DO
   END DO

   gam0 = R2/R1
   
!!! dislocation density calculation
   
   rt1=tau(1)**(1.5d0-stressexp)*ABS(gam(1)*gam0)**(1.5d0/stressexp)
   rt2=tau(2)**(1.5d0-stressexp)*ABS(gam(2)*gam0)**(1.5d0/stressexp)
   rt3=tau(3)**(1.5d0-stressexp)*ABS(gam(3)*gam0)**(1.5d0/stressexp)

   rt(i) = rt1*exp(-lambda*rt1**2)+ &
           rt2*exp(-lambda*rt2**2)+ &
           rt3*exp(-lambda*rt3**2)

!!! calculation of the rotation rate:
   if (alpha .eq. 0 .and. setDiffusionCreep .eqv. .true.) then
      call fse_decompose(Fij,rot,ex)
      rot(3) = (lx(2,1)-lx(1,2))*0.5d0 + rot(3)
      rot(2) = (lx(1,3)-lx(3,1))*0.5d0 + rot(2)
      rot(1) = (lx(3,2)-lx(2,3))*0.5d0 + rot(1)
   else
      rot(3) = (lx(2,1)-lx(1,2))*0.5d0-(g(2,1)-g(1,2))*0.5d0*gam0
      rot(2) = (lx(1,3)-lx(3,1))*0.5d0-(g(1,3)-g(3,1))*0.5d0*gam0
      rot(1) = (lx(3,2)-lx(2,3))*0.5d0-(g(3,2)-g(2,3))*0.5d0*gam0
   end if

!!! derivative of the matrix of direction cosine

   dotacs(i,:,:) = 0d0 

   DO i1 = 1 , 3 ; DO i2 = 1 , 3 ; DO i3 = 1 , 3 ; DO i4 = 1 , 3
      dotacs(i,i1,i2)=dotacs(i,i1,i2)+alt(i2,i3,i4)*acsi(i,i1,i4)*rot(i3)
   END DO ; END DO ; END DO ; END DO
   
!!! grain boundary sliding for small grains
   !IF (odfi(i) .LT. chi/REAL(size)) THEN
   !   dotacs(i,:,:) = 0d0
   !   rt(i) = 0d0
   !END IF

   END DO

!!! Volume averaged energy
   Emean = SUM(odfi*rt)

!!! Change of volume fraction by grain boundary migration
   DO i = 1 , size
      dotodf(i) = Xol * Mob * odfi(i) * (Emean-rt(i))
   END DO

!!! ENSTATITE

   DO i=1,size

!!! Calculate slip tensor

   g = 0d0

   DO i1 = 1,3 ; DO i2 = 1,3

      g(i1,i2) = 2d0*acsi_ens(i,3,i1)*acsi_ens(i,1,i2)

   ENDDO ; ENDDO

!!! calculation of strain rate

   R1 = 0d0 ; R2 = 0d0

   DO j= 1 , 3
      i2 = j + 2
      IF (i2 .GT. 3) i2 = i2 - 3

      R1 = R1 - (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j))
      R2 = R2 - (g(j,i2)-g(i2,j))*(lx(j,i2)-lx(i2,j))

      DO k = 1 , 3

         R1 = R1 + 2d0*g(j,k)*g(j,k)
         R2 = R2 + 2d0*lx(j,k)*g(j,k)

      END DO
   END DO

   gam0 = R2/R1

! weight factor between olivine and enstatite

   gam0 = gam0*(1d0/tau_ens)**stressexp

!!! dislocation density calculation

   rt0_ens=tau_ens**(1.5d0-stressexp)*ABS(gam0)**(1.5d0/stressexp)

   rt_ens(i) = rt0_ens*exp(-lambda*rt0_ens**2)
   
!!! calculation of the rotation rate: 
   
   rot(3) = (lx(2,1)-lx(1,2))/2d0-(g(2,1)-g(1,2))/2d0*gam0

   rot(2) = (lx(1,3)-lx(3,1))/2d0-(g(1,3)-g(3,1))/2d0*gam0

   rot(1) = (lx(3,2)-lx(2,3))/2d0-(g(3,2)-g(2,3))/2d0*gam0

   dotacs_ens(i,:,:) = 0d0

   DO i1 = 1 , 3 ; DO i2 = 1 , 3 ; DO i3 = 1 , 3 ; DO i4 = 1 , 3
      dotacs_ens(i,i1,i2)=dotacs_ens(i,i1,i2)+alt(i2,i3,i4)*acsi_ens(i,i1,i4)*rot(i3)
   END DO ; END DO ; END DO ; END DO

!!! grain boundary sliding for small grains
   !IF (odfi_ens(i) .LT. chi/REAL(size)) THEN
   !   dotacs_ens(i,:,:) = 0d0
   !   rt_ens(i) = 0d0
   !END IF

   END DO

!!! Volume averaged energy
   Emean_ens = SUM(odfi_ens*rt_ens)

!!! Change of volume fraction by grain boundary migration
   DO i = 1 , size
      dotodf_ens(i) = (1d0-Xol) * Mob * odfi_ens(i) * (Emean_ens-rt_ens(i))
   END DO

   RETURN

   END SUBROUTINE deriv
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine voigt - Calculates elastic tensor cav_{ijkl} for an olivine !!!
!!! aggregate using Voigt averaging                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE voigt

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,j,k,ll,n,nu,p,q,r,ss

   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C0,Cav2

   C0 = 0d0 ; Cav = 0d0 ; Sav = 0d0

!!! Single-xl elastic tensors c0_{ijkl}
   DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
      C0(i,j,k,ll) = S0(ijkl(i,j),ijkl(k,ll))
   END DO ; END DO ; END DO ; END DO

   DO nu = 1 , size

     Cav2 = 0d0

      DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
         DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
            Cav2(i,j,k,ll) = Cav2(i,j,k,ll) + &
   acs(nu,p,i)*acs(nu,q,j)*acs(nu,r,k)*acs(nu,ss,ll)*C0(p,q,r,ss)
         END DO ; END DO ; END DO ; END DO
         Cav(i,j,k,ll) = Cav(i,j,k,ll) + Cav2(i,j,k,ll)*odf(nu)*Xol
      END DO ; END DO ; END DO ; END DO

   END DO

!!! Enstatite-xl elastic tensors c0_{ijkl}
   DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
      C0(i,j,k,ll) = S0_ens(ijkl(i,j),ijkl(k,ll))
   END DO ; END DO ; END DO ; END DO

   DO nu = 1 , size

      Cav2 = 0d0

      DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
         DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
            Cav2(i,j,k,ll) = Cav2(i,j,k,ll) + &
   acs_ens(nu,p,i)*acs_ens(nu,q,j)*acs_ens(nu,r,k)*acs_ens(nu,ss,ll)*C0(p,q,r,ss)
           END DO ; END DO ; END DO ; END DO
         Cav(i,j,k,ll) = Cav(i,j,k,ll) + Cav2(i,j,k,ll)*odf_ens(nu)*(1d0-Xol)
      END DO ; END DO ; END DO ; END DO

   END DO

!!! Average stiffness matrix

   DO i = 1 , 6 ; DO j = 1 , 6
      Sav(i,j) = Cav(l1(i),l2(i),l1(j),l2(j))
   END DO ; END DO

   RETURN

   END SUBROUTINE voigt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine for elastic decomposition
!  into transverse isotropy
!
!  DECMOD   :   module of global parameters
!  DECSYM   :   decomposition into transverse
!               isotropy symmetry
!
!  Library of independent subroutines
!  for elasticity and tensor
!
!  FULLSYM6 :   upper right 6x6 matrix into
!               full symmetric 6x6 matrix
!  NDELTA   :   Kronecker function
!  TENS4    :   transform 6x6 matrix
!               into 4th order tensor
!  JACOBI   :   get eigenvectors and eigenvalues
!               of a real symmetric matrix
!  EIGSRT   :   order eigenvectors with
!               increasing eigenvalues
!  MAT6     :   transform 4th order tensor
!               into 6x6 matrix
!  PERMUT   :   permutation of index 123
!  ROT4     :   rotation of 4th order tensor
!  fse_decompose : calculates rotation of a material line based on the finite-strain elipsoid
!
!  Library of subroutines
!  for symmetry decomposition
!  (elastic coefficients in GPa)
!
!  SCCA       : form symmetry cartesian system
!  V21D       : form 21D vector from 6x6 matrix (independent)
!  PROJECTI   : transverse isotropy projector   (independent)
!
!****************************************************************
!

      module DECMOD

      double precision, dimension(3,3) :: SCC
      double precision, dimension(6,6) :: CE1
      double precision, dimension(3,3,3,3) :: EL1
      double precision, dimension(21) :: XEC
      double precision :: XN,ANIS

      end module DECMOD

!
!****************************************************************
!

      subroutine DECSYM(CED,PERC,INCLTI,EPSRAD,EPSPERC)

      use DECMOD

      implicit none

      double precision, dimension(6,6) :: CED
      double precision :: PERC,INCLTI,DC5,PI
      double precision, dimension(3) :: TIAXIS
      double precision :: EPSRAD, toprad, botrad !Parameters to calculate the radial anisotropy
      double precision :: EPSPERC, perc1, perc2 !Parameters to calculate the percentage of azimuthal anisotropy

      PI=acos(-1D0)
      CE1=CED
      EL1=0D0
      call FULLSYM6(CE1)

      toprad=1D0/8D0*(CE1(1,1)+CE1(2,2))-1D0/4D0*CE1(1,2)+1D0/2D0*CE1(6,6)
      botrad=1D0/2D0*(CE1(4,4)+CE1(5,5))
      EPSRAD=toprad/botrad

      perc1=1D0/2D0*(CE1(5,5)-CE1(4,4))
      perc2=CE1(5,4)
      EPSPERC=SQRT(perc1**2+perc2**2)
      
      call TENS4(CE1,EL1)
      call V21D(CE1,XEC)
      XN=sqrt(dot_product(XEC,XEC))

      call SCCA
      call PROJECTI(XEC,DC5)

      PERC=(ANIS-DC5)/XN*100D0

      TIAXIS=SCC(3,:)
      TIAXIS=TIAXIS/sqrt(sum(TIAXIS*TIAXIS))
      INCLTI=asin(TIAXIS(3))

      return

      end  subroutine DECSYM

!
!****************************************************************
!

      subroutine FULLSYM6(C)

      implicit none

      double precision, dimension(6,6) :: C

      C(3,2)=C(2,3)
      C(3,1)=C(1,3)
      C(2,1)=C(1,2)

      C(4,1)=C(1,4)
      C(5,1)=C(1,5)
      C(6,1)=C(1,6)
      C(4,2)=C(2,4)
      C(5,2)=C(2,5)
      C(6,2)=C(2,6)
      C(4,3)=C(3,4)
      C(5,3)=C(3,5)
      C(6,3)=C(3,6)

      C(6,5)=C(5,6)
      C(6,4)=C(4,6)
      C(5,4)=C(4,5)

      return

      end subroutine FULLSYM6

!
!****************************************************************
!

      subroutine TENS4(C,C4)

      implicit none

      integer :: i,j,k,l
      integer :: p,q
      integer :: NDELTA
      double precision, dimension(6,6) :: C
      double precision, dimension(3,3,3,3) :: C4

      C4=0D0

      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3

         p=NDELTA(i,j)*i+(1-NDELTA(i,j))*(9-i-j)
         q=NDELTA(k,l)*k+(1-NDELTA(k,l))*(9-k-l)
         C4(i,j,k,l)=C(p,q)

      end do
      end do
      end do
      end do

      end subroutine TENS4

!
!****************************************************************
!

      function NDELTA(i,j)

      implicit none

      integer :: i,j
      integer :: NDELTA

      NDELTA=0
      if (i==j) NDELTA=1

      end function NDELTA

!
!****************************************************************
!

      subroutine JACOBI(a,n,np,d,v,nrot)

      ! Jacobi algorithm for real symmetric matrix
      ! Gives eigenvalues and orthonormalized eigenvectors
      ! Half of the input matrix is crushed

      implicit none

      integer :: n,np,nrot
      integer, parameter :: NMAX=500,IDP=kind(1D0)

      double precision, dimension(np,np) :: a,v
      double precision, dimension(np) :: d
      double precision, dimension(NMAX) :: b,z

      integer :: i,ip,iq,j
      double precision :: c,g,h,s,sm,t,tau,theta,tresh

      do ip=1,n
         do iq=1,n
            v(ip,iq)=0D0
         enddo
         v(ip,ip)=1D0
      enddo
      do ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=0D0
      enddo
      nrot=0
      do i=1,50
         sm=0D0
         do ip=1,n-1
         do iq=ip+1,n
            sm=sm+abs(a(ip,iq))
         enddo
         enddo
         if (sm==0D0) return
         if (i<4) then
            tresh=0.2D0*sm/real(n,IDP)**2D0
         else
            tresh=0D0
         endif
         do ip=1,n-1
         do iq=ip+1,n
            g=100D0*abs(a(ip,iq))
            if ((i>4).and.(abs(d(ip))+                          &
            g==abs(d(ip))).and.(abs(d(iq))+g==abs(d(iq)))) then
               a(ip,iq)=0D0
            else if (abs(a(ip,iq))>tresh) then
               h=d(iq)-d(ip)
               if (abs(h)+g==abs(h)) then
                  t=a(ip,iq)/h
               else
                  theta=0.5D0*h/a(ip,iq)
                  t=1D0/(abs(theta)+sqrt(1D0+theta**2D0))
                  if (theta<0D0) t=-t
               endif
               c=1D0/sqrt(1D0+t**2D0)
               s=t*c
               tau=s/(1D0+c)
               h=t*a(ip,iq)
               z(ip)=z(ip)-h
               z(iq)=z(iq)+h
               d(ip)=d(ip)-h
               d(iq)=d(iq)+h
               a(ip,iq)=0D0
               do j=1,ip-1
                  g=a(j,ip)
                  h=a(j,iq)
                  a(j,ip)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo
               do j=ip+1,iq-1
                  g=a(ip,j)
                  h=a(j,iq)
                  a(ip,j)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo
               do j=iq+1,n
                  g=a(ip,j)
                  h=a(iq,j)
                  a(ip,j)=g-s*(h+g*tau)
                  a(iq,j)=h+s*(g-h*tau)
               enddo
               do j=1,n
                  g=v(j,ip)
                  h=v(j,iq)
                  v(j,ip)=g-s*(h+g*tau)
                  v(j,iq)=h+s*(g-h*tau)
               enddo
               nrot=nrot+1
            endif
         enddo
         enddo
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0D0
         enddo
      enddo
      write(6,'(''Too many iterations in JACOBI'')')

      return

      end subroutine JACOBI

!
!****************************************************************
!

      subroutine EIGSRT(d,v,n,np)

      ! Order eigenvalues and eigenvectors
      ! 1 : max
      ! 2 : mid
      ! 3 : min

      implicit none

      integer :: np,n
      integer :: i,j,k
      double precision, dimension(np) :: d
      double precision, dimension(np,np) :: v
      double precision :: p

      do i=1,n-1
         k=i
         p=d(i)
         do j=i+1,n
            if (d(j)>=p) then
               k=j
               p=d(j)
            end if
         end do
         if (k/=i) then
            d(k)=d(i)
            d(i)=p
            do j=1,n
               p=v(j,i)
               v(j,i)=v(j,k)
               v(j,k)=p
            end do
         end if
      end do

      return

      end subroutine EIGSRT

!
!****************************************************************
!

      subroutine MAT6(C4,C)

      implicit none

      integer :: i
      double precision, dimension(6,6) :: C
      double precision, dimension(3,3,3,3) :: C4

      C = 0D0

      do i=1,3
         C(i,i)=C4(i,i,i,i)
      end do
      do i=2,3
         C(1,i)=(C4(1,1,i,i)+C4(i,i,1,1))/2D0
         C(i,1)=C(1,i)
      end do
      C(2,3)=(C4(2,2,3,3)+C4(3,3,2,2))/2D0
      C(3,2)=C(2,3)

      do i=1,3
         C(i,4)=(C4(i,i,2,3)+C4(i,i,3,2)+                       &
                 C4(2,3,i,i)+C4(3,2,i,i))/4D0
         C(4,i)=C(i,4)
      end do
      do i=1,3
         C(i,5)=(C4(i,i,1,3)+C4(i,i,3,1)+                       &
                 C4(1,3,i,i)+C4(3,1,i,i))/4D0
         C(5,i)=C(i,5)
      end do
      do i=1,3
         C(i,6)=(C4(i,i,1,2)+C4(i,i,2,1)+                       &
                 C4(1,2,i,i)+C4(2,1,i,i))/4D0
         C(6,i)=C(i,6)
      end do

      C(4,4)=(C4(2,3,2,3)+C4(2,3,3,2)+                          &
              C4(3,2,2,3)+C4(3,2,3,2))/4D0
      C(5,5)=(C4(1,3,1,3)+C4(1,3,3,1)+                          &
              C4(3,1,1,3)+C4(3,1,3,1))/4D0
      C(6,6)=(C4(2,1,2,1)+C4(2,1,1,2)+                          &
              C4(1,2,2,1)+C4(1,2,1,2))/4D0
      C(4,5)=(C4(2,3,1,3)+C4(2,3,3,1)+                          &
              C4(3,2,1,3)+C4(3,2,3,1)+                          &
              C4(1,3,2,3)+C4(1,3,3,2)+                          &
              C4(3,1,2,3)+C4(3,1,3,2))/8D0

      C(5,4)=C(4,5)
      C(4,6)=(C4(2,3,1,2)+C4(2,3,2,1)+                          &
              C4(3,2,1,2)+C4(3,2,2,1)+                          &
              C4(1,2,2,3)+C4(1,2,3,2)+                          &
              C4(2,1,2,3)+C4(2,1,3,2))/8D0
      C(6,4)=C(4,6)
      C(5,6)=(C4(1,3,1,2)+C4(1,3,2,1)+                          &
              C4(3,1,1,2)+C4(3,1,2,1)+                          &
              C4(1,2,1,3)+C4(1,2,3,1)+                          &
              C4(2,1,1,3)+C4(2,1,3,1))/8D0
      C(6,5)=C(5,6)

      return

      end subroutine MAT6

!
!****************************************************************
!

      subroutine PERMUT(INDEX,PERM)

      implicit none

      integer :: INDEX
      integer, dimension(3) :: PERM

      if (INDEX==1) then
         PERM(1)=1
         PERM(2)=2
         PERM(3)=3
      end if
      if (INDEX==2) then
         PERM(1)=2
         PERM(2)=3
         PERM(3)=1
      end if
      if (INDEX==3) then 
         PERM(1)=3
         PERM(2)=1
         PERM(3)=2
      endif

      return

      end subroutine PERMUT

!
!****************************************************************
!

      subroutine ROT4(C4,R,C4C)

      implicit none

      integer :: i1,i2,i3,i4,j1,j2,j3,j4
      double precision, dimension(3,3,3,3) :: C4,C4C
      double precision, dimension(3,3) :: R

      C4C = 0D0

      do i1=1,3
      do i2=1,3
      do i3=1,3
      do i4=1,3

         do j1=1,3
         do j2=1,3
         do j3=1,3
         do j4=1,3

            C4C(i1,i2,i3,i4) = C4C(i1,i2,i3,i4) +               &
            R(i1,j1)*R(i2,j2)*R(i3,j3)*R(i4,j4)*C4(j1,j2,j3,j4)

         end do
         end do
         end do
         end do

      end do
      end do
      end do
      end do

      return

      end subroutine ROT4

!
!****************************************************************
!

      subroutine SCCA

      use DECMOD

      implicit none

      integer :: i,NROT,i1,i2,NDVC
      integer, dimension(3) :: IHS
      double precision, dimension(3) :: EGDI,EGVO
      double precision, dimension(3,3) :: DI,VO,VECDI,VECVO
      double precision, dimension(6,6) :: CEC
      double precision, dimension(3,3,3,3) :: ELC
      double precision, dimension(21) :: XH,XD
      double precision :: SDV,ADV,ADVC,SCN,DEV,K,G

      DI=0D0
      VO=0D0
      K=0D0
      G=0D0

      do i=1,3
         DI(1,1)=CE1(1,i)+DI(1,1)
         DI(2,2)=CE1(2,i)+DI(2,2)
         DI(3,3)=CE1(3,i)+DI(3,3)
         DI(2,1)=CE1(6,i)+DI(2,1)
         DI(3,1)=CE1(5,i)+DI(3,1)
         DI(3,2)=CE1(4,i)+DI(3,2)
      end do
      DI(1,2)=DI(2,1)
      DI(1,3)=DI(3,1)
      DI(2,3)=DI(3,2)

      VO(1,1)=CE1(1,1)+CE1(6,6)+CE1(5,5)
      VO(2,2)=CE1(6,6)+CE1(2,2)+CE1(4,4)
      VO(3,3)=CE1(5,5)+CE1(4,4)+CE1(3,3)
      VO(2,1)=CE1(1,6)+CE1(2,6)+CE1(4,5)
      VO(1,2)=VO(2,1)
      VO(3,1)=CE1(1,5)+CE1(3,5)+CE1(4,6)
      VO(1,3)=VO(3,1)
      VO(3,2)=CE1(2,4)+CE1(3,4)+CE1(5,6)
      VO(2,3)=VO(3,2)

      do i=1,3
         K=K+DI(i,i)
         G=G+VO(i,i)
      end do
      K=K/9D0
      G=G/10D0-3D0*K/10D0

      ! Anisotropy

      ANIS=0D0
      XH=0D0
      XD=0D0
      XH(1) = K + 4D0*G/3D0
      XH(2) = K + 4D0*G/3D0
      XH(3) = K + 4D0*G/3D0
      XH(4) = sqrt(2D0)*(K-2D0*G/3D0)
      XH(5) = sqrt(2D0)*(K-2D0*G/3D0)
      XH(6) = sqrt(2D0)*(K-2D0*G/3D0)
      XH(7) = 2D0*G
      XH(8) = 2D0*G
      XH(9) = 2D0*G

      XD=XEC-XH
      ANIS=sqrt(dot_product(XD,XD))

      ! Dil. and Voigt axes

      call JACOBI(DI,3,3,EGDI,VECDI,NROT)
      call EIGSRT(EGDI,VECDI,3,3)
      call JACOBI(VO,3,3,EGVO,VECVO,NROT)
      call EIGSRT(EGVO,VECVO,3,3)

      ! Search for SCCA directions

      do i1=1,3
         NDVC=0
         ADVC=10D0
         SCN=0D0
         do i2=1,3
            SDV=dot_product(VECDI(:,i1),VECVO(:,i2))
            if (abs(SDV)>=1D0) SDV=sign(1D0,SDV)
            ADV=acos(SDV)
            if (SDV<0D0) ADV=acos(-1D0)-ADV
            if (ADV<ADVC) then
               NDVC=sign(1D0,SDV)*i2
               ADVC=ADV
            endif
         end do

         VECDI(:,i1)=(VECDI(:,i1)+NDVC*VECVO(:,abs(NDVC)))/2D0
         SCN=sqrt(VECDI(1,i1)**2D0+VECDI(2,i1)**2D0+            &
                  VECDI(3,i1)**2D0)
         VECDI(:,i1)=VECDI(:,i1)/SCN
      end do

      ! Higher symmetry axis

      SCC = transpose(VECDI)
      ELC = 0D0
      SDV = XN
      do i=1,3
         call PERMUT(i,IHS)
         do i1=1,3
            VECDI(i1,:)=SCC(IHS(i1),:)
         end do
         call ROT4(EL1,VECDI,ELC)
         call MAT6(ELC,CEC)
         call V21D(CEC,XEC)
         call PROJECTI(XEC,DEV)
         if (DEV<SDV) then
             SDV=DEV
             NDVC=i
         endif
      end do

      VECDI=SCC
      call PERMUT(NDVC,IHS)
      do i1=1,3
         SCC(i1,:)=VECDI(IHS(i1),:)
      end do

      ! Rotate in SCCA

      call ROT4(EL1,SCC,ELC)
      call MAT6(ELC,CEC)
      call V21D(CEC,XEC)

      return

      end subroutine SCCA

!
!****************************************************************
!

      subroutine V21D(C,X)

      implicit none

      double precision, dimension(6,6) :: C
      double precision, dimension(21) :: X

      X  = 0D0
      X(1)  = C(1,1)
      X(2)  = C(2,2)
      X(3)  = C(3,3)
      X(4)  = sqrt(2D0)*C(2,3)
      X(5)  = sqrt(2D0)*C(1,3)
      X(6)  = sqrt(2D0)*C(1,2)
      X(7)  = 2D0*C(4,4)
      X(8)  = 2D0*C(5,5)
      X(9)  = 2D0*C(6,6)
      X(10) = 2D0*C(1,4)
      X(11) = 2D0*C(2,5)
      X(12) = 2D0*C(3,6)
      X(13) = 2D0*C(3,4)
      X(14) = 2D0*C(1,5)
      X(15) = 2D0*C(2,6)
      X(16) = 2D0*C(2,4)
      X(17) = 2D0*C(3,5)
      X(18) = 2D0*C(1,6)
      X(19) = 2D0*sqrt(2D0)*C(5,6)
      X(20) = 2D0*sqrt(2D0)*C(4,6)
      X(21) = 2D0*sqrt(2D0)*C(4,5)

      return

      end subroutine V21D

!
!***************************************************************
!

      subroutine PROJECTI(X,DEV)

      implicit none

      double precision :: DEV
      double precision, dimension(21) :: X,XH,XD

      XH=0D0
      XD=0D0
      DEV=0D0

      XH(1)=3D0/8D0*(X(1)+X(2))+X(6)/4D0/sqrt(2D0)+X(9)/4D0
      XH(2)=XH(1)
      XH(3)=X(3)
      XH(4)=(X(4)+X(5))/2D0
      XH(5)=XH(4)
      XH(6)=(X(1)+X(2))/4D0/sqrt(2D0)+3D0/4D0*X(6)              &
            -X(9)/2D0/sqrt(2D0)
      XH(7)=(X(7)+X(8))/2D0
      XH(8)=XH(7)
      XH(9)=(X(1)+X(2))/4D0-X(6)/2D0/sqrt(2D0)+X(9)/2D0

      XD=X-XH
      DEV=sqrt(dot_product(XD,XD))

      return

      end subroutine PROJECTI

!
!***************************************************************
!

      subroutine fse_decompose(fse,rot,Eij)

      use comvar 
      use DECMOD
      implicit none

      integer :: nrot
      ! number of rotations for the Jacobi
      integer :: i1,i2,i3,i4,n
      !! counters
      integer, dimension(3) :: n_2a, n_3a
      double precision :: int1, int2 !! Calculation intermediate
      double precision, dimension(3,3) :: fse,Eij !! finite strain ellipsoid
      ! dimensionless velocity gradient and strain rate tensors
      double precision, dimension(3,3) :: LSij
      ! left-strech tensor for FSE calculation
      double precision, dimension(3) :: evals
      double precision, dimension(3,3) :: evects
      ! eigen values and vectors in jacobi
      
      double precision, dimension(3) :: H
      !! ration H_n of the eigenvalues as in Ribe, 1992
      
      double precision, dimension(3) :: rot
      !! vector of induced spin omega_i, i=1,3
      
      int1 = 0.d0 ; 
      int2 = 0.d0 ; 
      H = 0.d0 ; 
      rot = 0.d0 ;

      n_2a(1) = 2
      n_2a(2) = 3
      n_2a(3) = 1
      
      n_3a(1) = 3
      n_3a(2) = 1
      n_3a(3) = 2

      !!! Left-stretch tensor for FSE calculation
      LSij = matmul(fse,transpose(fse)) !!!!already done above??
      call jacobi(LSij,3,3,evals,evects,nrot)
      
      !! Each column of evects is an eigenvector
      !! Its transpose is the cosines matrix
      evects = transpose(evects);
      
      !!! Calculation of eigenvalue rations
      if (evals(1) .eq. evals(2) .or. evals(1) .eq. evals(3) .or. evals(2) .eq. evals(3)) then
         H=0d0
      else
         H(1)=(evals(1))*(evals(2)-evals(3))/((evals(1)-evals(2))*(evals(1)-evals(3)))
         H(2)=(evals(2))*(evals(3)-evals(1))/((evals(2)-evals(3))*(evals(2)-evals(1)))
         H(3)=(evals(3))*(evals(1)-evals(2))/((evals(3)-evals(1))*(evals(3)-evals(2)))
      end if
      
      !! equation 33 of Ribe 1992 with n=n, i=i1, j=i2, k=i3, l=i4
      do i1 = 1 , 3  
         do n = 1 , 3     
            int2 = 0.d0 ;
            do i2 = 1 , 3 ; do i3 = 1 , 3 ;
               int1 = 0.d0 ;
               do i4 = 1 , 3 ;
                  int1 = int1 + alt(i1,i4,i3)*evects(n_3a(n),i2)*evects(n_3a(n),i4)
               end do
               int2 = int2 + (evects(n,i2)*evects(n_2a(n),i3)*evects(n_3a(n),i1)+int1)*Eij(i2,i3)     
            end do; end do ;
            rot(i1)=rot(i1) + H(n)*int2
         end do
      end do
      
      end subroutine fse_decompose
