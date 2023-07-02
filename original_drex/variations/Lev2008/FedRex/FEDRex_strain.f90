!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine STRAIN - Calculation of strain along pathlines              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Changed - Einat Lev, July 2007
!! This subroutine now only calculates the finite strain and ACS/ODF tensors 
!! for the current timestep and for one particle. In the original DREX 
!! code this used to be a loop.
!! Original code - Ed Kaminski, IPG-Paris, 2004


   SUBROUTINE strain(step_counter, particle, fse)

   USE comvar 
   USE GENERAL_VARS

   IMPLICIT NONE

   INTEGER :: j,j1,j2,nn,N_strain 
   ! loop counters

   INTEGER :: step_counter
   ! surrent timestep number

   INTEGER :: particle
   ! number of particle being processed. Needed for Xol value
   
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

!! Get the results of the pervious step for this particle
	acs = acs_particles(particle, :,:,:)
	odf = odf_particles(particle,:)
	
	acs_ens = acs_ens_particles(particle, :,:,:)
	odf_ens = odf_ens_particles(particle,:)

!! Lij and time from record
      l = stream_Lij(:,:,step_counter)
      t_loc = stream_dt(step_counter)
      epsnot = stream_e(step_counter)

!! strain rate tensor
      e(1,1) = l(1,1) ; e(3,3) = l(3,3) ; e(2,2) = l(2,2)
      e(3,1) = (l(1,3)+l(3,1))/2d0 ; e(1,3) = e(3,1)
      e(1,2) = (l(2,1)+l(1,2))/2d0 ; e(2,1) = e(1,2)
      e(2,3) = (l(3,2)+l(2,3))/2d0 ; e(3,2) = e(2,3)

!! time stepping for LPO calculation
      dt = MIN(t_loc,1d-2/epsnot)
!! number of iterations in the LPO loop (nint=nearest integer)
      N_strain = NINT(t_loc/dt)


!! LPO loop on the point on the pathline
  DO nn = 1 , N_strain
      fsei = fse
      odfi = odf ; acsi = acs
      odfi_ens = odf_ens ; acsi_ens = acs_ens
      CALL deriv(particle)

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

      CALL deriv(particle)

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

      CALL deriv(particle)

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

      CALL deriv(particle)

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
      END DO ; END DO ; END DO
 
      odf = odf/SUM(odf)
      odf_ens = odf_ens/SUM(odf_ens)

      END DO !! end N_strain loop

	acs_particles(particle, :,:,:) 		= acs
	odf_particles(particle,:) 		= odf
	acs_ens_particles(particle, :,:,:) 	= acs_ens
	odf_ens_particles(particle,:) 		= odf_ens

   RETURN

   END SUBROUTINE strain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
