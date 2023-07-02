
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine DERIV, calculation of the rotation vector and slip rate     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE deriv (particle)

   USE comvar
   USE GENERAL_VARS

   IMPLICIT NONE

   INTEGER :: i,i1,i2,i3,i4,j,k ! counters
   INTEGER :: imax,iint,imin,iinac ! dummies
   INTEGER :: particle		! number of particle being processed
   INTEGER , DIMENSION(1) :: ti ! reordering array

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

       DO j = 1 , 3
          i2 = j + 2
          IF (i2 .GT. 3) i2 = i2 - 3

          R1 = R1 - (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j))
          R2 = R2 - (g(j,i2)-g(i2,j))*(lx(j,i2)-lx(i2,j))

          DO k = 1 , 3
             R1 = R1 + 2d0*g(j,k)*g(j,k)
            R2 = R2 + 2d0*lx(j,k)*g(j,k)
          END DO !!i2
       END DO   !!j2

     gam0 = R2/R1
   
    !!! dislocation density calculation
   
      rt1=tau(1)**(1.5d0-stressexp)*ABS(gam(1)*gam0)**(1.5d0/stressexp)
      rt2=tau(2)**(1.5d0-stressexp)*ABS(gam(2)*gam0)**(1.5d0/stressexp)
      rt3=tau(3)**(1.5d0-stressexp)*ABS(gam(3)*gam0)**(1.5d0/stressexp)

      rt(i) = rt1*exp(-lambda*rt1**2)+ &
           rt2*exp(-lambda*rt2**2)+ &
           rt3*exp(-lambda*rt3**2)

    !!! calculation of the rotation rate:

     rot(3) = (lx(2,1)-lx(1,2))/2d0-(g(2,1)-g(1,2))/2d0*gam0
     rot(2) = (lx(1,3)-lx(3,1))/2d0-(g(1,3)-g(3,1))/2d0*gam0
     rot(1) = (lx(3,2)-lx(2,3))/2d0-(g(3,2)-g(2,3))/2d0*gam0

    !!! derivative of the matrix of direction cosine
       dotacs(i,:,:) = 0d0 

       DO i1 = 1 , 3 ; 
            DO i2 = 1 , 3 ; 
                DO i3 = 1 , 3 ; 
                    DO i4 = 1 , 3
                      dotacs(i,i1,i2)=dotacs(i,i1,i2)+alt(i2,i3,i4)*acsi(i,i1,i4)*rot(i3)
                    END DO ; 
                END DO 
            END DO 
        END DO

    !!! grain boundary sliding for small grains
     IF (odfi(i) .LT. chi/REAL(size)) THEN
          dotacs(i,:,:) = 0d0
          rt(i) = 0d0
      END IF

   END DO ! end loop over grains i=1,size


!!! Volume averaged energy
   Emean = SUM(odfi*rt)

!!! Change of volume fraction by grain boundary migration
   DO i = 1 , size
      dotodf(i) = Xol(particle) * Mob * odfi(i) * (Emean-rt(i)) 
   END DO


!!! ENSTATITE
!!-------------------

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

!!! weight factor between olivine and enstatite
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
   IF (odfi_ens(i) .LT. chi/REAL(size)) THEN
      dotacs_ens(i,:,:) = 0d0
      rt_ens(i) = 0d0
   END IF

   END DO

!!! Volume averaged energy
   Emean_ens = SUM(odfi_ens*rt_ens)

!!! Change of volume fraction by grain boundary migration
   DO i = 1 , size
      dotodf_ens(i) = (1d0-Xol(particle)) * Mob * odfi_ens(i) * (Emean_ens-rt_ens(i))
   END DO

   RETURN

   END SUBROUTINE deriv
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


