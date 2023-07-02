!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine voigt - Calculates elastic tensor cav_{ijkl} for an olivine !!!
!!! aggregate using Voigt averaging                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE voigt (particle)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,j,k,ll,n,nu,p,q,r,ss, particle

   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C0,Cav2

   acs = acs_particles(particle, :,:,:) 
   odf = odf_particles(particle,:) 		
   acs_ens = acs_ens_particles(particle, :,:,:) 	
   odf_ens = odf_ens_particles(particle,:) 	


  C0 = 0d0 ; Cav = 0d0 ; Sav = 0d0

!!! Single-crystal elastic tensors c0_{ijkl}
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
         Cav(i,j,k,ll) = Cav(i,j,k,ll) + Cav2(i,j,k,ll)*odf(nu)*Xol(particle)
      END DO ; END DO ; END DO ; END DO

  END DO !!! end of loop over grains nu=1,size

!!! Enstatite-crystal elastic tensors c0_{ijkl}
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
         Cav(i,j,k,ll) = Cav(i,j,k,ll) + Cav2(i,j,k,ll)*odf_ens(nu)*(1d0-Xol(particle))
      END DO ; END DO ; END DO ; END DO

   END DO !!! end of loop over grains nu=1,size

!!! Average stiffness matrix
   DO i = 1 , 6 
	DO j = 1 , 6
            Sav(i,j) = Cav(l1(i),l2(i),l1(j),l2(j))
	END DO 
   END DO

   Sav_particles(particle, :,:) = Sav
   RETURN

   END SUBROUTINE voigt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
