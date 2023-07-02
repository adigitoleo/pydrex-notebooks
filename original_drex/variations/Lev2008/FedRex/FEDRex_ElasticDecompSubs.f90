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
!!! Changed - Einat Lev, July 2007
!!! This is for the most part the original code of D-Rex's DECSYM. 
!!! However, the projection part was replaced with Thorsten Becker's modification
!!! which takes into account the decomposition into the various symmetries, 
!!! following the theory of Jules Browaeys.
!!! This enables us to handle cases that regular averaging leads to 
!!! TI axis = slow axis (instead of fast axis).
!!! See also a comment about this in FEDRex_main.f90

      subroutine DECSYM(CED,PERC,INCLTI, Aaxis_out, frac_symm)

      USE DECMOD

      IMPLICIT NONE

      double precision, dimension(6,6) :: CED
      double precision :: PERC,INCLTI,DC5,PI, perc_, inclti_, dev, devo
      DOUBLE PRECISION, DIMENSION(3) :: TIAXIS
      DOUBLE PRECISION, DIMENSION(3) :: Aaxis_out
      
      DOUBLE PRECISION, DIMENSION(21) :: xp, xd
      double precision, dimension(6) :: frac_symm
      
      
      PI=acos(-1D0)
      CE1=CED
      EL1=0D0
      frac_symm = 0d0

      call FULLSYM6(CE1)
      call TENS4(CE1,EL1)
      call V21D(CE1,XEC)
      XN=sqrt(dot_product(XEC,XEC)) !! norm of the whole vector

      call SCCA
      !!call PROJECTI(XEC,DC5)
      
      call project_various(XEC,xp,xd,dev,2)  	!! first get the isotropic projection
      frac_symm(1) = XN - dev 			!! norm of the isotropic projection
      devo = dev				!!devo = non-isotropic part == ANIS
      
      call project_various(xd, xp, xd, dev, 5)  !!get hexagonal projection, dev=non-iso,non-hex
      frac_symm(2) = devo-dev			!!devo-dev = hex only part
      devo = dev				!! devo = npn-iso, non-hex
         
      	perc_ =  (ANIS-dev)/XN*100D0 
      	!! dev is the norm of the non-iso,non-hex part, ANIS is the norm of the non-iso part, so 
      	!! this line calculates: only hex part / whole vector size, times 100
      
      call project_various(XEC, xp, xd, dev, 6) !! get tetra part, dev=non-iso,non-hex,non-tetra
      frac_symm(3) = devo-dev			!! devo-dev = tetra part only
      devo = dev				!! devo = non-iso, non-hex, non-tetra .....
              
      call project_various(XEC, xp, xd, dev, 9) !! get ortho part, dev=non-iso,non-hex,non-tetra,non-ortho
      frac_symm(4) = devo-dev			!! devo-dev = ortho part only
      devo = dev				!! devo = non-iso, non-hex, non-tetra,non-ortho
      
      call project_various(XEC, xp, xd, dev, 13)!! get ortho, dev=non-iso,non-hex,non-tetra,non-ortho,non-mono
      frac_symm(5) = devo-dev			!! devo-dev = mono part only
      devo = dev
      frac_symm(6) = dev			!! devo = non-everything, so only triclinic
      
      TIAXIS=SCC(3,:)
      TIAXIS=TIAXIS/sqrt(sum(TIAXIS*TIAXIS))

      !! INCLTI=asin(TIAXIS(3))
      Aaxis_out = TIaxis
      inclti_ = ATAN(TIAXIS(3)/TIAXIS(1))
      PERC = perc_
      INCLTI = inclti_
      frac_symm = frac_symm/XN*100d0
      
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
! MAT6 - convert from a 4th rank tensor to a 6 by 6 tensor using Kelvin notation
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
! ROT4 - Rotation of a 4th rank tensor
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
      ! XH is the isotropic projection. XD is thus the remaining - anisotropic - component of the vector XEC
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

      !!! Rotate in SCCA
      call ROT4(EL1,SCC,ELC)
      call MAT6(ELC,CEC)
      call V21D(CEC,XEC)
 
      return

      end subroutine SCCA

!
!****************************************************************
!
!
!****************************************************************
!
!  drex_d21v       : form 6x6 Voigt matrix from [21] 
!                    vector(independent)
!

subroutine D21V(x,c)

  implicit none
  ! input
  double precision, intent(in),dimension(21) :: X
  ! output
  double precision, intent(out),dimension(6,6) :: C
  double precision :: one_over_sqrt_two
  parameter (one_over_sqrt_two = 0.70710678118654752440084436210485d0)

  c = 0.0d0

  C(1,1) = x(1)
  C(2,2) = x(2)
  C(3,3) = x(3)
  c(2,3) = X(4)* one_over_sqrt_two 
  c(1,3) = X(5)* one_over_sqrt_two 
  c(1,2) = X(6)* one_over_sqrt_two 
  c(4,4) = x(7)/2d0
  c(5,5) = x(8)/2d0
  c(6,6) = x(9)/2d0
  c(1,4) = x(10)/2d0
  c(2,5) = x(11)/2d0
  c(3,6) = x(12)/2d0
  c(3,4) = x(13)/2d0
  c(1,5) = x(14)/2d0
  c(2,6) = x(15)/2d0
  c(2,4) = x(16)/2d0
  c(3,5) = x(17)/2d0
  c(1,6) = x(18)/2d0
  c(5,6) = x(19)/2d0* one_over_sqrt_two
  c(4,6) = x(20)/2d0* one_over_sqrt_two
  c(4,5) = x(21)/2d0* one_over_sqrt_two
  ! fill in other elements
  call FULLSYM6(c)

  return

end subroutine D21V

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
! compute isotropic projection and elasticity constants
! In TWB's code this is called drex_dec_iso
! Input: X
! output: DEV
!****
      subroutine PROJECTI(X,DEV)!!,XD)

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


!
!  subroutines from D-REX by E. Kaminski as of 01/11/2003
!
!  added code from Jules Browaeys as of 07/2005
!
!  minor changes by TWB
!
!  $Id: drex_decsym.f90,v 1.16 2006/08/02 17:39:33 becker Exp $
!
! Equivalent transverse isotropic media calculated using the theory of   
! Jules Browaeys. 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine for elastic decomposition
!  into transverse isotropy
!
!  drex_decsym   :   decomposition into transverse
!                     isotropy symmetry
!
!  Library of independent subroutines
!  for elasticity and tensor
!

!
!****************************************************************
!
!  DECSYM   :   decomposition into transverse
!               isotropy symmetry
subroutine drex_decsym_ftrn(c, kmod, gmod, symm_frac, tiaxis, &
     dc_iso, dc_hex, dc_tet, dc_ort, dc_mon, dc_tri, vel, &
     dc_orient, scc_rot)
  !
  ! input:  ced(6,6) elastic matrix in Voigt format
  !
  ! output: 
  !
  ! kmod, gmod: K and G moduli of isotropic component
  ! symm_frac[6]
  !
  ! fractional norms for:
  ! 1: aniisotropy total
  ! 2: anisotropy due to hex
  ! 3: anisotropy due to tet
  ! 4: anisotropy due to ort
  ! 5: anisotropy due to mon
  ! 6: anisotropy due to tri
  !
  ! tiaxis(6): normalized symmetry axis vector for the best fit hexagonal (TI) axis
  !            the first three are the best-fit, the last three the 
  !            worst fit (for slow symmetries)
  !
  ! dc_iso, dc_hex, dc_tet, dc_ort, dc_mon, dc_tri:
  !
  ! tensor components as (6,6) Voigt matrices, dc_orient is the c tensor in the hexagonal SCC system
  ! for isotropic, hexagonal, tetragonal, orthorhombic, monoclinic, and triclinic parts
  !
  ! vel(7) velocitiy factors
  !
  ! vel(1) C_11 for ani in hex sytem for p velocities
  ! vel(2) C_33
  !
  ! vel(3) = c_{33}^0,pseudo v_p^2\rho
  ! vel(4) = c_{44}^0 pseudo v_S^2\rho
  !
  !
  ! vel(5) = epsilon*c_33^0 = (C_11-C-33)/2, p wave anisotropy
  ! vel(6) = gamma *c_44^0  = (C_66-C_44)/2, s wave anisotropy
  ! vel(7) = delta*c_33^0   = (C_13 - C_33+2 C_44), ellipticity
  !
  ! vel(8), vel(9): vs1*sqrt(rho), vs2*sqrt(rho)
  !
  implicit none
  
  ! input
  double precision, intent(in),dimension(6,6) :: c
  
  ! output
  double precision, intent(out) :: kmod, gmod ! isotropic tensors
  double precision, intent(out),dimension(9) :: vel 

  double precision, intent(out),dimension(6) :: symm_frac 

  double precision, intent(out),dimension(6) :: tiaxis ! VTI axis, best, worst
  double precision, intent(out),dimension(6,6) :: dc_iso,dc_hex,dc_tet,dc_ort,dc_mon,&
       dc_tri,dc_orient
  double precision, intent(out), dimension (3,3) :: scc_rot ! rotation matrix
  !
  ! local
  double precision, dimension(21) :: x,xp,xd,xiso,xhex,xani
  double precision, dimension(6,6) :: ced_k
  double precision :: xn2,anis,xn,dev,devo,dev13,dev9,dev6,dev5,dev2
  double precision, dimension(3,3) :: di, vo,vvo,vdi
  double precision, dimension(3) :: dval, vval
  ! functions 
  double precision, dimension (6,6) :: drex_ct4,drex_kelc
  double precision, dimension (3,3,3,3) :: drex_t4c
  ! init
  symm_frac = 0d0
  !
  !
  call FULLSYM6(c)! fill in symmetric elements in [6,6] tensor
  !
  ! get the contraction of the Voigt matrix
  !
  call DEC_CONT(c,di,dval,vdi,vo,vval,vvo)
  !
  ! convert the Voigt matrix to [21] format
  call V21D(c,x)
  !
  ! norm of total tensor
  !
  xn = (sqrt(dot_product(x,x)))
  !
  ! get isotropic component xiso[21] and anisotropic vector xd[21]
  ! symm frac will be the anisotropic norm (of xd[])
  ! 
  call drex_dec_iso(x,xiso,xd,symm_frac(1),kmod,gmod,di,vo)
  
  call D21V(xiso,dc_iso)
  !
  ! rotate into Scc system, this takes the voigt matrix
  !
  call drex_dec_scca(c,vdi,vvo,scc_rot,x)
  ! original tensor in rotated frame
  call D21V(x,dc_orient)
  !
  ! fast VTI axis
  !
  tiaxis(1:3)=scc_rot(3,:)          !best fit 
  tiaxis(4:6)=scc_rot(1,:)          !worst fit
  ! normalize
  tiaxis(1:3)=tiaxis(1:3)/sqrt(sum(tiaxis(1:3)*tiaxis(1:3)))
  tiaxis(4:6)=tiaxis(4:6)/sqrt(sum(tiaxis(4:6)*tiaxis(4:6)))
  !
  ! remove isotropic (which does not depend on this rotation)
  xani = x - xiso               ! anisotropic tensor in hex system
  !
  ! hexagonal projection
  !
  call drex_dec_proj(xani,xp,xd,dev,5)
  symm_frac(2) = symm_frac(1) - dev
  call D21V(xp,dc_hex)
  xhex = xiso + xp ! best fit hexagonal = isotropic + hex component
  !
  ! use best-fit hexagonal to compute eps, gamma, delta factors
  !
  !
  ! fast and slow P velocities from hexagonal system tensor
  !
  ! normally, vp2 > vp1, but could be reversed if the symmetry axis 
  ! is slow
  !
  vel(1) = sqrt(xhex(3))                  ! vp1*sqrt(dens)
  vel(2) = sqrt(xhex(1))                  ! vp2*sqrt(dens)
  
  !
  ! shear wave
  !
  vel(8)=sqrt(xhex(7)/2.0d0)           ! vs1 = sqrt(c_hex(4,4)/dens)
  vel(9)=sqrt(xhex(9)/2.0d0)           ! vs2 = sqrt(c_hex(6,6)/dens)
  ! 
  ! reference values, based on best-fit hex
  !
  ! mean of c11 and c11, c33ref, P velocity
  ! (better to use C^iso_11)
  !
  vel(3) = (xhex(1)+xhex(3))/2.0d0 ! c_33^0 = (c_11+c_33)/2; c_11=x(1);   c_33=x(3)
  ! mean of c44 and c66, S velocity
  ! (better to use C^iso_44
  !
  vel(4) = (xhex(7)+xhex(9))/4.0d0 ! c_44^0 = (c_44+c_66)/2; c_44=x(7)/2; c_66=x(9)/2
 
  !
  ! epsilon, gamma, delta
  ! 
  ! epsilon = (C_11-C_33)/2C_33^0; vel(3) is eps*c_33^0
  ! gamma   = (C_66-C_44)/2C_44^0; vel(4) is gamma*c_44^0
  ! delta   = (C_13-C_33+2 C_44)/(C_33^0); C_13 = x(5)/sqrt(2); vel(5)=delta*c_33^0
  !
  vel(5) = (xhex(1)-xhex(3))/2.0d0 ! eps
  vel(6) = (xhex(9)-xhex(7))/4.0d0 ! gamma
  vel(7) = (xhex(5)/sqrt(2.0d0)-xhex(3)+xhex(7)) ! delta
  !
  ! tetragonal fraction
  devo = dev
  call drex_dec_proj(xd,xp,xd,dev,6)
  symm_frac(3) = devo - dev
  call D21V(xp,dc_tet)
  ! orthorhombic
  devo=dev
  call drex_dec_proj(xd,xp,xd,dev,9)
  symm_frac(4) = devo-dev
  call D21V(xp,dc_ort)
  ! monoclinic
  devo=dev
  call drex_dec_proj(xd,xp,xd,dev,13)
  symm_frac(5) = devo-dev
  call D21V(xp,dc_mon)
  ! triclinic
  symm_frac(6) = dev
  call D21V(xd,dc_tri)

  !
  ! normalize the norms, these are now symmetry fractions
  symm_frac = symm_frac/xn
  
  !call drex_dec_proj(x,xp,xd,dev13,13) !monoclinic
  !call drex_dec_proj(x,xp,xd,dev9,9) !orth
  !call drex_dec_proj(x,xp,xd,dev6,6) !tetra
  !call drex_dec_proj(x,xp,xd,dev5,5) !hex
  !call drex_dec_proj(x,xp,xd,dev2,2) !iso
  !write(*,*) ' monoclinic % =',(dev9-dev13)/xn
  !write(*,*) ' orthorhombic % =',(dev6-dev9)/xn
  !write(*,*) ' tetragonal % =',(dev5-dev6)/xn
  !write(*,*) ' hexagonal % =',(dev2-dev5)/xn
  !write(*,*) ' anisotropic % =',dev2/xn

end  subroutine drex_decsym_ftrn




!
!****************************************************************
!
!  EIGSRT   :   order eigenvectors with
!               increasing eigenvalues

subroutine drex_eigsrt(d,v,n,np)

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

end subroutine DREX_EIGSRT

!
!****************************************************************
!
!  PERMUT   :   permutation of index 123

subroutine PERMUT(INDEX,PERM)

  implicit none

  integer,intent(in) :: INDEX
  integer, intent(out), dimension(3) :: PERM

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



subroutine DEC_CONT(c,di,dval,vdi,vo,vval,vvo)
  !
  ! Calculates contractions of elastic cij voigt matrix
  ! i.e. dilatationnal and Voigt tensors
  ! and produces eigenvectors and eigenvalues
  !
  !use drex_nrmod
  implicit none
  integer :: i,nrot
  double precision, intent(out),dimension(3) :: dval,vval
  double precision, intent(out),dimension(3,3) :: di,vo,vdi,vvo
  double precision, intent(in), dimension(6,6) :: C
  double precision, dimension(3,3) :: dicopy,vocopy
  !
  ! Dilatationnal 3x3 tensor
  !
  DI=0D0
  do i=1,3
     DI(1,1)=C(1,i)+DI(1,1)
     DI(2,2)=C(2,i)+DI(2,2)
     DI(3,3)=C(3,i)+DI(3,3)
     DI(2,1)=C(6,i)+DI(2,1)
     DI(3,1)=C(5,i)+DI(3,1)
     DI(3,2)=C(4,i)+DI(3,2)
  end do
  DI(1,2)=DI(2,1)
  DI(1,3)=DI(3,1)
  DI(2,3)=DI(3,2)
  !Voigt 3x3 tensor
  VO=0D0
  VO(1,1)=C(1,1)+C(6,6)+C(5,5)
  VO(2,2)=C(6,6)+C(2,2)+C(4,4)
  VO(3,3)=C(5,5)+C(4,4)+C(3,3)
  VO(2,1)=C(1,6)+C(2,6)+C(4,5)
  VO(1,2)=VO(2,1)
  VO(3,1)=C(1,5)+C(3,5)+C(4,6)
  VO(1,3)=VO(3,1)
  VO(3,2)=C(2,4)+C(3,4)+C(5,6)
  VO(2,3)=VO(3,2)
  ! Diagonalization
  ! make copies, jacobi is destructive
  dicopy=di
  vocopy=vo
  ! get eigensystem
  
      call JACOBI(di,3,3,dval,VDI,nrot)
      call EIGSRT(dval,VDI,3,3)
      call JACOBI(vocopy,3,3,vval,VVO,nrot)
      call EIGSRT(vval,VVO,3,3)
  
end subroutine DEC_CONT

!
! compute isotropic projection and elasticity constants
!
subroutine drex_dec_iso(x,xp,xd,dev,kiso,giso,di,vo)
  implicit none
  !Isotropic parameters and anisotropic percentage
  integer :: i
  double precision, intent(in),dimension(21) :: x
  double precision, intent(in), dimension(3,3) :: di,vo
  double precision, intent(out) :: kiso,giso,dev
  double precision, intent(out),dimension(21) :: xp,xd
  double precision :: vsum, dsum
  !
  ! isotropic parameters
  !
  vsum=0d0
  dsum=0d0
  do i=1,3
     dsum=dsum+di(i,i)
     vsum=vsum+vo(i,i)
  end do
  kiso=dsum/9d0
  giso=0.1d0*vsum - dsum/30.0d0

  !Isotropic projection vector
  xp=0d0
  xp(1)=KISO+4D0*GISO/3D0
  XP(2)=xp(1)
  XP(3)=xp(1)
  XP(4)=sqrt(2D0)*(KISO-2D0*GISO/3D0)
  XP(5)=xp(4)
  XP(6)=xp(4)
  XP(7)=2D0*GISO
  XP(8)=xp(7)
  XP(9)=xp(7)
  ! difference 
  xd = x - xp
  dev = (sqrt(dot_product(xd,xd)))

  
end subroutine DREX_DEC_ISO
!
! get scc system, rotate, and obtain hex fraction
!
subroutine drex_dec_scca(c,vdi,vvo,scc,x21r)
  implicit none
  !
  ! Calculates SCC directions as bissectrix
  ! of each close pair of dilatationnal and
  ! Voigt eigenvectors
  ! Higher symmetry axis = 3-axis
  !
  !t   = angles VO(1st index)/DI(2nd index)
  !bss = bissectrix direction intermediate
  !SCT = SCC test
  !dvm = array deviation for permutations
  !
  double precision, parameter :: pi=3.14159265358979323846d0
  !
  ! input voigt matrix
  !
  double precision, intent(in),dimension(6,6) :: c
  !
  ! input eigensystem of contractions
  !
  double precision, intent(in),dimension(3,3) :: vdi,vvo
  ! output
  double precision, intent(out),dimension(3,3) :: scc
  double precision, intent(out),dimension(21) :: x21r
  ! local
  integer :: i,j
  integer, dimension(1) :: pos
  integer, dimension(3) :: npos
  double precision :: dev
  double precision, dimension(3) :: bss,dvm
  double precision, dimension(3,3) :: t,SCT
  double precision, dimension(6,6) :: cr
  double precision, dimension(21) :: xp,xd
  !
  !VO/DI eigenvectors angle matrix
  !
  do i=1,3
     do j=1,3
        t(j,i) = dot_product(vdi(:,i),vvo(:,j))
     end do
  end do
  
  where (abs(t)>=1D0) t=sign(1D0,t)
  t=acos(t)
  where (t>0.5D0*pi) t=t-pi
  !Indice position
  !Association VO eigenvector to DI
  npos=minloc(abs(t),dim=1)
  !
  ! Calculates bissectrix vector & SCC
  !
  do i=1,3
     bss=VDI(:,i)+sign(1D0,t(npos(i),i))*vvo(:,npos(i))
     scc(:,i)=bss/(sqrt(dot_product(bss,bss)))
  end do
  !Transpose : basis transfer
  scc=transpose(scc)
  !Permutation 1,2,3
  npos=(/1,2,3/)
  !
  ! rotate c and convert to [21]
  call ROT6(c,scc,cr)
  call V21D(cr,x21r)
  ! project 
  call drex_dec_proj(x21r,xp,xd,dvm(1),5)
  !
  !Permutations 2,3,1 & 3,1,2
  !
  do j=2,3
     npos=cshift(npos,shift=1)
     do i=1,3
        SCT(i,:)=SCC(npos(i),:)
     end do
     call ROT6(c,sct,cr)
     call V21D(cr,x21r)
     ! project this permutation
     call drex_dec_proj(x21r,xp,xd,dvm(j),5)
  end do
  !
  ! choose the permutation for minimum deviation
  !
  pos = minloc(dvm)
  ! rest of transverse isotropic fraction
  dev = dvm(pos(1))
  !SCC coordinates system
  scc=cshift(scc,shift=pos(1)-1,dim=1)
  !
  !21-D elastic vector in SCC axes
  !
  call ROT6(c,scc,cr)
  call V21D(cr,x21r)
end subroutine drex_dec_scca




!
! compute symmetry projections
!
subroutine drex_dec_proj(x,xp,xd,dev,nsym)
  implicit none
  !
  !Symmetry projectors
  !
  !X    =  input 21-D vector
  !Xp   =  projected vector
  !XD   =  deviation vector
  !xdn  =  norm of deviation
  !NSYM = 13 monoclinic
  !     =  9 orthorhombic
  !     =  6 tetragonal
  !     =  5 hexagonal (transverse isotropic)
  !     =  2 isotropic
  integer, intent(in) :: nsym
  double precision, intent(in),dimension(21) :: x
  double precision, intent(out),dimension(21) :: xp,xd
  double precision, intent(out) :: dev
  ! local
  double precision :: sq2,isq2,i15
  parameter (sq2 = 1.41421356237309504d0) ! sqrt(2)
  isq2=1D0/sq2
  i15=1D0/15D0
  xp=0d0
  if (NSYM==13) then            ! monoclinic
     XP=X
     XP(10:11)=0D0
     XP(13:14)=0D0
     XP(16:17)=0D0
     XP(19:20)=0D0
  endif
  if (NSYM<=9) then
     XP(1:9)=X(1:9)
  endif
  if (NSYM==6) then
     XP(1)=0.5D0*(X(1)+X(2))
     XP(2)=XP(1)
     XP(4)=0.5D0*(X(4)+X(5))
     XP(5)=XP(4)
     XP(7)=0.5D0*(X(7)+X(8))
     XP(8)=XP(7)
  endif
  if (NSYM==5) then
     XP(1)=0.375D0*(X(1)+X(2))+0.25D0*X(6)*isq2+0.25D0*X(9)
     XP(2)=XP(1)
     XP(4)=0.5D0*(X(4)+X(5))
     XP(5)=XP(4)
     XP(6)=0.25D0*(X(1)+X(2))*isq2+0.75D0*X(6)-0.5D0*X(9)*isq2
     XP(7)=0.5D0*(X(7)+X(8))
     XP(8)=XP(7)
     XP(9)=0.25D0*(X(1)+X(2))-0.5D0*X(6)*isq2+0.5D0*X(9)
  endif
  if (NSYM==2) then
     XP(1)=3D0*(X(1)+X(2)+X(3))+sq2*(X(4)+X(5)+X(6))
     XP(1)=XP(1)+2D0*(X(7)+X(8)+X(9))
     XP(1)=XP(1)*i15
     XP(2)=XP(1)
     XP(3)=XP(1)
     XP(4)=sq2*(X(1)+X(2)+X(3))+4D0*(X(4)+X(5)+X(6))
     XP(4)=XP(4)-sq2*(X(7)+X(8)+X(9))
     XP(4)=XP(4)*i15
     XP(5)=XP(4)
     XP(6)=XP(4)
     XP(7)=2D0*(X(1)+X(2)+X(3))-sq2*(X(4)+X(5)+X(6))
     XP(7)=XP(7)+3D0*(X(7)+X(8)+X(9))
     XP(7)=XP(7)*i15
     XP(8)=XP(7)
     XP(9)=XP(7)
  endif
  xd = x - xp                       ! deviation
  ! norm of rest of projection
  dev = (sqrt(dot_product(xd,xd)))
end subroutine DREX_DEC_PROJ




!
!
! rotate a sav(6,6) Voigt matrix by rotation matrix rot
!
!
subroutine ROT6(sav,rot,savr)
  implicit none
  double precision, intent(in), dimension(6,6) :: sav
  double precision, intent(out), dimension(6,6) :: savr
  double precision, intent(in), dimension(3,3) :: rot
  double precision, dimension(3,3,3,3) :: cv,cvr
  
  call FULLSYM6(sav)        ! fill in symmetric elements in [6,6] matrix
  call TENS4(sav,cv)       !convert to 4th order tensor
  call ROT4(cv,rot,cvr)         !rotate 4th order tensor
  call MAT6(cvr,savr)        !convert 4th order tensor into 6,6 voigt

end subroutine ROT6
 
 
 
 !
! compute symmetry projections
!
subroutine project_various(x,xp,xd,dev,nsym)
  implicit none
  !
  !Symmetry projectors
  !
  !X    =  input 21-D vector
  !Xp   =  projected vector
  !XD   =  deviation vector
  !xdn  =  norm of deviation
  !NSYM = 13 monoclinic
  !     =  9 orthorhombic
  !     =  6 tetragonal
  !     =  5 hexagonal (transverse isotropic)
  !     =  2 isotropic
  integer, intent(in) :: nsym
  double precision, intent(in),dimension(21) :: x
  double precision, intent(out),dimension(21) :: xp,xd
  double precision, intent(out) :: dev
  ! local
  double precision :: sq2,isq2,i15
  parameter (sq2 = 1.41421356237309504d0) ! sqrt(2)
  isq2=1D0/sq2
  i15=1D0/15D0
  xp=0d0
  if (NSYM==13) then            ! monoclinic
     XP=X
     XP(10:11)=0D0
     XP(13:14)=0D0
     XP(16:17)=0D0
     XP(19:20)=0D0
  endif
  if (NSYM<=9) then
     XP(1:9)=X(1:9)
  endif
  if (NSYM==6) then
     XP(1)=0.5D0*(X(1)+X(2))
     XP(2)=XP(1)
     XP(4)=0.5D0*(X(4)+X(5))
     XP(5)=XP(4)
     XP(7)=0.5D0*(X(7)+X(8))
     XP(8)=XP(7)
  endif
  if (NSYM==5) then
     XP(1)=0.375D0*(X(1)+X(2))+0.25D0*X(6)*isq2+0.25D0*X(9)
     XP(2)=XP(1)
     XP(4)=0.5D0*(X(4)+X(5))
     XP(5)=XP(4)
     XP(6)=0.25D0*(X(1)+X(2))*isq2+0.75D0*X(6)-0.5D0*X(9)*isq2
     XP(7)=0.5D0*(X(7)+X(8))
     XP(8)=XP(7)
     XP(9)=0.25D0*(X(1)+X(2))-0.5D0*X(6)*isq2+0.5D0*X(9)
  endif
  if (NSYM==2) then
     XP(1)=3D0*(X(1)+X(2)+X(3))+sq2*(X(4)+X(5)+X(6))
     XP(1)=XP(1)+2D0*(X(7)+X(8)+X(9))
     XP(1)=XP(1)*i15
     XP(2)=XP(1)
     XP(3)=XP(1)
     XP(4)=sq2*(X(1)+X(2)+X(3))+4D0*(X(4)+X(5)+X(6))
     XP(4)=XP(4)-sq2*(X(7)+X(8)+X(9))
     XP(4)=XP(4)*i15
     XP(5)=XP(4)
     XP(6)=XP(4)
     XP(7)=2D0*(X(1)+X(2)+X(3))-sq2*(X(4)+X(5)+X(6))
     XP(7)=XP(7)+3D0*(X(7)+X(8)+X(9))
     XP(7)=XP(7)*i15
     XP(8)=XP(7)
     XP(9)=XP(7)
  endif
  xd = x - xp                       ! deviation
  ! norm of REST of projection
  dev = (sqrt(dot_product(xd,xd)))
return
end subroutine PROJECT_VARIOUS
