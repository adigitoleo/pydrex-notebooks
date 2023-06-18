!!! Fortran program based on the DRex code to produce the data plotted in figure 5
!!! of the 2001 paper: https://doi.org/10.1016%2Fs0012-821x%2801%2900356-9
!!! This version does not rely on reading any external files.
!!! The jacobi algorithm for calculating eigenvalues/eigenvectors was copied directly.
!!! The strain algorithm was modified to use fixed timesteps and no location information.
!!! The deriv algorithm was modified to remove code for enstatite.
!!! FIXME: Fix dimension mismatch issue in dotodf L488 and add output writing.
!!! FIXME: We can use the `parameter` attribute to assign constants and avoid some hassle.

module parameters
    implicit none

    ! Constant velocity gradient with dv_y/dx = 2.
    real, dimension(3,3) :: velocity_gradient = reshape( &
        [0,2,0,0,0,0,0,0,0], shape(velocity_gradient) &
    )

    ! Number of grains in the polycrystal, and its cubed root (both integers).
    integer :: n_grains_3sq = 15  ! Actual n_grains = n_grains_3sq ** 3
    integer :: n_grains = 3375

    ! Timestamps for CPO integration.
    integer, private :: i
    integer :: n_steps = 200
    integer, dimension(200) :: timestamps = [(i, i=1, 200)]

    ! Olivine stiffness tensor.
    real, dimension(6,6) :: olivine_stiffness = reshape( &
        [ &
            320.71, 69.84, 71.22, 0.0, 0.0, 0.0, &
            69.84, 197.25, 74.8, 0.0, 0.0, 0.0, &
            71.22, 74.8, 234.32, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 63.77, 0.0, 0.0, &
            0.0, 0.0, 0.0, 0.0, 77.67, 0.0, &
            0.0, 0.0, 0.0, 0.0, 0.0, 78.36 &
        ], shape(olivine_stiffness) &
    )

    ! 3x3x3 permutation symbol (Levi-Civita).
    real, dimension(3,3,3) :: permutation_symbol = reshape( &
        [0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0], &
        shape(permutation_symbol) &
    )

    ! Olivine A-type CRSS values.
    double precision, dimension(4) :: crss = [1d0, 2d0, 3d0, 1d60]

    ! DRex parameters for pure olivine.
    real :: stress_exponent = 3.5
    real, dimension(3) :: gbm_mobility = [0, 50, 200]
    real :: gbs_threshold = 0
    real :: nucleation_efficiency = 5

end module parameters


module outputs
    double precision, dimension(:,:,:), allocatable :: deformation_gradients
    double precision, dimension(:,:,:,:), allocatable :: orientations
    double precision, dimension(:,:), allocatable :: fractions
end module outputs


module local_arrays
    double precision, dimension(3,3) :: strain_rate
    double precision, dimension(3) :: evals
    double precision, dimension(3,3) :: evects
    double precision :: max_strain_rate
    double precision, dimension(3,3) :: fsei, kfse1, kfse2, kfse3, kfse4
    double precision, dimension(:), allocatable :: odfi, kodf1, kodf2, kodf3, kodf4, dotodf
    double precision, dimension(:,:,:), allocatable :: acsi, kac1, kac2, kac3, kac4, dotacs
    double precision, dimension(:), allocatable :: rt
end module local_arrays


program drex_2001_test
    use parameters
    use outputs
    implicit none

    ! Allocate output arrays.
    allocate(deformation_gradients(3,3,n_steps))
    allocate(orientations(n_grains,3,3,n_steps))
    allocate(fractions(n_grains,n_steps))

    ! Initial deformation gradient is the identity.
    deformation_gradients(:,:,1) = 0d0
    deformation_gradients(1,1,1) = 1d0
    deformation_gradients(2,2,1) = 1d0
    deformation_gradients(3,3,1) = 1d0

    ! Initial orientations are random.
    call init_random_orientations

    ! Initial fractions are uniform.
    fractions(:,n_grains) = 1d0/n_grains

    ! Main routine to calculate strain and update mineral orientations.
    call strain
end program drex_2001_test


subroutine init_random_orientations
    !!! Sets the orientation matrices for the first timestep to random orientations.
    use parameters, only: n_grains_3sq, n_grains
    use outputs, only: orientations
    implicit none

    ! Intermediate values and euler angles for orientation initialisation.
    double precision, dimension(:), allocatable :: rand, xe1, xe2, xe3
    double precision :: phi1, theta, phi2

    integer ri, r1, r2, r3
    allocate(rand(3*n_grains), xe1(n_grains), xe2(n_grains), xe3(n_grains))
    call random_number(rand)

    ri = 1
    do r1=1, n_grains_3sq; do r2=1, n_grains_3sq; do r3=1, n_grains_3sq
        xe1(ri) = (real(r1)-rand(ri))/real(n_grains_3sq)*acos(-1d0)
        xe2(ri) = acos(-1d0 + (real(r2)-rand(n_grains+ri))/real(n_grains_3sq)*2d0)
        xe3(ri) = (real(r3)-rand(ri+2*n_grains))/real(n_grains_3sq)*acos(-1d0)
        ri = ri + 1
    end do; end do; end do

    do ri=1, n_grains
        phi1 = xe1(ri) ; theta = xe2(ri) ; phi2 = xe3(ri)

        ! Construct direction cosine matrices.
        orientations(ri,1,1,1)=cos(phi2)*cos(phi1)-cos(theta)*sin(phi1)*sin(phi2)
        orientations(ri,1,2,1)=cos(phi2)*sin(phi1)+cos(theta)*cos(phi1)*sin(phi2)
        orientations(ri,1,3,1)=sin(phi2)*sin(theta)

        orientations(ri,2,1,1)=-sin(phi2)*cos(phi1)-cos(theta)*sin(phi1)*cos(phi2)
        orientations(ri,2,2,1)=-sin(phi2)*sin(phi1)+cos(theta)*cos(phi1)*cos(phi2)
        orientations(ri,2,3,1)=cos(phi2)*sin(theta)

        orientations(ri,3,1,1)=sin(theta)*sin(phi1)
        orientations(ri,3,2,1)=-sin(theta)*cos(phi1)
        orientations(ri,3,3,1)=cos(theta)
    end do
end subroutine init_random_orientations


subroutine strain
    use parameters
    use outputs
    use local_arrays
    implicit none

    integer :: i, j, j1, j2, nrot
    real :: dt

    allocate( &
        odfi(n_grains), kodf1(n_grains), kodf2(n_grains), kodf3(n_grains), kodf4(n_grains), &
        dotodf(n_grains) &
    )
    allocate( &
        acsi(n_grains,3,3), kac1(n_grains,3,3), kac2(n_grains,3,3), kac3(n_grains,3,3), &
        kac4(n_grains,3,3), dotacs(n_grains,3,3) &
    )
    allocate(rt(n_grains))

    ! Constant strain rate calculated from constant velocity gradient.
    strain_rate = (velocity_gradient + transpose(velocity_gradient)) / 2
    call jacobi(strain_rate,3,3,evals,evects,nrot)
    max_strain_rate = maxval(abs(evals))

    do i=1, n_steps-1
        dt = timestamps(i+1) - timestamps(i)

        fsei = deformation_gradients(:,:,i)
        odfi = fractions(:,i) ; acsi = orientations(:,:,:,i)

        call deriv

        kfse1 = matmul(velocity_gradient,fsei)*dt
        kodf1 = dotodf*dt*max_strain_rate
        kac1 = dotacs*dt*max_strain_rate

        fsei = deformation_gradients(:,:,i+1) + 0.5d0*kfse1
        odfi = fractions(:,i+1) + 0.5d0*kodf1
        acsi = orientations(:,:,:,i+1) + 0.5d0*kac1

        do j = 1 , n_grains ; do j1 = 1 , 3 ; do j2 = 1, 3
            if (acsi(j,j1,j2) .gt. 1d0) acsi(j,j1,j2) = 1d0
            if (acsi(j,j1,j2) .lt. -1d0) acsi(j,j1,j2) = -1d0
        end do ; end do ; end do
        do j = 1 , n_grains
            if (odfi(j) .le. 0 ) odfi(j) = 0d0
        end do
        odfi = odfi/sum(odfi)

        call deriv

        kfse2 = matmul(velocity_gradient,fsei)*dt
        kodf2 = dotodf*dt*max_strain_rate
        kac2 = dotacs*dt*max_strain_rate

        fsei = deformation_gradients(:,:,i+1) + 0.5d0*kfse1
        odfi = fractions(:,i+1) + 0.5d0*kodf1
        acsi = orientations(:,:,:,i+1) + 0.5d0*kac1

        do j = 1 , n_grains ; do j1 = 1 , 3 ; do j2 = 1, 3
            if (acsi(j,j1,j2) .gt. 1d0) acsi(j,j1,j2) = 1d0
            if (acsi(j,j1,j2) .lt. -1d0) acsi(j,j1,j2) = -1d0
        end do ; end do ; end do
        do j = 1 , n_grains
            if (odfi(j) .le. 0 ) odfi(j) = 0d0
        end do
        odfi = odfi/sum(odfi)

        call deriv

        kfse3 = matmul(velocity_gradient,fsei)*dt
        kodf3 = dotodf*dt*max_strain_rate
        kac3 = dotacs*dt*max_strain_rate

        fsei = deformation_gradients(:,:,i+1) + 0.5d0*kfse1
        odfi = fractions(:,i+1) + 0.5d0*kodf1
        acsi = orientations(:,:,:,i+1) + 0.5d0*kac1

        do j = 1 , n_grains ; do j1 = 1 , 3 ; do j2 = 1, 3
            if (acsi(j,j1,j2) .gt. 1d0) acsi(j,j1,j2) = 1d0
            if (acsi(j,j1,j2) .lt. -1d0) acsi(j,j1,j2) = -1d0
        end do ; end do ; end do
        do j = 1 , n_grains
            if (odfi(j) .le. 0 ) odfi(j) = 0d0
        end do
        odfi = odfi/sum(odfi)

        call deriv

        kfse4 = matmul(velocity_gradient,fsei)*dt
        kodf4 = dotodf*dt*max_strain_rate
        kac4 = dotacs*dt*max_strain_rate

        deformation_gradients(:,:,i+1) = (kfse1/2d0+kfse2+kfse3+kfse4/2d0)/3d0
        orientations(:,:,:,i+1) = (kac1/2d0+kac2+kac3+kac4/2d0)/3d0
        fractions(:,i+1) = (kodf1/2d0+kodf2+kodf3+kodf4/2d0)/3d0

        do j = 1 , n_grains ; do j1 = 1 , 3 ; do j2 = 1, 3
            if (orientations(j,j1,j2,i+1) .gt. 1d0) orientations(j,j1,j2,i+1) = 1d0
            if (orientations(j,j1,j2,i+1) .lt. -1d0) orientations(j,j1,j2,i+1) = -1d0
        end do ; end do ; end do

        fractions(:,i+1) = fractions(:,i+1)/sum(fractions(:,i+1))
    end do
end subroutine strain


subroutine jacobi(a,n,np,d,v,nrot)
    ! jacobi algorithm for real symmetric matrix
    ! gives eigenvalues and orthonormalized eigenvectors
    ! half of the input matrix is crushed
    implicit none

    integer :: n,np,nrot
    integer, parameter :: nmax=500,idp=kind(1d0)

    double precision, dimension(np,np) :: a,v
    double precision, dimension(np) :: d
    double precision, dimension(nmax) :: b,z

    integer :: i,ip,iq,j
    double precision :: c,g,h,s,sm,t,tau,theta,tresh

    do ip=1,n
        do iq=1,n
            v(ip,iq)=0d0
        enddo
        v(ip,ip)=1d0
    enddo
    do ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0d0
    enddo
    nrot=0
    do i=1,50
        sm=0d0
        do ip=1,n-1
            do iq=ip+1,n
                sm=sm+abs(a(ip,iq))
            enddo
        enddo
        if (sm==0d0) return
        if (i<4) then
            tresh=0.2d0*sm/real(n,idp)**2d0
        else
            tresh=0d0
        endif
        do ip=1,n-1
            do iq=ip+1,n
                g=100d0*abs(a(ip,iq))
                if ((i>4).and.(abs(d(ip))+                          &
                    g==abs(d(ip))).and.(abs(d(iq))+g==abs(d(iq)))) then
                    a(ip,iq)=0d0
                else if (abs(a(ip,iq))>tresh) then
                    h=d(iq)-d(ip)
                    if (abs(h)+g==abs(h)) then
                        t=a(ip,iq)/h
                    else
                        theta=0.5d0*h/a(ip,iq)
                        t=1d0/(abs(theta)+sqrt(1d0+theta**2d0))
                        if (theta<0d0) t=-t
                    endif
                    c=1d0/sqrt(1d0+t**2d0)
                    s=t*c
                    tau=s/(1d0+c)
                    h=t*a(ip,iq)
                    z(ip)=z(ip)-h
                    z(iq)=z(iq)+h
                    d(ip)=d(ip)-h
                    d(iq)=d(iq)+h
                    a(ip,iq)=0d0
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
            z(ip)=0d0
        enddo
    enddo
    write(6,'(''too many iterations in jacobi'')')

    return
end subroutine jacobi


subroutine deriv
    use parameters
    use local_arrays
    implicit none

    integer :: i,i1,i2,i3,i4,j,k ! counters
    integer :: imax,iint,imin,iinac ! dummies
    integer , dimension(1) :: ti ! reordering array

    double precision :: emean,rt1,rt2,rt3
    !! surface averaged aggregate nrj
    !! dislocation density for each slip system

    double precision :: gam0
    ! slip rate on the softest slip system

    double precision :: r1,r2
    double precision :: qint,qmin,sn1,rat
    !!! dummies

    double precision, dimension(4) :: bigi,q,qab ! intermediates for g calc

    double precision, dimension(4) :: gam
    ! ratios of strain between softest slip system and slip system s for olivine

    double precision, dimension(3) :: rot
    !! rotation rate vector

    double precision, dimension(3,3) :: g
    ! slip tensor

    double precision, dimension(3,3) :: lx,ex
    ! dimensionless velocity gradient and strain rate tensors

    !!! dimensionless strain rate and velocity gradient tensors
    lx = velocity_gradient/max_strain_rate ; ex = strain_rate/max_strain_rate

    !!! plastic deformation + dynamic recrystallization
    do i=1,n_grains

    !!! calculate invariants e_{pr} t_{pr} for the four slip systems of olivine

    bigi=0d0 ; gam = 0d0 ; g = 0d0

    do i1 = 1,3 ; do i2 = 1,3
        bigi(1) = bigi(1)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,2,i2)
        bigi(2) = bigi(2)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,3,i2)
        bigi(3) = bigi(3)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,2,i2)
        bigi(4) = bigi(4)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,1,i2)
    enddo ; enddo

    !!! quotients I/crss
    q = bigi/crss

    !!! reorder quotients I/crss according to absolute magnitude
    qab = abs(q)
    ti = maxloc(qab) ; imax = ti(1) ; qab(imax)=-1d0
    ti = maxloc(qab) ; iint = ti(1) ; qab(iint)=-1d0
    ti = maxloc(qab) ; imin = ti(1) ; qab(imin)=-1d0
    ti = maxloc(qab) ; iinac = ti(1)

    !!! calculate weighting factors gam_s relative to value gam_i for which
    !!! i/crss is largest
    gam(imax)=1d0

    rat = crss(imax)/bigi(imax)
    qint = rat*bigi(iint)/crss(iint)
    qmin = rat*bigi(imin)/crss(imin)
    sn1 = stress_exponent-1d0

    gam(iint)=qint*(abs(qint))**sn1
    gam(imin)=qmin*(abs(qmin))**sn1
    gam(iinac)=0d0
    !!! calculation of g tensor

    do i1 = 1,3 ; do i2 = 1,3
        g(i1,i2)=2d0*(gam(1)*acsi(i,1,i1)*acsi(i,2,i2) + &
                        gam(2)*acsi(i,1,i1)*acsi(i,3,i2) + &
                        gam(3)*acsi(i,3,i1)*acsi(i,2,i2) + &
                        gam(4)*acsi(i,3,i1)*acsi(i,1,i2))
    end do ; end do
    !!! calculation of strain rate on the softest slip system

    r1 = 0d0 ; r2 = 0d0

    do j= 1 , 3
        i2 = j + 2
        if (i2 .gt. 3) i2 = i2 - 3
        r1 = r1 - (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j))
        r2 = r2 - (g(j,i2)-g(i2,j))*(lx(j,i2)-lx(i2,j))

        do k = 1 , 3
            r1 = r1 + 2d0*g(j,k)*g(j,k)
            r2 = r2 + 2d0*lx(j,k)*g(j,k)
        end do
    end do

    gam0 = r2/r1

    !!! dislocation density calculation
    rt1=crss(1)**(1.5d0-stress_exponent)*abs(gam(1)*gam0)**(1.5d0/stress_exponent)
    rt2=crss(2)**(1.5d0-stress_exponent)*abs(gam(2)*gam0)**(1.5d0/stress_exponent)
    rt3=crss(3)**(1.5d0-stress_exponent)*abs(gam(3)*gam0)**(1.5d0/stress_exponent)

    rt(i) = rt1*exp(-nucleation_efficiency*rt1**2)+ &
            rt2*exp(-nucleation_efficiency*rt2**2)+ &
            rt3*exp(-nucleation_efficiency*rt3**2)
    !!! calculation of the rotation rate:

    rot(3) = (lx(2,1)-lx(1,2))/2d0-(g(2,1)-g(1,2))/2d0*gam0
    rot(2) = (lx(1,3)-lx(3,1))/2d0-(g(1,3)-g(3,1))/2d0*gam0
    rot(1) = (lx(3,2)-lx(2,3))/2d0-(g(3,2)-g(2,3))/2d0*gam0
    !!! derivative of the matrix of direction cosine

    dotacs(i,:,:) = 0d0

    do i1 = 1 , 3 ; do i2 = 1 , 3 ; do i3 = 1 , 3 ; do i4 = 1 , 3
        dotacs(i,i1,i2)=dotacs(i,i1,i2)+permutation_symbol(i2,i3,i4)*acsi(i,i1,i4)*rot(i3)
    end do ; end do ; end do ; end do

    !!! grain boundary sliding for small grains
    if (odfi(i) .lt. gbs_threshold/real(n_grains)) then
        dotacs(i,:,:) = 0d0
        rt(i) = 0d0
    end if

    end do

    !!! volume averaged energy
    emean = sum(odfi*rt)

    !!! change of volume fraction by grain boundary migration
    do i = 1 , n_grains
        dotodf(i) = gbm_mobility * odfi(i) * (emean-rt(i))
    end do

    return

end subroutine deriv
