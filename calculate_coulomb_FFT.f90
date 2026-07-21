module calculate_coulomb_FFT_mod
    use iso_fortran_env, only: real64
    use constant_mod
    use common_func_mod
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
    private

    integer, parameter :: dp = real64
    contains

        subroutine calculate_coulomb(rho, R_pot, hx, V)
            implicit none
            real(dp), intent(in) :: rho(:,:,:)
            real(dp), intent(in) :: R_pot
            real(dp), intent(in) :: hx
            real(dp), allocatable :: V_long(:,:,:)
            real(dp), allocatable :: rho_short(:,:,:)
            real(dp), allocatable :: rho_long(:,:,:)
            real(dp), allocatable :: V_short(:,:,:)
            real(dp), intent(out) :: V(:,:,:)
            real(dp), allocatable :: phi_l(:,:)
            real(dp), allocatable :: rho_l3(:,:)
            real(dp), allocatable :: rho_l(:,:)
            real(dp), allocatable :: a_lm(:,:)
            real(dp), allocatable :: q_lm(:,:),integral_rho_r(:)
            real(dp), allocatable :: difference(:,:,:)
            real(dp), allocatable :: P(:) ! lagrange multiplier for the linear system
            real(dp) :: r, theta, phi
            real(C_DOUBLE), allocatable :: in(:,:,:)
            complex(C_DOUBLE_COMPLEX), allocatable :: out(:,:,:)
            type(C_PTR) :: plan
            integer :: info
            integer, allocatable :: ipiv(:)
            integer :: l, m, n, n_middle
            integer :: i, j, k, o
            integer :: nx, ny, nz, n_max
            nx = size(rho, 1)
            ny = size(rho, 2)
            nz = size(rho, 3)
            n_middle = (nx + 1) / 2  ! Assuming nx is odd, this gives the middle index
            n_max = nx - n_middle     ! Number of points from the middle to the edge
            allocate(rho_l(0:n_max, 0:2))
            allocate(a_lm(0:2, -2:2))
            allocate(q_lm(0:2, -2:2))
            allocate(integral_rho_r(0:2))
            


            !---------------------------
            ! q_lm
            !---------------------------
            q_lm = 0.0_dp
            do m = -2, 2
                do l = 0, 2
                    do k = 1, nz
                        do j = 1, ny
                            do i = 1, nx
                                r = sqrt(real(i - n_middle, dp)**2 + real(j - n_middle, dp)**2 + real(k - n_middle, dp)**2) * hx
                                theta = acos(real(k - n_middle, dp) * hx / r)
                                phi = atan2(real(j - n_middle, dp) * hx, real(i - n_middle, dp) * hx)
                                q_lm(l,m) = q_lm(l,m) + rho(i,j,k) * spherical_harmonic(l, m, theta, phi) * r**l * hx**3
                            end do
                        end do
                    end do
                end do
            end do

            !---------------------------
            ! a_lm
            !---------------------------
            integral_rho_r = 0.0_dp
            do l = 0, 2
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            r = sqrt(real(i - n_middle, dp)**2 + real(j - n_middle, dp)**2 + real(k - n_middle, dp)**2) * hx
                            integral_rho_r(l) = integral_rho_r(l) + rho(i,j,k) * r**(2*l+1) * hx**3
                        end do
                    end do
                end do
            end do

            a_lm = 0.0_dp
            do m = -2, 2
                do l = 0, 2
                    a_lm(l,m) = q_lm(l,m) / integral_rho_r(l)/4.0_dp / pi
                end do
            end do

            !---------------------------
            ! rho_l
            !---------------------------
            rho_l = 0.0_dp
            do l = 0, 2
               do i = 0, n_max
                    r = real(i, dp) * hx
                    if ( r < R_pot ) then
                        if (l == 0) then
                            rho_l(i,l) = 1
                        else
                            rho_l(i,l) = 2.0_dp * r**l /(R_pot**l + r**l)
                        end if
                    else
                        rho_l(i,l) = 0.0_dp
                    end if
                    rho_l3(i,l) = - rho_l(i,l) * r**3 * 4.0_dp * pi *e2 
               end do
            end do
            
            !---------------------------
            ! phi_l
            !---------------------------
            allocate(phi_l(0:n_max, 0:2))
            allocate(difference(1:n_max, 1:n_max, 0:2))
            phi_l(0,:) = 0.0_dp

            !---------------------------
            ! difference
            !---------------------------
            difference = 0.0_dp
            do l = 0, 2
                do j = 1, n_max
                    do i = 1, n_max
                        r = j * hx
                        if ( i == j ) then
                            difference(i,j,l) = -2 * r**2 / hx**2 - l*(l+1)
                        else if ( i == j - 1 ) then
                            difference(i,j,l) = r**2 / hx**2 
                        else if ( i == j + 1 ) then
                            difference(i,j,l) = r**2 / hx**2 
                        else
                            difference(i,j,l) = 0.0_dp
                        end if
                    end do
                end do
            end do

            !---------------------------
            ! Solve the linear system for phi_l
            !---------------------------
            do l = 0, 2
                call DGESV(n_max, 1, difference(1:n_max, 1:n_max, l), n_max, ipiv, rho_l3(1:n_max, l), n_max, info)
                if (info /= 0) then
                    print *, "Error in solving linear system for phi_l, info =", info
                    stop
                end if
                phi_l(1:n_max, l) = rho_l3(1:n_max, l)
            end do

            !---------------------------
            ! Calculate V_long
            !---------------------------
            allocate(V_long(nx, ny, nz))
            allocate(P(nx))
            P = 0.0_dp

            V_long = 0.0_dp
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        r = sqrt(real(i - n_middle, dp)**2 + real(j - n_middle, dp)**2 + real(k - n_middle, dp)**2) * hx
                        theta = acos(real(k - n_middle, dp) * hx / r)
                        phi = atan2(real(j - n_middle, dp) * hx, real(i - n_middle, dp) * hx)
                        do l = 0, 2
                            do m = -l, l
                                do n = 1, nx
                                    do o = 1, nx
                                        if ( n == o ) then
                                            cycle
                                        else
                                        P(n) = P(n) * (r - o*hx)/(n*hx - o*hx)
                                        end if
                                    end do
                                    P(n) = P(n) * phi_l(n, l)
                                    V_long(i,j,k) = V_long(i,j,k) + a_lm(l,m) * P(n) * spherical_harmonic(l, m, theta, phi)
                                end do 
                            end do
                        end do
                    end do
                end do
            end do
            V = V_long
            deallocate(V_long)
            deallocate(P)
            deallocate(phi_l)
            deallocate(rho_l)
            deallocate(rho_l3)
            deallocate(a_lm)
            deallocate(q_lm)


            !---------------------------
            ! Calculate rho_long
            !---------------------------
            allocate(rho_long(nx, ny, nz))
            rho_long = 0.0_dp
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        r = sqrt(real(i - n_middle, dp)**2 + real(j - n_middle, dp)**2 + real(k - n_middle, dp)**2) * hx
                        theta = acos(real(k - n_middle, dp) * hx / r)
                        phi = atan2(real(j - n_middle, dp) * hx, real(i - n_middle, dp) * hx)
                        do l = 0, 2
                            do m = -l, l
                                rho_long(i,j,k) = rho_long(i,j,k) + a_lm(l,m) * spherical_harmonic(l, m, theta, phi) * 2 * r**l / (R_pot**l + r**l)
                            end do
                        end do
                    end do
                end do
            end do


            !---------------------------
            ! Calculate rho_short
            !---------------------------
            allocate(rho_short(nx, ny, nz))
            rho_short = rho - rho_long  
            deallocate(rho_long)
            allocate(in(nx, ny, nz))
            allocate(out(nx, ny, nz/2 + 1))

            !---------------------------
            ! Perform FFT on rho_short
            !---------------------------
            in = rho_short
            plan = fftw_plan_dft_r2c_3d(nx, ny, nz, in, out, FFTW_ESTIMATE)
            call fftw_execute_dft_r2c(plan, in, out)
            call fftw_destroy_plan(plan)   
            


            


                        


        end subroutine calculate_coulomb
    

end module calculate_coulomb_FFT_mod