module calculate_coulomb_FFT_mod
    use iso_fortran_env, only: real64
    use common_func_mod
    implicit none
    private

    integer, parameter :: dp = real64
    contains

        subroutine calculate_coulomb_long(rho, R_pot, hx, V_long, V_short)
            implicit none
            real(dp), intent(in) :: rho(:,:,:)
            real(dp), intent(in) :: R_pot
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: V_long(:,:,:)
            real(dp), intent(out) :: V_short(:,:,:)
            real(dp), allocatable :: phi_l(:,:,:,:)
            real(dp), allocatable :: rho_l(:,:)
            real(dp), allocatable :: a_lm(:,:)
            real(dp), allocatable :: q_lm(:,:),integral_rho_r(:)
            real(dp) :: r, theta, phi,pi
            integer :: l, m, n, n_middle
            integer :: i, j, k
            integer :: nx, ny, nz
            pi = 4.0_dp * atan(1.0_dp)
            nx = size(rho, 1)
            ny = size(rho, 2)
            nz = size(rho, 3)
            n_middle = (nx + 1) / 2  ! Assuming nx is odd, this gives the middle index
            allocate(rho_l(n_middle:nx, 0:2))
            allocate(a_lm(n_middle:nx, 0:2))
            allocate(q_lm(n_middle:nx, 0:2))
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
            ! 
            !---------------------------



        end subroutine calculate_coulomb_long


        function calculate_rho_l(R_pot, nx, hx, i,j,k) result(rho_l)
            implicit none
            real(dp), intent(in) :: R_pot
            integer, intent(in) :: nx
            real(dp), intent(in) :: hx
            integer, intent(in) :: i, j, k
            real(dp) :: rho_l
            integer :: l, n_middle
            real(dp) :: r
            n_middle = (nx + 1) / 2  ! Assuming nx is odd, this gives the middle index
            rho_l = 0.0_dp
            r = sqrt(real(i - n_middle, dp)**2 + real(j - n_middle, dp)**2 + real(k - n_middle, dp)**2) * hx
            if (r <= R_pot) then
                do l = 0, 2
                    rho_l = 1.0_dp * 2.0_dp * r**2 /(R_pot**2 + r**2)
                end do
            end if
            

        end function calculate_rho_l






    

end module calculate_coulomb_FFT_mod