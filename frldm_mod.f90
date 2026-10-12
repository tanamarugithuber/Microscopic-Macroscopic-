module frldm_mod
    !$ use omp_lib
    use iso_fortran_env, only: real64
    use constant_mod, only: e2, pi
    use nucleus_mod, only: nucleus_property
    use grid_mod, only: grid_type
    use shape_mod, only: surface_shape, gauss_legendre
    use hockney_mod, only: hockney_coulomb, hockney_yukawa, hockney_exp
    implicit none
    private

    integer, parameter :: dp = real64

    !---------------------------
    ! FRLDM constants
    ! Units:
    !   Energy : MeV
    !   Length : fm
    !---------------------------

    !----------------------------
    ! First category constants
    !----------------------------
    real(dp), parameter, public :: M_H = 7.289034_dp ! MeV, mass excess of hydrogen atom
    real(dp), parameter, public :: M_n = 8.071431_dp ! MeV, mass excess of neutron

    !----------------------------
    ! Second category constants
    !----------------------------
    real(dp), parameter, public :: a_Yukawa = 0.68_dp ! fm, range of Yukawa potential
    real(dp), parameter, public :: r_p = 0.80_dp ! fm, proton root square radius
    real(dp), parameter, public :: a_den = 0.70_dp ! fm, diffuseness of nuclear density distribution
    real(dp), parameter, public :: a_el = 1.433e-5_dp ! MeV, electronic-binding constant

    !----------------------------
    ! Third category constants
    !----------------------------
    real(dp), parameter, public :: r_mac = 4.80_dp ! MeV, average pairing gap constant
    real(dp), parameter, public :: h_np = 6.6_dp ! MeV, neutron-proton interaction strength
    real(dp), parameter, public :: Wigner = 30.0_dp ! MeV, strength of Wigner term
    real(dp), parameter, public :: a_d = 0.90_dp ! Wigner damping constant
    
    !----------------------------
    ! Fourth category constants
    !----------------------------
    real(dp), parameter, public :: a_v = 16.022835d0 ! volume energy constant in MeV
    real(dp), parameter, public :: k_v = 1.927910d0  ! volume asymmetry constant in MeV
    real(dp), parameter, public :: a_s = 21.26946d0  ! surface energy constant in MeV
    real(dp), parameter, public :: k_s = 2.388587d0  ! surface asymmetry constant in MeV
    real(dp), parameter, public :: a_0 = 2.649971d0  ! a0 energy constant in MeV
    real(dp), parameter, public :: c_a = 0.055673d0  ! charge asymmetry constant in MeV

    !----------------------------
    ! Variables for FRLDM calculations
    !----------------------------
    type :: frldm_variables
        real(dp) :: I_n
        real(dp) :: c_1 
        real(dp) :: c_4 
        real(dp) :: f_kfrp 
        real(dp) :: k_f
        real(dp) :: delta_pn_ave ! neutron-proton pairing correction
        real(dp) :: delta_n_ave ! neutron pairing correction
        real(dp) :: delta_p_ave ! proton pairing correction
        real(dp) :: B_s
        real(dp) :: B_3
        real(dp) :: B_1
        ! yukawa potential 
        real(dp), allocatable :: B_1pot_y(:)
        real(dp), allocatable :: B_3pot_y(:)
        real(dp), allocatable :: B_1pot_y3d(:,:,:)
        real(dp), allocatable :: B_3pot_y3d(:,:,:)
        ! coulomb potential
        real(dp), allocatable :: B_3pot_c(:)
        real(dp), allocatable :: B_3pot_c3d(:,:,:)
        ! exponential factors for the surface and Coulomb energy
        real(dp), allocatable :: B_1exp(:,:,:)
        real(dp), allocatable :: B_3exp(:,:,:)


        real(dp) :: B_w
        logical :: verbose = .true.    ! print the terms in calculate_frldm_energy
        real(dp) :: E_frldm            ! MeV, macroscopic contribution to the mass excess
        ! terms of Eq. (62), MeV
        real(dp) :: E_mass = 0.0_dp, E_v = 0.0_dp, E_s = 0.0_dp, E_a0 = 0.0_dp, E_c = 0.0_dp
        real(dp) :: E_cex = 0.0_dp, E_pff = 0.0_dp, E_ca = 0.0_dp, E_W = 0.0_dp, E_pair = 0.0_dp, E_el = 0.0_dp

        contains
            procedure :: calculate_frldm_energy
            procedure :: calculate_b_1_b_3
            procedure :: calculate_b_1_b_3_surface
            procedure :: calculate_exp
    end type frldm_variables

    public :: frldm_variables

    contains
        subroutine calculate_exp(this, g_mod)
            !$use omp_lib
            implicit none
            class(frldm_variables), intent(inout) :: this
            
            type(grid_type), intent(in) :: g_mod

            integer :: i, j, k, l,m,n,o, n_x, n_y, n_z
            real(dp) :: a_Yukawa_inv = 1.0_dp / a_Yukawa
            real(dp) :: a_den_inv = 1.0_dp / a_den
            real(dp) :: x, y, z, r
            print *, "Calculating exponential factors for the surface and Coulomb energy..."
            

            allocate(this%B_1exp(g_mod%n_x_min*2:g_mod%n_x_max*2, g_mod%n_y_min*2:g_mod%n_y_max*2, g_mod%n_z_min*2:g_mod%n_z_max*2))
            allocate(this%B_3exp(g_mod%n_x_min*2:g_mod%n_x_max*2, g_mod%n_y_min*2:g_mod%n_y_max*2, g_mod%n_z_min*2:g_mod%n_z_max*2))

            
            !$omp parallel do collapse(3) default(none) &
            !$omp shared(this, g_mod, a_Yukawa_inv, a_den_inv) &
            !$omp private(i, j, k, x, y, z, r) schedule(static)
            do k = g_mod%n_z_min*2, g_mod%n_z_max*2
                do j = g_mod%n_y_min*2, g_mod%n_y_max*2
                    do i = g_mod%n_x_min*2, g_mod%n_x_max*2

                        x = real(i, dp) * g_mod%h_x
                        y = real(j, dp) * g_mod%h_y
                        z = real(k, dp) * g_mod%h_z

                        r = sqrt(x*x + y*y + z*z)

                        this%B_1exp(i,j,k) = exp(-r * a_Yukawa_inv)
                        this%B_3exp(i,j,k) = exp(-r * a_den_inv)

                    end do
                end do
            end do
            !$omp end parallel do
            print *, "Exponential factors calculated."
        end subroutine calculate_exp

        ! subroutine calculate_b_1_b_3(this, nucleus,g_mod)
        !     class(frldm_variables), intent(inout) :: this
        !     type(nucleus_property), intent(in) :: nucleus
        !     type(grid_type), intent(in) :: g_mod
            
        !     integer :: i, j, k, l,m,n,o,p,q
        !     integer :: n_max
        !     real(dp) :: denominator1 
        !     real(dp) :: denominator3
        !     real(dp) :: inv_a_den = 1.0_dp / a_den
        !     real(dp), allocatable :: pi_rho(:)
        !     print *, "Calculating B_1 and B_3 for the FRLDM..."

        !     denominator1 = 8.0_dp * pi**2 * nucleus%R0**2 * a_Yukawa**4
        !     denominator3 = 32.0_dp * pi**2 * nucleus%R0**5
        !     allocate(this%B_1pot_y(g_mod%n_points))
        !     allocate(this%B_3pot_y(g_mod%n_points))
        !     allocate(this%B_3pot_c(g_mod%n_points))

        !     ! 4pi*rho
        !     allocate(pi_rho(g_mod%n_points))
        !     pi_rho = g_mod%density_index* (4.0_dp * pi )


        !     allocate(helmholtz_matrix(g_mod%n_x_points**3, g_mod%n_x_points**3))
        !     call initialize_helmholtz_matrix(g_mod%n_x_points, g_mod%h_x, a_Yukawa, helmholtz_matrix)
        !     call CG_method(-helmholtz_matrix, pi_rho, this%B_1pot_y)
        !     deallocate(helmholtz_matrix)

        !     allocate(helmholtz_matrix(g_mod%n_x_points**3, g_mod%n_x_points**3))

        !     call initialize_helmholtz_matrix(g_mod%n_x_points, g_mod%h_x, a_den, helmholtz_matrix)
        !     call CG_method(-helmholtz_matrix, pi_rho, this%B_3pot_y)
        !     deallocate(helmholtz_matrix)

        !     allocate(poisson_matrix(g_mod%n_x_points**3, g_mod%n_x_points**3))
        !     call initialize_poisson_matrix(g_mod%n_x_points, g_mod%h_x, poisson_matrix)
        !     call CG_method(-poisson_matrix, pi_rho, this%B_3pot_c)
        !     deallocate(poisson_matrix)

        !     call this%calculate_exp(g_mod)
        !     this%B_1 = 0.0_dp
        !     this%B_3 = 0.0_dp

        !     do k = 1 , g_mod%n_z_points
        !         do j = 1, g_mod%n_y_points
        !             do i = 1, g_mod%n_x_points
        !                 o = (k-1)*g_mod%n_x_points*g_mod%n_y_points + (j-1)*g_mod%n_x_points + i
        !                 do n = 1, g_mod%n_z_points
        !                     do m = 1, g_mod%n_y_points
        !                         do l = 1, g_mod%n_x_points
        !                             p = (n-1)*g_mod%n_x_points*g_mod%n_y_points + (m-1)*g_mod%n_x_points + l
        !                             this%B_1 = this%B_1 - this%B_1exp(i-l,j-m,k-n) * g_mod%dV**2 * g_mod%density_index(p) &
        !                             * g_mod%density_index(o) 
        !                             this%B_3 = this%B_3 - this%B_3exp(i-l,j-m,k-n) * g_mod%dV**2 * 0.5_dp * inv_a_den &
        !                             * g_mod%density_index(p) * g_mod%density_index(o) 
        !                         end do
        !                     end do
        !                 end do
        !                 this%B_1 = this%B_1 + 2.0_dp * this%B_1pot_y(o)* g_mod%dV * g_mod%density_index(o) * a_Yukawa
        !                 this%B_3 = this%B_3 + (this%B_3pot_c(o) - this%B_3pot_y(o)  )*g_mod%density_index(o)* g_mod%dV
        !             end do
        !         end do
        !     end do


        !     this%B_1 = this%B_1 / denominator1 * nucleus%A**(-2.0_dp / 3.0_dp)
        !     this%B_3 = this%B_3 / denominator3 * 15.0_dp * nucleus%A**(-5.0_dp / 3.0_dp)

        !     deallocate(pi_rho)
        !     deallocate(this%B_1pot_y)
        !     deallocate(this%B_3pot_y)
        !     deallocate(this%B_3pot_c)
        !     deallocate(this%B_1exp)
        !     deallocate(this%B_3exp)
        !     print *, "B_1 and B_3 calculated."

        ! end subroutine calculate_b_1_b_3

        subroutine calculate_b_1_b_3(this, nucleus,g_mod)
            !$ use omp_lib
            implicit none
            class(frldm_variables), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            type(grid_type), intent(in) :: g_mod
            
            integer :: i, j, k, l,m,n,o,p,q
            integer :: n_max
            real(dp) :: denominator1 
            real(dp) :: denominator3
            real(dp) :: a_Yukawa_inv = 1.0_dp / a_Yukawa
            real(dp) :: a_den_inv = 1.0_dp / a_den
            real(dp) :: B1_sum, B3_sum
            real(dp), allocatable :: pi_rho(:)
            real(dp), allocatable :: pi_rho3d(:,:,:)
            real(dp), allocatable :: W_1(:,:,:), W_3(:,:,:)
            print *, "Calculating B_1 and B_3 for the FRLDM..."

            denominator1 = 8.0_dp * pi**2 * nucleus%R0**2 * a_Yukawa**4
            denominator3 = 32.0_dp * pi**2 * nucleus%R0**5
            ! allocate(this%B_1pot_y(g_mod%n_points))
            ! allocate(this%B_3pot_y(g_mod%n_points))
            ! allocate(this%B_3pot_c(g_mod%n_points))
            allocate(this%B_1pot_y3d(g_mod%n_x_points, g_mod%n_y_points, g_mod%n_z_points))
            allocate(this%B_3pot_y3d(g_mod%n_x_points, g_mod%n_y_points, g_mod%n_z_points))
            allocate(this%B_3pot_c3d(g_mod%n_x_points, g_mod%n_y_points, g_mod%n_z_points))




            ! call CG_method_helmholtz_EQ(pi_rho, this%B_1pot_y, g_mod%n_x_points, g_mod%h_x, a_Yukawa_inv)
            ! call CG_method_helmholtz_EQ(pi_rho, this%B_3pot_y, g_mod%n_x_points, g_mod%h_x, a_den_inv)
            ! call CG_method_poisson(g_mod%density_index, this%B_3pot_c, g_mod%n_x_points, g_mod%h_x)

            ! Yukawa and Coulomb potentials of the sharp-surface density by Hockney's method
            ! (free-space FFT convolution, exact boundary condition at infinity).
            ! These replace CG_method_helmholtz_EQ3D / CG_method_poisson3D, which imposed a
            ! Dirichlet condition on the box, and the O(N^2) double sum over B_1exp / B_3exp.
            ! density_index3D is constant in each cell (smooth = .false.).
            !   B_1 kernel: (2 - r/a) exp(-r/a) / (r/a) = 2a exp(-r/a)/r - exp(-r/a)
            !   B_3 kernel: 1/r - (1 + r/(2a)) exp(-r/a)/r = 1/r - exp(-r/a)/r - exp(-r/a)/(2a)
            ! Both parts of each kernel use the same cell-averaged quadrature, so their large
            ! cancellation (the kernels integrate to zero in the bulk) is not spoiled.
            allocate(W_1(g_mod%n_x_points, g_mod%n_y_points, g_mod%n_z_points))
            allocate(W_3(g_mod%n_x_points, g_mod%n_y_points, g_mod%n_z_points))
            call hockney_yukawa(g_mod%density_index3D, g_mod%h_x, a_Yukawa, this%B_1pot_y3d)
            call hockney_exp(g_mod%density_index3D, g_mod%h_x, a_Yukawa, W_1)
            call hockney_yukawa(g_mod%density_index3D, g_mod%h_x, a_den, this%B_3pot_y3d)
            call hockney_exp(g_mod%density_index3D, g_mod%h_x, a_den, W_3)
            call hockney_coulomb(g_mod%density_index3D, g_mod%h_x, this%B_3pot_c3d)

            B1_sum = sum((2.0_dp * a_Yukawa * this%B_1pot_y3d - W_1) * g_mod%density_index3D) * g_mod%dV
            B3_sum = sum((this%B_3pot_c3d - this%B_3pot_y3d - 0.5_dp * a_den_inv * W_3) &
                         * g_mod%density_index3D) * g_mod%dV

            this%B_1 = B1_sum / denominator1 * nucleus%A**(-2.0_dp / 3.0_dp)
            this%B_3 = B3_sum / denominator3 * 15.0_dp * nucleus%A**(-5.0_dp / 3.0_dp)

            ! deallocate(pi_rho)
            ! deallocate(this%B_1pot_y)
            ! deallocate(this%B_3pot_y)
            ! deallocate(this%B_3pot_c)
            deallocate(this%B_1pot_y3d)
            deallocate(this%B_3pot_y3d)
            deallocate(this%B_3pot_c3d)
            deallocate(W_1, W_3)
            print *, "B_1 and B_3 calculated."

        end subroutine calculate_b_1_b_3

        subroutine calculate_frldm_energy(this, nucleus, shape, E_frldm, g_mod)
            !----------------------------
            ! FRLDM macroscopic energy, Eq. (62) of FRDM2012, as a contribution to the
            ! atomic mass excess (it includes M_H Z + M_n N):
            !   E_frldm = M_H Z + M_n N - a_v (1 - k_v I^2) A + a_s (1 - k_s I^2) B_1 A^(2/3) + a_0 B_W
            !           + c_1 Z^2/A^(1/3) B_3 - c_4 Z^(4/3)/A^(1/3) + f(k_F r_p) Z^2/A - c_a (N - Z)
            !           + W (|I| B_W + 1/A [Z = N odd]) + average pairing - a_el Z^2.39
            ! B_1 and B_3 are evaluated for `shape` by double surface integrals
            ! (calculate_b_1_b_3_surface); if g_mod is present they are taken from the grid
            ! instead (Hockney convolution, calculate_b_1_b_3), which is slower and has an O(h)
            ! surface error in B_1. B_s from the same shape (Eq. 70). B_W = 1 (single-body shapes).
            !----------------------------
            implicit none
            class(frldm_variables), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            class(surface_shape), intent(in) :: shape
            real(dp), intent(out) :: E_frldm
            type(grid_type), intent(in), optional :: g_mod
            real(dp) :: Z, N, A, x_0, y_0, B1_0, B3_0

            Z = real(nucleus%Z, dp)
            N = real(nucleus%N, dp)
            A = real(nucleus%A, dp)
            this%I_n = (N - Z) / A
            this%c_1 = 0.6_dp * e2 / nucleus%R0                                       ! Eq. (65)
            this%c_4 = 1.25_dp * (1.5_dp / pi)**(2.0_dp / 3.0_dp) * this%c_1         ! Eq. (66)
            this%k_f = (9.0_dp * pi * Z / (4.0_dp * A))**(1.0_dp / 3.0_dp) / nucleus%R0  ! Eq. (68)
            this%f_kfrp = -0.125_dp * r_p**2 * e2 / nucleus%R0**3 &                    ! Eq. (67)
                          * (145.0_dp / 48.0_dp - 327.0_dp / 2880.0_dp * (this%k_f * r_p)**2 &
                             + 1527.0_dp / 1209600.0_dp * (this%k_f * r_p)**4)
            this%B_s = shape%relative_area()                                    ! Eq. (70)
            this%B_w = 1.0_dp
            this%delta_n_ave = r_mac * this%B_s / N**(1.0_dp / 3.0_dp)               ! Eq. (63)
            this%delta_p_ave = r_mac * this%B_s / Z**(1.0_dp / 3.0_dp)               ! Eq. (64)
            this%delta_pn_ave = h_np / (this%B_s * A**(2.0_dp / 3.0_dp))             ! Eq. (65)

            if (present(g_mod)) then
                call this%calculate_b_1_b_3(nucleus, g_mod)
            else
                call this%calculate_b_1_b_3_surface(nucleus, shape)
            end if

            this%E_mass = M_H * Z + M_n * N
            this%E_v = -a_v * (1.0_dp - k_v * this%I_n**2) * A
            this%E_s = a_s * (1.0_dp - k_s * this%I_n**2) * A**(2.0_dp / 3.0_dp) * this%B_1
            this%E_a0 = a_0 * this%B_w
            this%E_c = this%c_1 * Z**2 / A**(1.0_dp / 3.0_dp) * this%B_3
            this%E_cex = -this%c_4 * Z**(4.0_dp / 3.0_dp) / A**(1.0_dp / 3.0_dp)
            this%E_pff = this%f_kfrp * Z**2 / A
            this%E_ca = -c_a * (N - Z)
            this%E_W = Wigner * abs(this%I_n) * this%B_w
            if (nucleus%N == nucleus%Z .and. mod(nucleus%N, 2) == 1) this%E_W = this%E_W + Wigner / A
            if (mod(nucleus%N, 2) == 1 .and. mod(nucleus%Z, 2) == 1) then
                this%E_pair = this%delta_n_ave + this%delta_p_ave - this%delta_pn_ave
            else if (mod(nucleus%N, 2) == 1) then
                this%E_pair = this%delta_n_ave
            else if (mod(nucleus%Z, 2) == 1) then
                this%E_pair = this%delta_p_ave
            else
                this%E_pair = 0.0_dp
            end if
            this%E_el = -a_el * Z**2.39_dp

            E_frldm = this%E_mass + this%E_v + this%E_s + this%E_a0 + this%E_c + this%E_cex + this%E_pff &
                    + this%E_ca + this%E_W + this%E_pair + this%E_el
            this%E_frldm = E_frldm

            ! spherical reference values, Eq. (74)
            x_0 = nucleus%R0 * A**(1.0_dp / 3.0_dp) / a_Yukawa
            y_0 = nucleus%R0 * A**(1.0_dp / 3.0_dp) / a_den
            B1_0 = 1.0_dp - 3.0_dp / x_0**2 + (1.0_dp + x_0) * (2.0_dp + 3.0_dp / x_0 + 3.0_dp / x_0**2) &
                 * exp(-2.0_dp * x_0)
            B3_0 = 1.0_dp - 5.0_dp / y_0**2 * (1.0_dp - 15.0_dp / (8.0_dp * y_0) + 21.0_dp / (8.0_dp * y_0**3) &
                 - 0.75_dp * (1.0_dp + 4.5_dp / y_0 + 7.0_dp / y_0**2 + 3.5_dp / y_0**3) * exp(-2.0_dp * y_0))

            if (.not. this%verbose) return
            print '(a, 2f10.6, a, 2f10.6)', " B_1, B_3        : ", this%B_1, this%B_3, &
                  "   sphere (Eq. 74): ", B1_0, B3_0
            print '(a, f10.6)', " B_s: ", this%B_s
            print '(a, 4f12.4)', " E_v, E_s, E_c, E_c,exch [MeV]   : ", this%E_v, this%E_s, this%E_c, this%E_cex
            print '(a, 5f10.4)', " a0, f(kf rp), c_a, Wigner, pair  : ", this%E_a0, this%E_pff, this%E_ca, &
                  this%E_W, this%E_pair
            print '(a, f12.4)', " FRLDM macroscopic mass excess [MeV]: ", E_frldm
        end subroutine calculate_frldm_energy

        subroutine calculate_b_1_b_3_surface(this, nucleus, shape, n_quad)
            !----------------------------
            ! B_1 and B_3 of a sharp-surface shape by double surface integrals.
            ! For a kernel f(sigma) with nabla^2 g = f,
            !   int_V int_V f(|r - r'|) d^3r d^3r' = - oint oint (dS . dS') g(|r - r'|)
            ! (nabla_r . nabla_r' g = - nabla^2 g, Gauss theorem twice; constants in g drop
            ! out because oint dS = 0). With s = sigma / a:
            !   B_1: f = (2 - s) e^-s / s,               g = - a^2 e^-s
            !   B_3: f = 1/sigma - (1 + s/2) e^-s/sigma, g = a [s/2 - e^-s/2 + 2 (1 - e^-s)/s]
            ! Quadrature: Gauss-Legendre in u (theta), trapezoid in v (phi); the outer point
            ! runs over one octant only (all shapes here have D2h symmetry). The B_3 integrand
            ! is smooth (g = 3a/2 + O(s^2)); the |sigma| cusp of e^-s makes the B_1 error
            ! ~ n^-3, which is removed by Richardson extrapolation from n and 3n/2 points.
            ! Sphere, A = 208: B_1 error 3e-5, B_3 error < 1e-9 for n_quad = 64.
            !----------------------------
            implicit none
            class(frldm_variables), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            class(surface_shape), intent(in) :: shape
            integer, intent(in), optional :: n_quad
            integer :: n1, n2
            real(dp) :: s1a, s3a, s1b, s3b, A

            n1 = 64
            if (present(n_quad)) n1 = n_quad
            n1 = 4 * ((n1 + 3) / 4)
            n2 = 4 * ((3 * n1 / 2 + 3) / 4)
            call surface_sums(shape, n1, s1a, s3a)
            call surface_sums(shape, n2, s1b, s3b)
            s1a = (real(n2, dp)**3 * s1b - real(n1, dp)**3 * s1a) / (real(n2, dp)**3 - real(n1, dp)**3)

            A = real(nucleus%A, dp)
            this%B_1 = a_Yukawa**2 * s1a / (8.0_dp * pi**2 * nucleus%R0**2 * a_Yukawa**4) * A**(-2.0_dp / 3.0_dp)
            this%B_3 = -a_den * s3b * 15.0_dp / (32.0_dp * pi**2 * nucleus%R0**5) * A**(-5.0_dp / 3.0_dp)
        end subroutine calculate_b_1_b_3_surface

        subroutine surface_sums(shape, nu, S1, S3)
            ! S1 = oint oint dS.dS' e^-s (a = a_Yukawa), S3 = oint oint dS.dS' G(s) (a = a_den)
            class(surface_shape), intent(in) :: shape
            integer, intent(in) :: nu
            real(dp), intent(out) :: S1, S3
            real(dp), allocatable :: x(:, :), dS(:, :), tq(:), wq(:)
            real(dp) :: sig, d, r(3), nrm(3), s
            integer :: nphi, i, j, m, p, q

            nphi = 2 * nu
            allocate(tq(nu), wq(nu), x(3, nu * nphi), dS(3, nu * nphi))
            call gauss_legendre(nu, tq, wq)
            m = 0
            do j = 1, nphi
                do i = 1, nu
                    m = m + 1
                    call shape%eval(0.5_dp * pi * (tq(i) + 1.0_dp), 2.0_dp * pi * (j - 0.5_dp) / nphi, r, nrm)
                    x(:, m) = r
                    dS(:, m) = nrm * (0.5_dp * pi * wq(i)) * (2.0_dp * pi / nphi)
                end do
            end do
            S1 = 0.0_dp
            S3 = 0.0_dp
            !$omp parallel do default(none) shared(x, dS, m) private(q, sig, d, s) &
            !$omp reduction(+:S1, S3) schedule(dynamic)
            do p = 1, m
                if (any(x(:, p) < 0.0_dp)) cycle
                do q = 1, m
                    sig = norm2(x(:, p) - x(:, q))
                    d = dot_product(dS(:, p), dS(:, q))
                    S1 = S1 + d * exp(-sig / a_Yukawa)
                    s = sig / a_den
                    if (s > 1.0e-4_dp) then
                        S3 = S3 + d * (0.5_dp * s - 0.5_dp * exp(-s) + 2.0_dp * (1.0_dp - exp(-s)) / s)
                    else
                        S3 = S3 + d * (1.5_dp + s**2 / 12.0_dp)
                    end if
                end do
            end do
            !$omp end parallel do
            S1 = 8.0_dp * S1
            S3 = 8.0_dp * S3
        end subroutine surface_sums

end module frldm_mod