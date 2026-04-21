module frldm_mod
    use iso_fortran_env, only: real64
    use constant_mod, only: e2, pi
    use nucleus_mod, only: nucleus_property
    use grid_mod, only: grid_type
    use three_D_derivative_mod
    use array_conversion_mod
    use CG_method_mod, only: CG_method
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
    real(dp), parameter, public :: a_el = 1.43e-5_dp ! fm, diffuseness of charge distribution

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
        real(dp), allocatable :: B_1pot(:)
        real(dp), allocatable :: B_3pot(:)
        real(dp) :: B_w
        real(dp) :: E_frldm

        contains
            procedure :: calculate_frldm_energy
            procedure :: calculate_b_1_b_3
    end type frldm_variables

    public :: frldm_variables

    contains
        subroutine calculate_frldm_energy(this, nucleus, g_mod, E_frldm)
            class(frldm_variables), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            type(grid_type), intent(in) :: g_mod
            real(dp), intent(out) :: E_frldm
            real(dp) :: E_v, E_s, E_c, E_cec, E_pffc, E_ca, E_W, E_pair, E_be
            allocate(this%B_1pot(g_mod%n_x_points * g_mod%n_y_points * g_mod%n_z_points))
            allocate(this%B_3pot(g_mod%n_x_points * g_mod%n_y_points * g_mod%n_z_points))


            !----------------------------
            ! clculate the FRLDM variables
            !----------------------------
            this%I_n =real(nucleus%N - nucleus%Z, dp) / real(nucleus%A, dp)

            this%c_1 = 0.6_dp * e2 / nucleus%R0

            this%c_4 = 0.8_dp *(1.5_dp / pi)**(2.0_dp / 3.0_dp) * this%c_1

            this%k_f = (9.0_dp * pi * real(nucleus%Z, dp) / 4.0_dp / real(nucleus%A, dp))**(1.0_dp / 3.0_dp) / nucleus%R0

            this%f_kfrp = -0.125_dp * r_p**2 * e2**2 * (145.0_dp/48.0_dp - 327.0_dp*(this%k_f * r_p)**2 / 2880.0_dp + 1527.0d0/1209600.0_dp * (this%k_f * r_p)**4) / nucleus%R0**3

            this%B_s = nucleus%surface_area * real(nucleus%A, dp)**(-2.0_dp / 3.0_dp)/(4.0_dp * pi * nucleus%R0**2)

            this%delta_pn_ave = h_np /(this%B_s * real(nucleus%A, dp)**(2.0_dp / 3.0_dp))

            this%delta_n_ave = r_mac * this%B_s / real(nucleus%N, dp)** (1.0_dp / 3.0_dp)

            this%delta_p_ave = r_mac * this%B_s / real(nucleus%Z, dp)** (1.0_dp / 3.0_dp)

            this%B_w = 1.0_dp
            
            call this%calculate_b_1_b_3(nucleus,g_mod, this%B_1, this%B_3)

            !----------------------------
            ! volume energy
            !----------------------------
            E_v = - a_v * (1.0_dp - k_v * this%I_n**2) * real(nucleus%A, dp)

            !----------------------------
            ! surface energy
            !----------------------------
            E_s = a_s * (1.0_dp - k_s * this%I_n**2) * real(nucleus%A, dp)**(2.0_dp / 3.0_dp) * this%B_1
            
            !----------------------------
            ! Coulomb energy
            !----------------------------
            E_c = this%c_1 * real(nucleus%Z, dp)**2 / real(nucleus%A, dp)**(1.0_dp / 3.0_dp)*this%B_1

            !----------------------------
            ! charge exchange correction energy
            !----------------------------
            E_cec = this%c_4 * real(nucleus%Z, dp)**(4.0_dp / 3.0_dp) / real(nucleus%A, dp)**(1.0_dp / 3.0_dp)

            !----------------------------
            ! proton form-factor correction energy
            !----------------------------
            E_pffc = this%f_kfrp * real(nucleus%Z, dp)**2 / real(nucleus%A, dp)

            !----------------------------
            ! charge asymmetry energy
            !----------------------------
            E_ca = c_a * (real(nucleus%Z, dp) - real(nucleus%N, dp))

            !----------------------------
            ! Wigner energy
            !----------------------------
            If (nucleus%N == nucleus%Z .and. mod(nucleus%N, 2) == 1) then
                E_W = (abs(this%I_n)* this%B_w + 1.0_dp/real(nucleus%A, dp))* Wigner
            else
                E_W = abs(this%I_n)* this%B_w * Wigner
            end if

            !----------------------------
            ! pairing energy
            !----------------------------
            if (mod(nucleus%N, 2) == 0 .and. mod(nucleus%Z, 2) == 0) then
                E_pair = 0
            else if (mod(nucleus%N, 2) == 1 .and. mod(nucleus%Z, 2) == 1) then
                E_pair = this%delta_n_ave + this%delta_p_ave - this%delta_pn_ave
            else if (mod(nucleus%N, 2) == 1 .and. mod(nucleus%Z, 2) == 0) then
                E_pair = this%delta_n_ave
            else
                E_pair = this%delta_p_ave
            end if

            !----------------------------
            ! binding energy of electron
            !----------------------------
            E_be = a_el * real(nucleus%Z, dp)**(2.39_dp)



            E_frldm = E_v + E_s + E_c + E_cec + E_pffc + E_ca + E_W + E_pair + E_be

        end subroutine calculate_frldm_energy



        subroutine calculate_b_1_b_3(this, nucleus,g_mod, B_1, B_3)
            class(frldm_variables), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            type(grid_type), intent(in) :: g_mod
            real(dp), intent(out) :: B_1, B_3
            
            integer :: i
            integer :: n_max
            n_max = g_mod%n_x_points * g_mod%n_y_points * g_mod%n_z_points
            ! call CG_method(g_mod%n_x_points, g_mod%n_y_points, g_mod%n_z_points, g_mod%h_x, a_Yukawa, this%B_1pot)
            ! call CG_method(g_mod%n_x_points, g_mod%n_y_points, g_mod%n_z_points, g_mod%h_x, a_den, this%B_3pot)

            do i = 1, n_max
                B_1 = B_1 + this%B_1pot(i) * g_mod%dV
                B_3 = B_3 + this%B_3pot(i) * g_mod%dV
            end do



            

            !----------------------------
            ! Calculate the integrands for B_1 and B_3
            !----------------------------



        end subroutine calculate_b_1_b_3




end module frldm_mod