module frldm_mod
    use iso_fortran_env, only: real64
    use constant_mod, only: e2, pi
    use nucleus_mod, only: nucleus_property
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
        real(dp) :: f_kfr_p 
        real(dp) :: k_f
        real(dp) :: delta_pn_ave ! neutron-proton pairing correction
        real(dp) :: delta_n_ave ! neutron pairing correction
        real(dp) :: delta_p_ave ! proton pairing correction
        real(dp) :: B_s

        contains
            procedure :: calculate_frldm_variables
    end type frldm_variables

    public :: frldm_variables

    contains
        subroutine calculate_frldm_variables(this, nucleus)
            class(frldm_variables), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus

            this%I_n =real(nucleus%N - nucleus%Z, dp) / real(nucleus%A, dp)

            this%c_1 = 0.6_dp * e2 / nucleus%R0

            this%c_4 = 0.8_dp *(1.5_dp / pi)**(2.0_dp / 3.0_dp) * this%c_1

            this%k_f = (9.0_dp * pi * real(nucleus%Z, dp) / 4.0_dp / real(nucleus%A, dp))**(1.0_dp / 3.0_dp) / nucleus%R0

            this%f_kfr_p = -0.125_dp * r_p**2 * e2**2 * (145.0_dp/48.0_dp - 327.0_dp*(this%k_f * r_p)**2 / 2880.0_dp + 1527.0d0/1209600.0_dp * (this%k_f * r_p)**4) / nucleus%R0**3

            this%B_s = nucleus%surface_area * real(nucleus%A, dp)**(-2.0_dp / 3.0_dp)/(4.0_dp * pi * nucleus%R0**2)

            this%delta_pn_ave = h_np /(this%B_s * real(nucleus%A, dp)**(2.0_dp / 3.0_dp))

            this%delta_n_ave = r_mac * this%B_s / real(nucleus%N, dp)** (1.0_dp / 3.0_dp)

            this%delta_p_ave = r_mac * this%B_s / real(nucleus%Z, dp)** (1.0_dp / 3.0_dp)
        end subroutine calculate_frldm_variables




end module frldm_mod