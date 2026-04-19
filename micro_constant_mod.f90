module micro_constant_mod
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use nucleus_mod, only: nucleus_property
    implicit none
    private

    integer, parameter :: dp = real64

    !---------------------------
    ! Constants for Microscopic Methods
    ! Units:
    !   Energy : MeV
    !   Length : fm
    !---------------------------

    !---------------------------
    ! second category constants
    !---------------------------
    real(dp), parameter, public :: C_cur = 41.0_dp      ! basis curvature constant in MeV
    real(dp), parameter, public :: V_s = 52.5_dp      ! symmetric potential-depth constant in MeV
    real(dp), parameter, public :: V_a = 48.7_dp      ! asymmetric potential-depth constant in MeV
    real(dp), parameter, public :: A_den = 0.82_dp ! potential radius correction constant in fm
    real(dp), parameter, public :: B_den = 0.56_dp ! potential radius curvature-correction constant in fm**2
    real(dp), parameter, public :: a_pot = 0.8_dp     ! potential diffuseness constant in fm
    real(dp), parameter, public :: k_p = 0.025_dp     ! proton spin–orbit A coefficient 
    real(dp), parameter, public :: l_p = 28.0_dp      ! proton spin–orbit constant
    real(dp), parameter, public :: k_n = 0.01875_dp   ! neutron spin–orbit A coefficient
    real(dp), parameter, public :: l_n = 31.5_dp      !  neutron spin–orbit constant

    !---------------------------
    ! third category constants
    !---------------------------
    real(dp), parameter, public :: N_bas = 12.0_dp ! number of basis functions
    real(dp), parameter, public :: p = 8.0_dp      ! order of Strutinsky shell correction
    real(dp), parameter, public :: C_s = 1.0_dp     ! Strutinsky range coefficient

    !---------------------------
    ! fourth category constants
    !---------------------------
    real(dp), parameter, public :: r_mic = 3.2_dp  ! LN effective-interaction pairing-gap constant
    real(dp), parameter, public :: frK = 0.2475_dp ! Zero-point energy constant
    real(dp), parameter, public :: a_2 = 22.00_dp  ! surface-energy constant in MeV
    real(dp), parameter, public :: J = 35.0_dp  ! symmetry-energy constant in MeV
    real(dp), parameter, public :: L = 99.0_dp  ! density-symmetry constant in MeV
    real(dp), parameter, public :: Q = 25.0_dp  ! effective surface-stiffness constant in MeV
    real(dp), parameter, public :: K = 300.0_dp  ! compressibility constant in MeV

    !---------------------------
    ! Variables for Microscopic calculations
    !---------------------------
    type :: microscopic_variables
        real(dp) :: delta_ave
        real(dp) :: epsilon_ave
        contains
            procedure :: calculate_pairing_gap
    end type microscopic_variables

    public :: microscopic_variables

    contains
        subroutine calculate_pairing_gap(this, nucleus)
            class(microscopic_variables), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus

            this%delta_ave = ((nucleus%N - nucleus%Z) / nucleus%A * 0.0112_dp * nucleus%Z**2/ nucleus%A**(5.0_dp / 3.0_dp))/(1.0_dp + 3.15_dp * nucleus%A**(-1.0_dp / 3.0_dp))
            this%epsilon_ave = -0.147_dp*nucleus%A**(1.0_dp / 3.0_dp) + 0.330_dp * this%delta_ave**2 + 0.00248_dp * nucleus%Z**2 / nucleus%A**(4.0_dp / 3.0_dp)
        end subroutine calculate_pairing_gap



end module micro_constant_mod