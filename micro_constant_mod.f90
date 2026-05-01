module micro_constant_mod
    use iso_fortran_env, only: real64
    use constant_mod, only: pi, e2
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
    real(dp), parameter, public :: mi_A_den = 0.82_dp ! potential radius correction constant in fm
    real(dp), parameter, public :: mi_B_den = 0.56_dp ! potential radius curvature-correction constant in fm**2
    real(dp), parameter, public :: a_pot = 0.8_dp     ! potential diffuseness constant in fm
    real(dp), parameter, public :: k_p = 0.025_dp     ! proton spin–orbit A coefficient 
    real(dp), parameter, public :: l_p = 28.0_dp      ! proton spin–orbit constant
    real(dp), parameter, public :: k_n = 0.01875_dp   ! neutron spin–orbit A coefficient
    real(dp), parameter, public :: l_n = 31.5_dp      !  neutron spin–orbit constant

    !---------------------------
    ! third category constants
    !---------------------------
    real(dp), parameter, public :: N_bas = 12.0_dp ! number of basis functions
    real(dp), parameter, public :: mi_p = 8.0_dp      ! order of Strutinsky shell correction
    real(dp), parameter, public :: C_s = 1.0_dp     ! Strutinsky range coefficient

    !---------------------------
    ! fourth category constants
    !---------------------------
    real(dp), parameter, public :: r_mic = 3.2_dp  ! LN effective-interaction pairing-gap constant
    real(dp), parameter, public :: frK = 0.2475_dp ! Zero-point energy constant
    real(dp), parameter, public :: a_2 = 22.00_dp  ! surface-energy constant in MeV
    real(dp), parameter, public :: mi_J = 35.0_dp  ! symmetry-energy constant in MeV
    real(dp), parameter, public :: mi_L = 99.0_dp  ! density-symmetry constant in MeV
    real(dp), parameter, public :: mi_Q = 25.0_dp  ! effective surface-stiffness constant in MeV
    real(dp), parameter, public :: mi_K = 300.0_dp  ! compressibility constant in MeV

    !---------------------------
    ! Variables for Microscopic calculations
    !---------------------------
    type :: microscopic_variables
        real(dp) :: delta_ave
        real(dp) :: epsilon_ave
        real(dp) :: R_den
        real(dp) :: R_pot
        real(dp) :: e_rho_c
        real(dp) :: hbar_omega0
        real(dp), allocatable :: mi_density_index(:) ! density_index(i,j,k) gives the wheather the grid point (i,j,k) is inside the nucleus (1) or outside the nucleus (0). And it can be used as a density if the density is assumed to be constant inside the nucleus and zero outside the nucleus.
        contains
            procedure :: calculate_microscopic_variables
    end type microscopic_variables

    public :: microscopic_variables

    contains
        subroutine calculate_microscopic_variables(this, nucleus)
            class(microscopic_variables), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus

            this%delta_ave = ((nucleus%N - nucleus%Z) / nucleus%A * 0.0112_dp * nucleus%Z**2&
                                / nucleus%A**(5.0_dp / 3.0_dp))/(1.0_dp + 3.15_dp * nucleus%A**(-1.0_dp / 3.0_dp))
            this%epsilon_ave = -0.147_dp*nucleus%A**(1.0_dp / 3.0_dp) + 0.330_dp * this%delta_ave**2 + 0.00248_dp &
                                * nucleus%Z**2 / nucleus%A**(4.0_dp / 3.0_dp)
            this%R_den = nucleus%R0 * nucleus%A**(1.0_dp / 3.0_dp) * (1.0_dp + this%epsilon_ave)
            this%R_pot = this%R_den + mi_A_den - mi_B_den / this%R_den
            this%e_rho_c = nucleus%Z * e2 * 3.0_dp / (4.0_dp * pi * this%R_pot**3)
            this%hbar_omega0 = C_cur / nucleus%A**(1.0_dp / 3.0_dp) 

        end subroutine calculate_microscopic_variables

end module micro_constant_mod