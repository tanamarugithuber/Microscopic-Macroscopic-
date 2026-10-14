module micro_constant_mod
    !$use omp_lib
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
    real(dp), parameter, public :: r_mic = 3.2_dp  ! MeV, LN effective-interaction pairing-gap constant
    real(dp), parameter, public :: mi_r_mac = 4.80_dp ! MeV, average pairing-gap constant of the macroscopic
                                                      ! model (= r_mac of frldm_mod), Eqs. (63), (64)
    real(dp), parameter, public :: frK = 0.2475_dp ! Zero-point energy constant
    real(dp), parameter, public :: a_2 = 22.00_dp  ! surface-energy constant in MeV
    real(dp), parameter, public :: mi_J = 35.0_dp  ! symmetry-energy constant in MeV
    real(dp), parameter, public :: mi_L = 99.0_dp  ! density-symmetry constant in MeV
    real(dp), parameter, public :: mi_Q = 25.0_dp  ! effective surface-stiffness constant in MeV
    real(dp), parameter, public :: mi_K = 300.0_dp  ! compressibility constant in MeV

    !---------------------------
    ! Quantities of the FRDM2012 single-particle potential, Eqs. (82)-(93)
    !---------------------------
    ! 【日本語】一粒子ポテンシャル (FRDM2012 式 80-93) を決める量
    !   delta_ave   : 平均の体積非対称度 delta-bar (液滴模型, 式 86)
    !   epsilon_ave : 核内部の密度の平均相対ずれ epsilon-bar (式 87)
    !   R_den       : 密度分布の半径 r0 A^(1/3) (1 + epsilon-bar) [fm]           (式 83)
    !   R_pot       : ポテンシャルの半径 R_den + A_den - B_den / R_den [fm]     (式 82)
    !   e_rho_c     : e x 電荷密度 = Z e^2 / (4/3 pi R_pot^3) [MeV/fm^2]        (式 92)
    !   hbar_omega0 : 振動子基底の hbar omega0 = C_cur / A^(1/3) [MeV]            (式 93)
    !   lambda_p, lambda_n : スピン軌道の強さ k A + l (無次元)                    (式 89, 90)
    !   V_p, V_n    : 陽子・中性子ポテンシャルの深さ V_s +- V_a delta-bar [MeV] (式 84, 85)
    type :: microscopic_variables
        real(dp) :: delta_ave      ! average bulk asymmetry, Eq. (86)
        real(dp) :: epsilon_ave    ! average relative deviation of the bulk density, Eq. (87)
        real(dp) :: R_den          ! fm, Eq. (83)
        real(dp) :: R_pot          ! fm, Eq. (82)
        real(dp) :: e_rho_c        ! MeV/fm^2, e times the charge density, Eq. (92)
        real(dp) :: hbar_omega0    ! MeV, C_cur / A^(1/3)
        real(dp) :: lambda_p       ! spin-orbit strengths, Eq. (89)
        real(dp) :: lambda_n
        real(dp) :: V_p            ! MeV, potential depths, Eqs. (84), (85)
        real(dp) :: V_n
        contains
            procedure :: calculate_microscopic_variables
    end type microscopic_variables

    public :: microscopic_variables

    contains
        subroutine calculate_microscopic_variables(this, nucleus)
            !---------------------------
            ! Droplet-model delta_ave and epsilon_ave (Eqs. 86, 87; numerically Eqs. 119, 120)
            ! and the radius, depths and spin-orbit strengths of the single-particle potential.
            !---------------------------
            ! 【日本語】
            ! 【目的】Z, N から一粒子ポテンシャルのパラメータ (上の型の各成分) を計算する.
            ! 【式】c1 = 3 e^2 / (5 r0) として
            !   delta-bar   = [I + 3/8 c1/Q Z^2/A^(5/3)] / [1 + 9/4 J/Q A^(-1/3)]
            !   epsilon-bar = [-2 a_2 A^(-1/3) + L delta-bar^2 + c1 Z^2/A^(4/3)] / K
            !   (J, L, Q, K, a_2 は液滴模型の定数 mi_J, mi_L, mi_Q, mi_K, a_2)
            implicit none
            class(microscopic_variables), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp) :: Z, N, A, I, c1

            Z = real(nucleus%Z, dp)
            N = real(nucleus%N, dp)
            A = real(nucleus%A, dp)
            I = (N - Z) / A
            c1 = 3.0_dp * e2 / (5.0_dp * nucleus%R0)
            this%delta_ave = (I + 3.0_dp / 8.0_dp * c1 / mi_Q * Z**2 / A**(5.0_dp / 3.0_dp)) &
                           / (1.0_dp + 9.0_dp / 4.0_dp * mi_J / mi_Q / A**(1.0_dp / 3.0_dp))
            this%epsilon_ave = (-2.0_dp * a_2 / A**(1.0_dp / 3.0_dp) + mi_L * this%delta_ave**2 &
                              + c1 * Z**2 / A**(4.0_dp / 3.0_dp)) / mi_K
            this%R_den = nucleus%R0 * A**(1.0_dp / 3.0_dp) * (1.0_dp + this%epsilon_ave)
            this%R_pot = this%R_den + mi_A_den - mi_B_den / this%R_den
            this%e_rho_c = Z * e2 / (4.0_dp * pi * this%R_pot**3 / 3.0_dp)
            this%hbar_omega0 = C_cur / A**(1.0_dp / 3.0_dp)
            this%lambda_p = k_p * A + l_p
            this%lambda_n = k_n * A + l_n
            this%V_p = V_s + V_a * this%delta_ave
            this%V_n = V_s - V_a * this%delta_ave
        end subroutine calculate_microscopic_variables

end module micro_constant_mod
