module constant_mod
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter :: dp = real64

    
        !---------------------------
        ! Fundamental constants
        ! Units:
        !   Energy : MeV
        !   Length : fm
        !---------------------------

        real(dp), parameter, public :: pi = 3.14159265358979323846_dp
        real(dp), parameter, public :: hbar_c = 197.32891_dp      ! MeV*fm
        real(dp), parameter, public :: hbar = 6.582119569e-22_dp  ! MeV*s
        real(dp), parameter, public :: c = 299792458._dp          ! m/s
        real(dp), parameter, public :: fm_to_m = 1e-15_dp         ! fm -> m
        real(dp), parameter, public :: m_nucleon = 938.90595_dp   ! MeV/c^2
        real(dp), parameter, public :: e2 = 1.4399764_dp          ! MeV*fm
        real(dp), parameter, public :: alpha = 1._dp / 137.035999084_dp
        complex(dp), parameter, public :: i_complex = (0.0_dp, 1.0_dp)


end module constant_mod