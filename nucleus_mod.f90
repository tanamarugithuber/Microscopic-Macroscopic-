module nucleus_mod
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    implicit none
    private

    integer, parameter :: dp = real64



    type :: nucleus_property
        !---------------------------
        ! Nucleus fundamental properties
        !---------------------------
        integer :: A  ! mass number
        integer :: Z  ! proton number
        integer :: N  ! neutron number
        real(dp) :: R0 = 1.16_dp ! fm, radius parameter
        real(dp) :: R_a ! fm, average radius

        !---------------------------
        ! Nucleus shape parameters, Z-axis is the symmetry axis
        !---------------------------
        real(dp) :: ecc ! eccentricity
        real(dp) :: semi1 ! fm, semi-major or semi-minor axis
        real(dp) :: semi2 ! fm, semi-major or semi-minor axis
        real(dp) :: volume ! fm^3, volume of the nucleus
        real(dp) :: surface_area ! fm^2, surface area of the nucleus

        contains
            procedure :: calculate_fundamental_properties

    end type nucleus_property

    public :: nucleus_property

    contains
        subroutine calculate_fundamental_properties(this)
            implicit none
            class(nucleus_property), intent(inout) :: this

            this%A = this%N + this%Z

            this%R_a = this%R0 * this%A**(1.0_dp / 3.0_dp)

            this%volume = (4.0_dp / 3.0_dp) * pi * this%R_a**3

            this%semi2 = this%volume / (4.0_dp / 3.0_dp * pi * this%semi1**2)

            if (this%semi1 > this%semi2) then
                this%ecc = sqrt(1.0_dp - (this%semi2 / this%semi1)**2)
                this%surface_area = 2.0_dp * pi * (this%semi1**2 + this%semi1*this%semi2 * asin(this%ecc) / this%ecc)
            else if (this%semi1 < this%semi2) then
                this%ecc = sqrt(1.0_dp - (this%semi1 / this%semi2)**2)
                this%surface_area = 2.0_dp * pi * (this%semi1**2 + this%semi2**2 * atanh(this%ecc) / this%ecc)
            else
                this%ecc = 0.0_dp
                this%surface_area = 4.0_dp * pi * this%semi1**2
            end if
            
        end subroutine calculate_fundamental_properties


end module nucleus_mod