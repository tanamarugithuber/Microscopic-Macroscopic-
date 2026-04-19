module grid_mod
    use iso_fortran_env, only: real64
    use nucleus_mod, only: nucleus_property

    implicit none
    private

    integer, parameter :: dp = real64

    type :: grid_type
        !---------------------------
        ! Grid properties
        !---------------------------
        real(dp) :: h_x = 0.1_dp ! grid spacing in fm
        real(dp) :: h_y = 0.1_dp ! grid spacing in fm
        real(dp) :: h_z = 0.1_dp ! grid spacing in fm
        real(dp) :: dV ! grid volume element in fm^3
        integer :: n_points ! total number of grid points
        integer :: n_x_points
        integer :: n_y_points
        integer :: n_z_points

        !---------------------------
        ! Grid boundaries in fm
        ! number of grid points from the center to the edge of the nucleus.   
        ! Defined by the size of nucleus.
        !---------------------------
        integer :: b_x_min
        integer :: b_x_max
        integer :: b_y_min
        integer :: b_y_max
        integer :: b_z_min
        integer :: b_z_max
    

        !---------------------------
        ! number of grid points from center to the box in fm
        !---------------------------
        integer :: n_x_min 
        integer :: n_x_max 
        integer :: n_y_min
        integer :: n_y_max
        integer :: n_z_min
        integer :: n_z_max

        integer :: n_times = 2
        ! This represents how many times the grid extends beyond the nucleus size. 
        ! For example, if n_times = 2, the grid extends to twice the size of the 
        ! nucleus in each direction.

        contains
            procedure :: initialize_grid
    end type grid_type


    public :: grid_type

    contains
        subroutine initialize_grid(this, nucleus)
            class(grid_type), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus

            ! Calculate grid boundaries based on the nucleus size
            this%b_x_min = - nint(nucleus%semi1 / this%h_x)
            this%b_x_max = nint(nucleus%semi1 / this%h_x)
            this%b_y_min = - nint(nucleus%semi1 / this%h_y)
            this%b_y_max = nint(nucleus%semi1 / this%h_y)
            this%b_z_min = - nint(nucleus%semi2 / this%h_z)
            this%b_z_max = nint(nucleus%semi2 / this%h_z)

            this%n_x_min = this%b_x_min * this%n_times
            this%n_x_max = this%b_x_max * this%n_times
            this%n_y_min = this%b_y_min * this%n_times
            this%n_y_max = this%b_y_max * this%n_times
            this%n_z_min = this%b_z_min * this%n_times
            this%n_z_max = this%b_z_max * this%n_times

            ! Calculate the volume element of the grid
            this%dV = this%h_x * this%h_y * this%h_z
            this%n_points = (this%n_x_max - this%n_x_min + 1) * (this%n_y_max - this%n_y_min + 1) * (this%n_z_max - this%n_z_min + 1)
            this%n_x_points = this%n_x_max - this%n_x_min + 1
            this%n_y_points = this%n_y_max - this%n_y_min + 1
            this%n_z_points = this%n_z_max - this%n_z_min + 1
        end subroutine initialize_grid


end module grid_mod