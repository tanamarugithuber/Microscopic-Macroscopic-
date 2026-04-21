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
        real(dp) :: h_x = 0.4_dp ! grid spacing in fm
        real(dp) :: h_y = 0.4_dp ! grid spacing in fm
        real(dp) :: h_z = 0.4_dp ! grid spacing in fm
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
        integer :: nu_x_min
        integer :: nu_x_max
        integer :: nu_y_min
        integer :: nu_y_max
        integer :: nu_z_min
        integer :: nu_z_max
    

        !---------------------------
        ! number of grid points from center to the box in fm
        !---------------------------
        integer :: n_x_min 
        integer :: n_x_max 
        integer :: n_y_min
        integer :: n_y_max
        integer :: n_z_min
        integer :: n_z_max

        integer :: n_times = 3
        ! This represents how many times the grid extends beyond the nucleus size. 
        ! For example, if n_times = 2, the grid extends to twice the size of the 
        ! nucleus in each direction.

        !---------------------------
        ! Grid points in 3D space can be calculated as (n_x_max - n_x_min + 1) * (n_y_max - n_y_min + 1) * (n_z_max - n_z_min + 1)
        !---------------------------
        real(dp), allocatable :: x(:,:,:)
        real(dp), allocatable :: density_index(:,:,:) ! density_index(i,j,k) gives the wheather the grid point (i,j,k) is inside the nucleus (1) or outside the nucleus (0). And it can be used as a density if the density is assumed to be constant inside the nucleus and zero outside the nucleus.



        contains
            procedure :: initialize_grid
            procedure :: inside_outside_nucleus
    end type grid_type


    public :: grid_type

    contains
        subroutine initialize_grid(this, nucleus)
            class(grid_type), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus

            ! Calculate grid boundaries based on the nucleus size
            if (nucleus%semi1 >= nucleus%semi2) then
                this%nu_x_min = - nint(nucleus%semi1 / this%h_x)
                this%nu_x_max = nint(nucleus%semi1 / this%h_x)
            else
                this%nu_x_min = - nint(nucleus%semi2 / this%h_x)
                this%nu_x_max = nint(nucleus%semi2 / this%h_x)
            end if
            this%nu_y_min = this%nu_x_min
            this%nu_y_max = this%nu_x_max
            this%nu_z_min = this%nu_x_min
            this%nu_z_max = this%nu_x_max

            this%n_x_min = this%nu_x_min * this%n_times
            this%n_x_max = this%nu_x_max * this%n_times
            this%n_y_min = this%nu_y_min * this%n_times
            this%n_y_max = this%nu_y_max * this%n_times
            this%n_z_min = this%nu_z_min * this%n_times
            this%n_z_max = this%nu_z_max * this%n_times

            ! Calculate the volume element of the grid
            this%dV = this%h_x * this%h_y * this%h_z

            ! Calculate the total number of grid points and allocate the grid arrays
            this%n_points = (this%n_x_max - this%n_x_min + 1) * (this%n_y_max - this%n_y_min + 1) * (this%n_z_max - this%n_z_min + 1)
            this%n_x_points = this%n_x_max - this%n_x_min + 1
            this%n_y_points = this%n_y_max - this%n_y_min + 1
            this%n_z_points = this%n_z_max - this%n_z_min + 1
        end subroutine initialize_grid

        subroutine inside_outside_nucleus(this, nucleus, index)
            
            class(grid_type), intent(in) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp), intent(out) :: index(this%n_x_points, this%n_y_points, this%n_z_points)
            integer :: i, j, k
            real(dp) :: x, y, z
            do k = 1, this%n_z_points
                z = (this%n_z_min + k - 1) * this%h_z
                do j = 1, this%n_y_points
                    y = (this%n_y_min + j - 1) * this%h_y
                    do i = 1, this%n_x_points
                        x = (this%n_x_min + i - 1) * this%h_x
                        if ((x/nucleus%semi1)**2 + (y/nucleus%semi2)**2 + (z/nucleus%semi2)**2 <= 1.0_dp) then
                            index(i,j,k) = 1.0_dp ! inside the nucleus
                        else
                            index(i,j,k) = 0.0_dp ! outside the nucleus
                        end if
                    end do
                end do
            end do

        end subroutine inside_outside_nucleus


end module grid_mod