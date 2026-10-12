module grid_mod
    !$ use omp_lib
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use nucleus_mod, only: nucleus_property
    use shape_mod, only: epsilon_shape, epsilon_shape_from

    implicit none
    private

    integer, parameter :: dp = real64

    type :: grid_type
        !---------------------------
        ! Grid properties
        !---------------------------
        real(dp) :: h_x  ! grid spacing in fm
        real(dp) :: h_y  ! grid spacing in fm
        real(dp) :: h_z  ! grid spacing in fm
        real(dp) :: dV ! grid volume element in fm^3
        logical :: cubic = .true. ! .true. : cubic box sized by the longest semi-axis (CG_method_mod needs n_x = n_y = n_z)
                                  ! .false.: box extent chosen separately for each axis

        !---------------------------
        ! Nuclear shape: epsilon parameterization (eps2, gamma, eps4, eps6) of shape_mod
        !---------------------------
        type(epsilon_shape) :: shape
        real(dp) :: semi(3) ! fm, largest extent of the shape along x, y, z (semi-axes for an ellipsoid)
        integer :: n_sub = 6 ! cells cut by the surface get the fraction of n_sub**3 sub-points inside
                             ! (1: only the cell centre is tested, density 0 or 1)
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

        integer :: n_times 
        ! This represents how many times the grid extends beyond the nucleus size. 
        ! For example, if n_times = 2, the grid extends to twice the size of the 
        ! nucleus in each direction.

        !---------------------------
        ! Grid points in 3D space can be calculated as (n_x_max - n_x_min + 1) * (n_y_max - n_y_min + 1) * (n_z_max - n_z_min + 1)
        !---------------------------
        !real(dp), allocatable :: r(:) ! r(i) gives the radial coordinate of the grid point i
        real(dp), allocatable :: density_index(:) ! density_index(i,j,k) gives the wheather the grid point (i,j,k) is inside the nucleus (1) or outside the nucleus (0). And it can be used as a density if the density is assumed to be constant inside the nucleus and zero outside the nucleus.
        real(dp), allocatable :: density_index3D(:,:,:) ! density_index3D(i,j,k) gives the wheather the grid point (i,j,k) is inside the nucleus (1) or outside the nucleus (0). And it can be used as a density if the density is assumed to be constant inside the nucleus and zero outside the nucleus.



        contains
            procedure :: initialize_grid
            procedure :: inside_outside_nucleus
            procedure :: inside_outside_nucleus3D
            procedure :: deformation_parameterization
            procedure :: cell_fraction
    end type grid_type


    public :: grid_type
    public :: semi_axes_from_eps_gamma

    contains
        function semi_axes_from_eps_gamma(R, eps, gamma) result(semi)
            !---------------------------
            ! Semi-axes of a triaxial ellipsoid from the Nilsson (eps, gamma) parameters.
            ! They are inversely proportional to the oscillator frequencies
            !   omega_x = omega0 [1 - 2/3 eps cos(gamma + 2pi/3)]
            !   omega_y = omega0 [1 - 2/3 eps cos(gamma - 2pi/3)]
            !   omega_z = omega0 [1 - 2/3 eps cos(gamma)]
            ! and normalised to R_x R_y R_z = R^3 (volume conservation).
            ! Same convention as ho_basis_mod and potential_gh_mod.
            ! gamma = 0: axially symmetric about z (eps > 0 prolate, eps < 0 oblate). gamma in radian.
            !---------------------------
            implicit none
            real(dp), intent(in) :: R, eps, gamma
            real(dp) :: semi(3)
            real(dp) :: f(3)

            f(1) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma + 2.0_dp * pi / 3.0_dp)
            f(2) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma - 2.0_dp * pi / 3.0_dp)
            f(3) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma)
            if (any(f <= 0.0_dp)) stop "semi_axes_from_eps_gamma: eps too large"
            semi = R / f * product(f)**(1.0_dp / 3.0_dp)
        end function semi_axes_from_eps_gamma

        subroutine initialize_grid(this, nucleus)
            !---------------------------
            ! The box is determined from the deformation (nucleus%epsilon2, gamma_deg, epsilon4, epsilon6):
            !   largest extents R_x, R_y, R_z of the shape  ->  nu_i = ceiling(R_i / h_i)
            !   box  = n_times * nucleus size
            ! cubic = .true. uses the longest extent on all three axes.
            !---------------------------
            implicit none
            class(grid_type), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            integer :: nu(3)

            print *, "Initializing grid based on the nucleus properties..."
            this%shape = epsilon_shape_from(nucleus%R_a, nucleus%epsilon2, nucleus%gamma_deg * pi / 180.0_dp, &
                                            nucleus%epsilon4, nucleus%epsilon6)
            this%semi = this%shape%extent()

            ! Calculate grid boundaries based on the nucleus size
            if (this%cubic) then
                nu = ceiling(maxval(this%semi) / this%h_x)
            else
                nu = ceiling(this%semi / [this%h_x, this%h_y, this%h_z])
            end if
            this%nu_x_min = - nu(1)
            this%nu_x_max = nu(1)
            this%nu_y_min = - nu(2)
            this%nu_y_max = nu(2)
            this%nu_z_min = - nu(3)
            this%nu_z_max = nu(3)

            this%n_x_min = this%nu_x_min * this%n_times
            this%n_x_max = this%nu_x_max * this%n_times
            this%n_y_min = this%nu_y_min * this%n_times
            this%n_y_max = this%nu_y_max * this%n_times
            this%n_z_min = this%nu_z_min * this%n_times
            this%n_z_max = this%nu_z_max * this%n_times

            ! Calculate the volume element of the grid
            this%dV = this%h_x * this%h_y * this%h_z

            ! Calculate the total number of grid points and allocate the grid arrays
            this%n_points = (this%n_x_max - this%n_x_min + 1) * (this%n_y_max - this%n_y_min + 1) & 
                            * (this%n_z_max - this%n_z_min + 1)
            this%n_x_points = this%n_x_max - this%n_x_min + 1
            this%n_y_points = this%n_y_max - this%n_y_min + 1
            this%n_z_points = this%n_z_max - this%n_z_min + 1

            print *, "eps2, gamma [deg], eps4, eps6: ", nucleus%epsilon2, nucleus%gamma_deg, nucleus%epsilon4, &
                     nucleus%epsilon6
            print *, "extent (x, y, z) [fm]: ", this%semi
            print *, "Grid boundaries calculated:"
            print *, "n_x_min: ", this%n_x_min, " n_x_max: ", this%n_x_max
            print *, "n_y_min: ", this%n_y_min, " n_y_max: ", this%n_y_max
            print *, "n_z_min: ", this%n_z_min, " n_z_max: ", this%n_z_max
            print *, "Total number of grid points calculated for x direction: ", this%n_x_points
            print *, "Total number of grid points calculated for y direction: ", this%n_y_points
            print *, "Total number of grid points calculated for z direction: ", this%n_z_points
            print *, "Total number of grid points: ", this%n_points

            print *, "Grid boundaries and size calculated:"
        end subroutine initialize_grid

        subroutine inside_outside_nucleus(this, nucleus, index)
            implicit none
            class(grid_type), intent(in) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp), intent(out) :: index(this%n_points)
            ! real(dp), intent(out) :: r(this%n_points)
            integer :: i, j, k, l
            real(dp) :: x, y, z
            print *, "Determining which grid points are inside or outside the nucleus..."
        
            do k = 1, this%n_z_points
                z = (this%n_z_min + k - 1) * this%h_z
                do j = 1, this%n_y_points
                    y = (this%n_y_min + j - 1) * this%h_y
                    do i = 1, this%n_x_points
                        x = (this%n_x_min + i - 1) * this%h_x
                        l = (k-1)*this%n_x_points*this%n_y_points + (j-1)*this%n_x_points + i
                        ! r(l) = sqrt(x**2 + y**2 + z**2)
                        index(l) = this%cell_fraction(x, y, z)
                    end do
                end do
            end do
            print *, "Number of grid points inside the nucleus: ", count(index == 1.0_dp), &
                     "  cut by the surface: ", count(index > 0.0_dp .and. index < 1.0_dp)
            print *, "Number of grid points:", this%n_points
            print *, "Grid points classified as inside or outside the nucleus."

        end subroutine inside_outside_nucleus

        subroutine inside_outside_nucleus3D(this, nucleus, index)
            implicit none
            class(grid_type), intent(in) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp), intent(out) :: index(this%n_x_points, this%n_y_points, this%n_z_points)
            ! real(dp), intent(out) :: r(this%n_points)
            integer :: i, j, k
            real(dp) :: x, y, z
            print *, "Determining which grid points are inside or outside the nucleus..."
        
            do k = 1, this%n_z_points
                z = (this%n_z_min + k - 1) * this%h_z
                do j = 1, this%n_y_points
                    y = (this%n_y_min + j - 1) * this%h_y
                    do i = 1, this%n_x_points
                        x = (this%n_x_min + i - 1) * this%h_x
                        ! l = (k-1)*this%n_x_points*this%n_y_points + (j-1)*this%n_x_points + i
                        ! r(l) = sqrt(x**2 + y**2 + z**2)
                        index(i,j,k) = this%cell_fraction(x, y, z)
                    end do
                end do
            end do
            print *, "Number of grid points inside the nucleus: ", count(index == 1.0_dp), &
                     "  cut by the surface: ", count(index > 0.0_dp .and. index < 1.0_dp)
            print *, "Number of grid points:", this%n_points
            print *, "Grid points classified as inside or outside the nucleus."

        end subroutine inside_outside_nucleus3D

        pure function cell_fraction(this, x, y, z) result(frac)
            !---------------------------
            ! fraction of the cell centred at (x, y, z) inside the nucleus:
            ! 1 or 0 if the centre and the 8 corners agree, otherwise the fraction of
            ! n_sub**3 sub-points inside. sum(index) * dV then reproduces the volume closely.
            !---------------------------
            implicit none
            class(grid_type), intent(in) :: this
            real(dp), intent(in) :: x, y, z
            real(dp) :: frac, h(3), p(3)
            integer :: a, b, c, n_in

            h = [this%h_x, this%h_y, this%h_z]
            n_in = merge(1, 0, this%shape%inside([x, y, z]))
            do c = -1, 1, 2
                do b = -1, 1, 2
                    do a = -1, 1, 2
                        n_in = n_in + merge(1, 0, this%shape%inside([x, y, z] + 0.5_dp * h * [a, b, c]))
                    end do
                end do
            end do
            if (this%n_sub <= 1 .or. n_in == 0 .or. n_in == 9) then
                frac = merge(1.0_dp, 0.0_dp, this%shape%inside([x, y, z]))
                return
            end if
            n_in = 0
            do c = 1, this%n_sub
                do b = 1, this%n_sub
                    do a = 1, this%n_sub
                        p = [x, y, z] + h * ([a, b, c] - 0.5_dp * (this%n_sub + 1)) / this%n_sub
                        if (this%shape%inside(p)) n_in = n_in + 1
                    end do
                end do
            end do
            frac = real(n_in, dp) / real(this%n_sub, dp)**3
        end function cell_fraction

        subroutine deformation_parameterization(this, nucleus, beta2, beta4)
            implicit none
            class(grid_type), intent(in) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp), intent(in) :: beta2, beta4
            print *, "Calculating deformation parameterization for the grid..."
            ! This subroutine can be used to calculate the deformation of the nucleus and adjust the grid accordingly.

            
        end subroutine deformation_parameterization



end module grid_mod