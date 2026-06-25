module three_D_derivative_cp_mod
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter :: dp = real64

    public :: laplacian_cp
    public :: d_x_cp, d_y_cp, d_z_cp
    public :: dd_x_cp, dd_y_cp, dd_z_cp

contains



    subroutine dd_x_cp(func, hx, second_derivative)
        implicit none
        complex(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hx
        complex(dp), intent(out) :: second_derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (nx < 7) stop "dd_x: size(func,1) must be >= 7"

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    if (i >= 5 .and. i <= nx-4) then

                        second_derivative(i,j,k) = ( &
                        -9.0_dp    * func(i-4,j,k) &
                        +128.0_dp  * func(i-3,j,k) &
                        -1008.0_dp * func(i-2,j,k) &
                        +8064.0_dp * func(i-1,j,k) &
                        -14350.0_dp* func(i  ,j,k) &
                        +8064.0_dp * func(i+1,j,k) &
                        -1008.0_dp * func(i+2,j,k) &
                        +128.0_dp  * func(i+3,j,k) &
                        -9.0_dp    * func(i+4,j,k) ) / (5040.0_dp*hx*hx)

                    else if (i == 1) then

                        second_derivative(i,j,k) = ( &
                        29531.0_dp*func(1,j,k) -138528.0_dp*func(2,j,k) &
                        +312984.0_dp*func(3,j,k) -448672.0_dp*func(4,j,k) &
                        +435330.0_dp*func(5,j,k) -284256.0_dp*func(6,j,k) &
                        +120008.0_dp*func(7,j,k) -29664.0_dp*func(8,j,k) &
                        +3267.0_dp*func(9,j,k) ) / (5040.0_dp*hx*hx)

                    else if (i == 2) then

                        second_derivative(i,j,k) = ( &
                        3267.0_dp*func(1,j,k) +128.0_dp*func(2,j,k) &
                        -20916.0_dp*func(3,j,k) +38556.0_dp*func(4,j,k) &
                        -37030.0_dp*func(5,j,k) +23688.0_dp*func(6,j,k) &
                        -9828.0_dp*func(7,j,k) +2396.0_dp*func(8,j,k) &
                        -261.0_dp*func(9,j,k) ) / (5040.0_dp*hx*hx)

                    else if (i == 3) then

                        second_derivative(i,j,k) = ( &
                        -261.0_dp*func(1,j,k) +5616.0_dp*func(2,j,k) &
                        -9268.0_dp*func(3,j,k) +1008.0_dp*func(4,j,k) &
                        +5670.0_dp*func(5,j,k) -4144.0_dp*func(6,j,k) &
                        +1764.0_dp*func(7,j,k) -432.0_dp*func(8,j,k) &
                        +47.0_dp*func(9,j,k) ) / (5040.0_dp*hx*hx)

                    else if (i == 4) then

                        second_derivative(i,j,k) = ( &
                        47.0_dp*func(1,j,k) -684.0_dp*func(2,j,k) &
                        +7308.0_dp*func(3,j,k) -13216.0_dp*func(4,j,k) &
                        +6930.0_dp*func(5,j,k) -252.0_dp*func(6,j,k) &
                        -196.0_dp*func(7,j,k) +72.0_dp*func(8,j,k) &
                        -9.0_dp*func(9,j,k) ) / (5040.0_dp*hx*hx)

                    else if (i == nx-3) then

                        second_derivative(i,j,k) = ( &
                        47.0_dp*func(nx,j,k) -684.0_dp*func(nx-1,j,k) &
                        +7308.0_dp*func(nx-2,j,k) -13216.0_dp*func(nx-3,j,k) &
                        +6930.0_dp*func(nx-4,j,k) -252.0_dp*func(nx-5,j,k) &
                        -196.0_dp*func(nx-6,j,k) +72.0_dp*func(nx-7,j,k) &
                        -9.0_dp*func(nx-8,j,k) ) / (5040.0_dp*hx*hx)

                    else if (i == nx-2) then

                        second_derivative(i,j,k) = ( &
                        -261.0_dp*func(nx,j,k) +5616.0_dp*func(nx-1,j,k) &
                        -9268.0_dp*func(nx-2,j,k) +1008.0_dp*func(nx-3,j,k) &
                        +5670.0_dp*func(nx-4,j,k) -4144.0_dp*func(nx-5,j,k) &
                        +1764.0_dp*func(nx-6,j,k) -432.0_dp*func(nx-7,j,k) &
                        +47.0_dp*func(nx-8,j,k) ) / (5040.0_dp*hx*hx)

                    else if (i == nx-1) then

                        second_derivative(i,j,k) = ( &
                        3267.0_dp*func(nx,j,k) +128.0_dp*func(nx-1,j,k) &
                        -20916.0_dp*func(nx-2,j,k) +38556.0_dp*func(nx-3,j,k) &
                        -37030.0_dp*func(nx-4,j,k) +23688.0_dp*func(nx-5,j,k) &
                        -9828.0_dp*func(nx-6,j,k) +2396.0_dp*func(nx-7,j,k) &
                        -261.0_dp*func(nx-8,j,k) ) / (5040.0_dp*hx*hx)

                    else if (i == nx) then

                        second_derivative(i,j,k) = ( &
                        29531.0_dp*func(nx,j,k) -138528.0_dp*func(nx-1,j,k) &
                        +312984.0_dp*func(nx-2,j,k) -448672.0_dp*func(nx-3,j,k) &
                        +435330.0_dp*func(nx-4,j,k) -284256.0_dp*func(nx-5,j,k) &
                        +120008.0_dp*func(nx-6,j,k) -29664.0_dp*func(nx-7,j,k) &
                        +3267.0_dp*func(nx-8,j,k) ) / (5040.0_dp*hx*hx)

                    end if

                end do
            end do
        end do
    end subroutine dd_x_cp

    subroutine dd_y_cp(func, hy, second_derivative)
        implicit none
        complex(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hy
        complex(dp), intent(out) :: second_derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (ny < 7) stop "dd_y: size(func,2) must be >= 7"

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    if (j >= 5 .and. j <= ny-4) then

                        second_derivative(i,j,k) = ( &
                        -9.0_dp    * func(i,j-4,k) &
                        +128.0_dp  * func(i,j-3,k) &
                        -1008.0_dp * func(i,j-2,k) &
                        +8064.0_dp * func(i,j-1,k) &
                        -14350.0_dp* func(i  ,j,k) &
                        +8064.0_dp * func(i,j+1,k) &
                        -1008.0_dp * func(i,j+2,k) &
                        +128.0_dp  * func(i,j+3,k) &
                        -9.0_dp    * func(i,j+4,k) ) / (5040.0_dp*hy*hy)

                    else if (j == 1) then

                        second_derivative(i,j,k) = ( &
                        29531.0_dp*func(i,1,k) -138528.0_dp*func(i,2,k) &
                        +312984.0_dp*func(i,3,k) -448672.0_dp*func(i,4,k) &
                        +435330.0_dp*func(i,5,k) -284256.0_dp*func(i,6,k) &
                        +120008.0_dp*func(i,7,k) -29664.0_dp*func(i,8,k) &
                        +3267.0_dp*func(i,9,k) ) / (5040.0_dp*hy*hy)

                    else if (j == 2) then

                        second_derivative(i,j,k) = ( &
                        3267.0_dp*func(i,1,k) +128.0_dp*func(i,2,k) &
                        -20916.0_dp*func(i,3,k) +38556.0_dp*func(i,4,k) &
                        -37030.0_dp*func(i,5,k) +23688.0_dp*func(i,6,k) &
                        -9828.0_dp*func(i,7,k) +2396.0_dp*func(i,8,k) &
                        -261.0_dp*func(i,9,k) ) / (5040.0_dp*hy*hy)

                    else if (j == 3) then

                        second_derivative(i,j,k) = ( &
                        -261.0_dp*func(i,1,k) +5616.0_dp*func(i,2,k) &
                        -9268.0_dp*func(i,3,k) +1008.0_dp*func(i,4,k) &
                        +5670.0_dp*func(i,5,k) -4144.0_dp*func(i,6,k) &
                        +1764.0_dp*func(i,7,k) -432.0_dp*func(i,8,k) &
                        +47.0_dp*func(i,9,k) ) / (5040.0_dp*hy*hy)

                    else if (j == 4) then

                        second_derivative(i,j,k) = ( &
                        47.0_dp*func(i,1,k) -684.0_dp*func(i,2,k) &
                        +7308.0_dp*func(i,3,k) -13216.0_dp*func(i,4,k) &
                        +6930.0_dp*func(i,5,k) -252.0_dp*func(i,6,k) &
                        -196.0_dp*func(i,7,k) +72.0_dp*func(i,8,k) &
                        -9.0_dp*func(i,9,k) ) / (5040.0_dp*hy*hy)

                    else if (j == ny-3) then

                        second_derivative(i,j,k) = ( &
                        47.0_dp*func(i,ny,k) -684.0_dp*func(i,ny-1,k) &
                        +7308.0_dp*func(i,ny-2,k) -13216.0_dp*func(i,ny-3,k) &
                        +6930.0_dp*func(i,ny-4,k) -252.0_dp*func(i,ny-5,k) &
                        -196.0_dp*func(i,ny-6,k) +72.0_dp*func(i,ny-7,k) &
                        -9.0_dp*func(i,ny-8,k) ) / (5040.0_dp*hy*hy)

                    else if (j == ny-2) then

                        second_derivative(i,j,k) = ( &
                        -261.0_dp*func(i,ny,k) +5616.0_dp*func(i,ny-1,k) &
                        -9268.0_dp*func(i,ny-2,k) +1008.0_dp*func(i,ny-3,k) &
                        +5670.0_dp*func(i,ny-4,k) -4144.0_dp*func(i,ny-5,k) &
                        +1764.0_dp*func(i,ny-6,k) -432.0_dp*func(i,ny-7,k) &
                        +47.0_dp*func(i,ny-8,k) ) / (5040.0_dp*hy*hy)

                    else if (j == ny-1) then

                        second_derivative(i,j,k) = ( &
                        3267.0_dp*func(i,ny,k) +128.0_dp*func(i,ny-1,k) &
                        -20916.0_dp*func(i,ny-2,k) +38556.0_dp*func(i,ny-3,k) &
                        -37030.0_dp*func(i,ny-4,k) +23688.0_dp*func(i,ny-5,k) &
                        -9828.0_dp*func(i,ny-6,k) +2396.0_dp*func(i,ny-7,k) &
                        -261.0_dp*func(i,ny-8,k) ) / (5040.0_dp*hy*hy)

                    else if (j == ny) then

                        second_derivative(i,j,k) = ( &
                        29531.0_dp*func(i,ny,k) -138528.0_dp*func(i,ny-1,k) &
                        +312984.0_dp*func(i,ny-2,k) -448672.0_dp*func(i,ny-3,k) &
                        +435330.0_dp*func(i,ny-4,k) -284256.0_dp*func(i,ny-5,k) &
                        +120008.0_dp*func(i,ny-6,k) -29664.0_dp*func(i,ny-7,k) &
                        +3267.0_dp*func(i,ny-8,k) ) / (5040.0_dp*hy*hy)

                    end if

                end do
            end do
        end do
    end subroutine dd_y_cp

    subroutine dd_z_cp(func, hz, second_derivative)
        implicit none
        complex(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hz
        complex(dp), intent(out) :: second_derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (nz < 7) stop "dd_z: size(func,3) must be >= 7"

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    if (k >= 5 .and. k <= nz-4) then

                        second_derivative(i,j,k) = ( &
                        -9.0_dp    * func(i,j,k-4) &
                        +128.0_dp  * func(i,j,k-3) &
                        -1008.0_dp * func(i,j,k-2) &
                        +8064.0_dp * func(i,j,k-1) &
                        -14350.0_dp* func(i,j,k) &
                        +8064.0_dp * func(i,j,k+1) &
                        -1008.0_dp * func(i,j,k+2) &
                        +128.0_dp  * func(i,j,k+3) &
                        -9.0_dp    * func(i,j,k+4) ) / (5040.0_dp*hz*hz)

                    else if (k == 1) then

                        second_derivative(i,j,k) = ( &
                        29531.0_dp*func(i,j,1) -138528.0_dp*func(i,j,2) &
                        +312984.0_dp*func(i,j,3) -448672.0_dp*func(i,j,4) &
                        +435330.0_dp*func(i,j,5) -284256.0_dp*func(i,j,6) &
                        +120008.0_dp*func(i,j,7) -29664.0_dp*func(i,j,8) &
                        +3267.0_dp*func(i,j,9) ) / (5040.0_dp*hz*hz)

                    else if (k == 2) then

                        second_derivative(i,j,k) = ( &
                        3267.0_dp*func(i,j,1) +128.0_dp*func(i,j,2) &
                        -20916.0_dp*func(i,j,3) +38556.0_dp*func(i,j,4) &
                        -37030.0_dp*func(i,j,5) +23688.0_dp*func(i,j,6) &
                        -9828.0_dp*func(i,j,7) +2396.0_dp*func(i,j,8) &
                        -261.0_dp*func(i,j,9) ) / (5040.0_dp*hz*hz)

                    else if (k == 3) then

                        second_derivative(i,j,k) = ( &
                        -261.0_dp*func(i,j,1) +5616.0_dp*func(i,j,2) &
                        -9268.0_dp*func(i,j,3) +1008.0_dp*func(i,j,4) &
                        +5670.0_dp*func(i,j,5) -4144.0_dp*func(i,j,6) &
                        +1764.0_dp*func(i,j,7) -432.0_dp*func(i,j,8) &
                        +47.0_dp*func(i,j,9) ) / (5040.0_dp*hz*hz)

                    else if (k == 4) then

                        second_derivative(i,j,k) = ( &
                        47.0_dp*func(i,j,1) -684.0_dp*func(i,j,2) &
                        +7308.0_dp*func(i,j,3) -13216.0_dp*func(i,j,4) &
                        +6930.0_dp*func(i,j,5) -252.0_dp*func(i,j,6) &
                        -196.0_dp*func(i,j,7) +72.0_dp*func(i,j,8) &
                        -9.0_dp*func(i,j,9) ) / (5040.0_dp*hz*hz)

                    else if (k == nz-3) then

                        second_derivative(i,j,k) = ( &
                        47.0_dp*func(i,j,nz) -684.0_dp*func(i,j,nz-1) &
                        +7308.0_dp*func(i,j,nz-2) -13216.0_dp*func(i,j,nz-3) &
                        +6930.0_dp*func(i,j,nz-4) -252.0_dp*func(i,j,nz-5) &
                        -196.0_dp*func(i,j,nz-6) +72.0_dp*func(i,j,nz-7) &
                        -9.0_dp*func(i,j,nz-8) ) / (5040.0_dp*hz*hz)

                    else if (k == nz-2) then

                        second_derivative(i,j,k) = ( &
                        -261.0_dp*func(i,j,nz) +5616.0_dp*func(i,j,nz-1) &
                        -9268.0_dp*func(i,j,nz-2) +1008.0_dp*func(i,j,nz-3) &
                        +5670.0_dp*func(i,j,nz-4) -4144.0_dp*func(i,j,nz-5) &
                        +1764.0_dp*func(i,j,nz-6) -432.0_dp*func(i,j,nz-7) &
                        +47.0_dp*func(i,j,nz-8) ) / (5040.0_dp*hz*hz)

                    else if (k == nz-1) then

                        second_derivative(i,j,k) = ( &
                        3267.0_dp*func(i,j,nz) +128.0_dp*func(i,j,nz-1) &
                        -20916.0_dp*func(i,j,nz-2) +38556.0_dp*func(i,j,nz-3) &
                        -37030.0_dp*func(i,j,nz-4) +23688.0_dp*func(i,j,nz-5) &
                        -9828.0_dp*func(i,j,nz-6) +2396.0_dp*func(i,j,nz-7) &
                        -261.0_dp*func(i,j,nz-8) ) / (5040.0_dp*hz*hz)

                    else if (k == nz) then

                        second_derivative(i,j,k) = ( &
                        29531.0_dp*func(i,j,nz) -138528.0_dp*func(i,j,nz-1) &
                        +312984.0_dp*func(i,j,nz-2) -448672.0_dp*func(i,j,nz-3) &
                        +435330.0_dp*func(i,j,nz-4) -284256.0_dp*func(i,j,nz-5) &
                        +120008.0_dp*func(i,j,nz-6) -29664.0_dp*func(i,j,nz-7) &
                        +3267.0_dp*func(i,j,nz-8) ) / (5040.0_dp*hz*hz)

                    end if

                end do
            end do
        end do
    end subroutine dd_z_cp

    subroutine laplacian_cp(func, hx, hy, hz, lap)
        implicit none
        complex(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hx, hy, hz
        complex(dp), intent(out) :: lap(size(func,1), size(func,2), size(func,3))

        complex(dp), allocatable :: d2x(:,:,:), d2y(:,:,:), d2z(:,:,:)
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        allocate(d2x(nx,ny,nz), d2y(nx,ny,nz), d2z(nx,ny,nz))

        call dd_x_cp(func, hx, d2x)
        call dd_y_cp(func, hy, d2y)
        call dd_z_cp(func, hz, d2z)

        lap = d2x + d2y + d2z

        deallocate(d2x, d2y, d2z)
    end subroutine laplacian_cp


    subroutine d_x_cp(func, hx, derivative)
        implicit none
        complex(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hx
        complex(dp), intent(out) :: derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (nx < 9) stop "d_x: size(func,1) must be >= 9"

        ! 9 point derivative in x-direction
        if (i >= 5 .and. i <= nx-4) then

            derivative(i,j,k) = ( &
                3.0_dp*func(i-4,j,k) - 32.0_dp*func(i-3,j,k) &
            + 168.0_dp*func(i-2,j,k) - 672.0_dp*func(i-1,j,k) &
            + 672.0_dp*func(i+1,j,k) - 168.0_dp*func(i+2,j,k) &
            + 32.0_dp*func(i+3,j,k) - 3.0_dp*func(i+4,j,k) ) / (840.0_dp*hx)

        else if (i == 1) then

            derivative(i,j,k) = ( &
            -2283.0_dp*func(1,j,k) + 6720.0_dp*func(2,j,k) &
            -11760.0_dp*func(3,j,k) + 15680.0_dp*func(4,j,k) &
            -14700.0_dp*func(5,j,k) + 9408.0_dp*func(6,j,k) &
            -3920.0_dp*func(7,j,k) + 960.0_dp*func(8,j,k) &
            -105.0_dp*func(9,j,k) ) / (840.0_dp*hx)

        else if (i == 2) then

            derivative(i,j,k) = ( &
            -105.0_dp*func(1,j,k) - 1338.0_dp*func(2,j,k) &
            +2940.0_dp*func(3,j,k) - 2940.0_dp*func(4,j,k) &
            +2450.0_dp*func(5,j,k) - 1470.0_dp*func(6,j,k) &
            +588.0_dp*func(7,j,k) - 140.0_dp*func(8,j,k) &
            +15.0_dp*func(9,j,k) ) / (840.0_dp*hx)

        else if (i == 3) then

            derivative(i,j,k) = ( &
            15.0_dp*func(1,j,k) - 240.0_dp*func(2,j,k) &
            -798.0_dp*func(3,j,k) + 1680.0_dp*func(4,j,k) &
            -1050.0_dp*func(5,j,k) + 560.0_dp*func(6,j,k) &
            -210.0_dp*func(7,j,k) + 48.0_dp*func(8,j,k) &
            -5.0_dp*func(9,j,k) ) / (840.0_dp*hx)

        else if (i == 4) then

            derivative(i,j,k) = ( &
            -5.0_dp*func(1,j,k) + 60.0_dp*func(2,j,k) &
            -420.0_dp*func(3,j,k) - 378.0_dp*func(4,j,k) &
            +1050.0_dp*func(5,j,k) - 420.0_dp*func(6,j,k) &
            +140.0_dp*func(7,j,k) - 30.0_dp*func(8,j,k) &
            +3.0_dp*func(9,j,k) ) / (840.0_dp*hx)

        else if (i == nx-3) then

            derivative(i,j,k) = -( &
            -5.0_dp*func(nx,j,k) + 60.0_dp*func(nx-1,j,k) &
            -420.0_dp*func(nx-2,j,k) - 378.0_dp*func(nx-3,j,k) &
            +1050.0_dp*func(nx-4,j,k) - 420.0_dp*func(nx-5,j,k) &
            +140.0_dp*func(nx-6,j,k) - 30.0_dp*func(nx-7,j,k) &
            +3.0_dp*func(nx-8,j,k) ) / (840.0_dp*hx)

        else if (i == nx-2) then

            derivative(i,j,k) = -( &
            15.0_dp*func(nx,j,k) - 240.0_dp*func(nx-1,j,k) &
            -798.0_dp*func(nx-2,j,k) + 1680.0_dp*func(nx-3,j,k) &
            -1050.0_dp*func(nx-4,j,k) + 560.0_dp*func(nx-5,j,k) &
            -210.0_dp*func(nx-6,j,k) + 48.0_dp*func(nx-7,j,k) &
            -5.0_dp*func(nx-8,j,k) ) / (840.0_dp*hx)

        else if (i == nx-1) then

            derivative(i,j,k) = -( &
            -105.0_dp*func(nx,j,k) - 1338.0_dp*func(nx-1,j,k) &
            +2940.0_dp*func(nx-2,j,k) - 2940.0_dp*func(nx-3,j,k) &
            +2450.0_dp*func(nx-4,j,k) - 1470.0_dp*func(nx-5,j,k) &
            +588.0_dp*func(nx-6,j,k) - 140.0_dp*func(nx-7,j,k) &
            +15.0_dp*func(nx-8,j,k) ) / (840.0_dp*hx)

        else if (i == nx) then

            derivative(i,j,k) = -( &
            -2283.0_dp*func(nx,j,k) + 6720.0_dp*func(nx-1,j,k) &
            -11760.0_dp*func(nx-2,j,k) + 15680.0_dp*func(nx-3,j,k) &
            -14700.0_dp*func(nx-4,j,k) + 9408.0_dp*func(nx-5,j,k) &
            -3920.0_dp*func(nx-6,j,k) + 960.0_dp*func(nx-7,j,k) &
            -105.0_dp*func(nx-8,j,k) ) / (840.0_dp*hx)

        end if
    end subroutine d_x_cp

    subroutine d_y_cp(func, hy, derivative)
        implicit none
        complex(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hy
        complex(dp), intent(out) :: derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (ny < 9) stop "d_y: size(func,2) must be >= 9"

        if (j >= 5 .and. j <= ny-4) then

            derivative(i,j,k) = ( &
                3.0_dp*func(i,j-4,k) - 32.0_dp*func(i,j-3,k) &
            + 168.0_dp*func(i,j-2,k) - 672.0_dp*func(i,j-1,k) &
            + 672.0_dp*func(i,j+1,k) - 168.0_dp*func(i,j+2,k) &
            + 32.0_dp*func(i,j+3,k) - 3.0_dp*func(i,j+4,k) ) / (840.0_dp*hy)

        else if (j == 1) then

            derivative(i,j,k) = ( &
            -2283.0_dp*func(i,1,k) + 6720.0_dp*func(i,2,k) &
            -11760.0_dp*func(i,3,k) + 15680.0_dp*func(i,4,k) &
            -14700.0_dp*func(i,5,k) + 9408.0_dp*func(i,6,k) &
            -3920.0_dp*func(i,7,k) + 960.0_dp*func(i,8,k) &
            -105.0_dp*func(i,9,k) ) / (840.0_dp*hy)

        else if (j == 2) then

            derivative(i,j,k) = ( &
            -105.0_dp*func(i,1,k) - 1338.0_dp*func(i,2,k) &
            +2940.0_dp*func(i,3,k) - 2940.0_dp*func(i,4,k) &
            +2450.0_dp*func(i,5,k) - 1470.0_dp*func(i,6,k) &
            +588.0_dp*func(i,7,k) - 140.0_dp*func(i,8,k) &
            +15.0_dp*func(i,9,k) ) / (840.0_dp*hy)

        else if (j == 3) then

            derivative(i,j,k) = ( &
            15.0_dp*func(i,j-4,k) - 240.0_dp*func(i,j-3,k) &
            -798.0_dp*func(i,j-2,k) + 1680.0_dp*func(i,j-1,k) &
            -1050.0_dp*func(i,j+1,k) + 560.0_dp*func(i,j+2,k) &
            -210.0_dp*func(i,j+3,k) + 48.0_dp*func(i,j+4,k) &
            -5.0_dp*func(i,j+5,k) ) / (840.0_dp*hy)

        else if (j == 4) then

            derivative(i,j,k) = ( &
            -5.0_dp*func(i,1,k) + 60.0_dp*func(i,2,k) &
            -420.0_dp*func(i,3,k) - 378.0_dp*func(i,4,k) &
            +1050.0_dp*func(i,5,k) - 420.0_dp*func(i,6,k) &
            +140.0_dp*func(i,7,k) - 30.0_dp*func(i,8,k) &
            +3.0_dp*func(i,9,k) ) / (840.0_dp*hy)

        else if (j == ny-3) then

            derivative(i,j,k) = -( &
            -5.0_dp*func(i,ny,k) + 60.0_dp*func(i,ny-1,k) &
            -420.0_dp*func(i,ny-2,k) - 378.0_dp*func(i,ny-3,k) &
            +1050.0_dp*func(i,ny-4,k) - 420.0_dp*func(i,ny-5,k) &
            +140.0_dp*func(i,ny-6,k) - 30.0_dp*func(i,ny-7,k) &
            +3.0_dp*func(i,ny-8,k) ) / (840.0_dp*hy)

        else if (j == ny-2) then

            derivative(i,j,k) = -( &
            15.0_dp*func(i,ny,k) - 240.0_dp*func(i,ny-1,k) &
            -798.0_dp*func(i,ny-2,k) + 1680.0_dp*func(i,ny-3,k) &
            -1050.0_dp*func(i,ny-4,k) + 560.0_dp*func(i,ny-5,k) &
            -210.0_dp*func(i,ny-6,k) + 48.0_dp*func(i,ny-7,k) &
            -5.0_dp*func(i,ny-8,k) ) / (840.0_dp*hy)

        else if (j == ny-1) then

            derivative(i,j,k) = -( &
            -105.0_dp*func(i,1,k) - 1338.0_dp*func(i,2,k) &
            +2940.0_dp*func(i,3,k) - 2940.0_dp*func(i,4,k) &
            +2450.0_dp*func(i,5,k) - 1470.0_dp*func(i,6,k) &
            +588.0_dp*func(i,7,k) - 140.0_dp*func(i,8,k) &
            +15.0_dp*func(i,9,k) ) / (840.0_dp*hy)

        else if (j == ny) then

            derivative(i,j,k) = -( &
            -2283.0_dp*func(i,ny,k) + 6720.0_dp*func(i,ny-1,k) &
            -11760.0_dp*func(i,ny-2,k) + 15680.0_dp*func(i,ny-3,k) &
            -14700.0_dp*func(i,ny-4,k) + 9408.0_dp*func(i,ny-5,k) &
            -3920.0_dp*func(i,ny-6,k) + 960.0_dp*func(i,ny-7,k) &
            -105.0_dp*func(i,ny-8,k) ) / (840.0_dp*hy)

        end if
    end subroutine d_y_cp

    subroutine d_z_cp(func, hz, derivative)
        implicit none
        complex(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hz
        complex(dp), intent(out) :: derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (nz < 9) stop "d_z: size(func,3) must be >= 9"

        if (k >= 5 .and. k <= nz-4) then

            derivative(i,j,k) = ( &
                3.0_dp*func(i,j,k-4) - 32.0_dp*func(i,j,k-3) &
            + 168.0_dp*func(i,j,k-2) - 672.0_dp*func(i,j,k-1) &
            + 672.0_dp*func(i,j,k+1) - 168.0_dp*func(i,j,k+2) &
            + 32.0_dp*func(i,j,k+3) - 3.0_dp*func(i,j,k+4) ) / (840.0_dp*hz)

        else if (k == 1) then

            derivative(i,j,k) = ( &
            -2283.0_dp*func(i,j,1) + 6720.0_dp*func(i,j,2) &
            -11760.0_dp*func(i,j,3) + 15680.0_dp*func(i,j,4) &
            -14700.0_dp*func(i,j,5) + 9408.0_dp*func(i,j,6) &
            -3920.0_dp*func(i,j,7) + 960.0_dp*func(i,j,8) &
            -105.0_dp*func(i,j,9) ) / (840.0_dp*hz)

        else if (k == 2) then

            derivative(i,j,k) = ( &
            -105.0_dp*func(i,j,1) - 1338.0_dp*func(i,j,2) &
            +2940.0_dp*func(i,j,3) - 2940.0_dp*func(i,j,4) &
            +2450.0_dp*func(i,j,5) - 1470.0_dp*func(i,j,6) &
            +588.0_dp*func(i,j,7) - 140.0_dp*func(i,j,8) &
            +15.0_dp*func(i,j,9) ) / (840.0_dp*hz)

        else if (k == 3) then

            derivative(i,j,k) = ( &
            15.0_dp*func(i,j,1) - 240.0_dp*func(i,j,2) &
            -798.0_dp*func(i,j,3) + 1680.0_dp*func(i,j,4) &
            -1050.0_dp*func(i,j,5) + 560.0_dp*func(i,j,6) &
            -210.0_dp*func(i,j,7) + 48.0_dp*func(i,j,8) &
            -5.0_dp*func(i,j,9) ) / (840.0_dp*hz)

        else if (k == 4) then

            derivative(i,j,k) = ( &
            -5.0_dp*func(i,j,1) + 60.0_dp*func(i,j,2) &
            -420.0_dp*func(i,j,3) - 378.0_dp*func(i,j,4) &
            +1050.0_dp*func(i,j,5) - 420.0_dp*func(i,j,6) &
            +140.0_dp*func(i,j,7) - 30.0_dp*func(i,j,8) &
            +3.0_dp*func(i,j,9) ) / (840.0_dp*hz)

        else if (k == nz-3) then

            derivative(i,j,k) = -( &
            -5.0_dp*func(i,j,nz) + 60.0_dp*func(i,j,nz-1) &
            -420.0_dp*func(i,j,nz-2) - 378.0_dp*func(i,j,nz-3) &
            +1050.0_dp*func(i,j,nz-4) - 420.0_dp*func(i,j,nz-5) &
            +140.0_dp*func(i,j,nz-6) - 30.0_dp*func(i,j,nz-7) &
            +3.0_dp*func(i,j,nz-8) ) / (840.0_dp*hz)

        else if (k == nz-2) then

            derivative(i,j,k) = -( &
            15.0_dp*func(i,j,nz) - 240.0_dp*func(i,j,nz-1) &
            -798.0_dp*func(i,j,nz-2) + 1680.0_dp*func(i,j,nz-3) &
            -1050.0_dp*func(i,j,nz-4) + 560.0_dp*func(i,j,nz-5) &
            -210.0_dp*func(i,j,nz-6) + 48.0_dp*func(i,j,nz-7) &
            -5.0_dp*func(i,j,nz-8) ) / (840.0_dp*hz)

        else if (k == nz-1) then

            derivative(i,j,k) = -( &
            -105.0_dp*func(i,j,nz) - 1338.0_dp*func(i,j,nz-1) &
            +2940.0_dp*func(i,j,nz-2) - 2940.0_dp*func(i,j,nz-3) &
            +2450.0_dp*func(i,j,nz-4) - 1470.0_dp*func(i,j,nz-5) &
      +588.0_dp*func(i,j,nz-6) - 140.0_dp*func(i,j,nz-7) &
      +15.0_dp*func(i,j,nz-8) ) / (840.0_dp*hz)

        else if (k == nz) then

            derivative(i,j,k) = -( &
            -2283.0_dp*func(i,j,nz) + 6720.0_dp*func(i,j,nz-1) &
            -11760.0_dp*func(i,j,nz-2) + 15680.0_dp*func(i,j,nz-3) &
            -14700.0_dp*func(i,j,nz-4) + 9408.0_dp*func(i,j,nz-5) &
            -3920.0_dp*func(i,j,nz-6) + 960.0_dp*func(i,j,nz-7) &
            -105.0_dp*func(i,j,nz-8) ) / (840.0_dp*hz)

        end if
    end subroutine d_z_cp

end module three_D_derivative_cp_mod