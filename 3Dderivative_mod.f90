module three_D_derivative_mod
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter :: dp = real64

    public :: laplacian

contains

    subroutine d_x(func, hx, derivative)
        implicit none
        real(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hx
        real(dp), intent(out) :: derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (nx < 7) stop "d_x: size(func,1) must be >= 7"

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (i >= 4 .and. i <= nx-3) then
                        derivative(i,j,k) = ( &
                             func(i+3,j,k) - 9.0_dp*func(i+2,j,k) + 45.0_dp*func(i+1,j,k) &
                           - 45.0_dp*func(i-1,j,k) + 9.0_dp*func(i-2,j,k) - func(i-3,j,k) ) / (60.0_dp*hx)
                    else if (i <= 3) then
                        derivative(i,j,k) = (func(i+1,j,k) - func(i,j,k)) / hx
                    else
                        derivative(i,j,k) = (func(i,j,k) - func(i-1,j,k)) / hx
                    end if
                end do
            end do
        end do
    end subroutine d_x

    subroutine d_y(func, hy, derivative)
        implicit none
        real(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hy
        real(dp), intent(out) :: derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (ny < 7) stop "d_y: size(func,2) must be >= 7"

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (j >= 4 .and. j <= ny-3) then
                        derivative(i,j,k) = ( &
                             func(i,j+3,k) - 9.0_dp*func(i,j+2,k) + 45.0_dp*func(i,j+1,k) &
                           - 45.0_dp*func(i,j-1,k) + 9.0_dp*func(i,j-2,k) - func(i,j-3,k) ) / (60.0_dp*hy)
                    else if (j <= 3) then
                        derivative(i,j,k) = (func(i,j+1,k) - func(i,j,k)) / hy
                    else
                        derivative(i,j,k) = (func(i,j,k) - func(i,j-1,k)) / hy
                    end if
                end do
            end do
        end do
    end subroutine d_y

    subroutine d_z(func, hz, derivative)
        implicit none
        real(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hz
        real(dp), intent(out) :: derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (nz < 7) stop "d_z: size(func,3) must be >= 7"

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (k >= 4 .and. k <= nz-3) then
                        derivative(i,j,k) = ( &
                             func(i,j,k+3) - 9.0_dp*func(i,j,k+2) + 45.0_dp*func(i,j,k+1) &
                           - 45.0_dp*func(i,j,k-1) + 9.0_dp*func(i,j,k-2) - func(i,j,k-3) ) / (60.0_dp*hz)
                    else if (k <= 3) then
                        derivative(i,j,k) = (func(i,j,k+1) - func(i,j,k)) / hz
                    else
                        derivative(i,j,k) = (func(i,j,k) - func(i,j,k-1)) / hz
                    end if
                end do
            end do
        end do
    end subroutine d_z

    subroutine dd_x(func, hx, second_derivative)
        implicit none
        real(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hx
        real(dp), intent(out) :: second_derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (nx < 7) stop "dd_x: size(func,1) must be >= 7"

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (i >= 4 .and. i <= nx-3) then
                        second_derivative(i,j,k) = ( &
                             2.0_dp*func(i+3,j,k) - 27.0_dp*func(i+2,j,k) + 270.0_dp*func(i+1,j,k) &
                           - 490.0_dp*func(i,j,k) + 270.0_dp*func(i-1,j,k) - 27.0_dp*func(i-2,j,k) &
                           + 2.0_dp*func(i-3,j,k) ) / (180.0_dp*hx*hx)
                    else if (i <= 3) then
                        second_derivative(i,j,k) = (func(i+2,j,k) - 2.0_dp*func(i+1,j,k) + func(i,j,k)) / (hx*hx)
                    else
                        second_derivative(i,j,k) = (func(i,j,k) - 2.0_dp*func(i-1,j,k) + func(i-2,j,k)) / (hx*hx)
                    end if
                end do
            end do
        end do
    end subroutine dd_x

    subroutine dd_y(func, hy, second_derivative)
        implicit none
        real(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hy
        real(dp), intent(out) :: second_derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (ny < 7) stop "dd_y: size(func,2) must be >= 7"

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (j >= 4 .and. j <= ny-3) then
                        second_derivative(i,j,k) = ( &
                             2.0_dp*func(i,j+3,k) - 27.0_dp*func(i,j+2,k) + 270.0_dp*func(i,j+1,k) &
                           - 490.0_dp*func(i,j,k) + 270.0_dp*func(i,j-1,k) - 27.0_dp*func(i,j-2,k) &
                           + 2.0_dp*func(i,j-3,k) ) / (180.0_dp*hy*hy)
                    else if (j <= 3) then
                        second_derivative(i,j,k) = (func(i,j+2,k) - 2.0_dp*func(i,j+1,k) + func(i,j,k)) / (hy*hy)
                    else
                        second_derivative(i,j,k) = (func(i,j,k) - 2.0_dp*func(i,j-1,k) + func(i,j-2,k)) / (hy*hy)
                    end if
                end do
            end do
        end do
    end subroutine dd_y

    subroutine dd_z(func, hz, second_derivative)
        implicit none
        real(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hz
        real(dp), intent(out) :: second_derivative(size(func,1), size(func,2), size(func,3))

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        if (nz < 7) stop "dd_z: size(func,3) must be >= 7"

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (k >= 4 .and. k <= nz-3) then
                        second_derivative(i,j,k) = ( &
                             2.0_dp*func(i,j,k+3) - 27.0_dp*func(i,j,k+2) + 270.0_dp*func(i,j,k+1) &
                           - 490.0_dp*func(i,j,k) + 270.0_dp*func(i,j,k-1) - 27.0_dp*func(i,j,k-2) &
                           + 2.0_dp*func(i,j,k-3) ) / (180.0_dp*hz*hz)
                    else if (k <= 3) then
                        second_derivative(i,j,k) = (func(i,j,k+2) - 2.0_dp*func(i,j,k+1) + func(i,j,k)) / (hz*hz)
                    else
                        second_derivative(i,j,k) = (func(i,j,k) - 2.0_dp*func(i,j,k-1) + func(i,j,k-2)) / (hz*hz)
                    end if
                end do
            end do
        end do
    end subroutine dd_z

    subroutine laplacian(func, hx, hy, hz, lap)
        implicit none
        real(dp), intent(in)  :: func(:,:,:)
        real(dp), intent(in)  :: hx, hy, hz
        real(dp), intent(out) :: lap(size(func,1), size(func,2), size(func,3))

        real(dp), allocatable :: d2x(:,:,:), d2y(:,:,:), d2z(:,:,:)
        integer :: nx, ny, nz

        nx = size(func,1)
        ny = size(func,2)
        nz = size(func,3)

        allocate(d2x(nx,ny,nz), d2y(nx,ny,nz), d2z(nx,ny,nz))

        call dd_x(func, hx, d2x)
        call dd_y(func, hy, d2y)
        call dd_z(func, hz, d2z)

        lap = d2x + d2y + d2z

        deallocate(d2x, d2y, d2z)
    end subroutine laplacian

end module three_D_derivative_mod