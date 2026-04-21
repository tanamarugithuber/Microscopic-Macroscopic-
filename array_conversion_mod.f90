module array_conversion_mod
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter :: dp = real64

    public :: reshape_3D_to_1D
    public :: reshape_1D_to_3D
    public :: reshape_2D_to_1D
    public :: reshape_1D_to_2D

contains
    subroutine reshape_3D_to_1D(func_3D, func_1D)
        implicit none
        real(dp), intent(in)  :: func_3D(:,:,:)
        real(dp), intent(out) :: func_1D(:)

        integer :: i, j, k
        integer :: nx, ny, nz
        integer :: index

        nx = size(func_3D, 1)
        ny = size(func_3D, 2)
        nz = size(func_3D, 3)

        if (size(func_1D) /= nx*ny*nz) stop "reshape_3D_to_1D: size of func_1D must be equal to nx*ny*nz"

        index = 0
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    index = index + 1
                    func_1D(index) = func_3D(i,j,k)
                end do
            end do
        end do
    end subroutine reshape_3D_to_1D

    subroutine reshape_1D_to_3D(func_1D, func_3D)
        implicit none
        real(dp), intent(in)  :: func_1D(:)
        real(dp), intent(out) :: func_3D(:,:,:)

        integer :: i, j, k
        integer :: nx, ny, nz
        integer :: index

        nx = size(func_3D, 1)
        ny = size(func_3D, 2)
        nz = size(func_3D, 3)

        if (size(func_1D) /= nx*ny*nz) stop "reshape_1D_to_3D: size of func_1D must be equal to nx*ny*nz"

        index = 0
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    index = index + 1
                    func_3D(i,j,k) = func_1D(index)
                end do
            end do
        end do
    end subroutine reshape_1D_to_3D

    subroutine reshape_2D_to_1D(func_2D, func_1D)
        implicit none
        real(dp), intent(in)  :: func_2D(:,:)
        real(dp), intent(out) :: func_1D(:)

        integer :: i, j
        integer :: nx, ny
        integer :: index

        nx = size(func_2D, 1)
        ny = size(func_2D, 2)

        if (size(func_1D) /= nx*ny) stop "reshape_2D_to_1D: size of func_1D must be equal to nx*ny"

        index = 0
        do j = 1, ny
            do i = 1, nx
                index = index + 1
                func_1D(index) = func_2D(i,j)
            end do
        end do
    end subroutine reshape_2D_to_1D

    subroutine reshape_1D_to_2D(func_1D, func_2D)
        implicit none
        real(dp), intent(in)  :: func_1D(:)
        real(dp), intent(out) :: func_2D(:,:)

        integer :: i, j
        integer :: nx, ny
        integer :: index

        nx = size(func_2D, 1)
        ny = size(func_2D, 2)

        if (size(func_1D) /= nx*ny) stop "reshape_1D_to_2D: size of func_1D must be equal to nx*ny"

        index = 0
        do j = 1, ny
            do i = 1, nx
                index = index + 1
                func_2D(i,j) = func_1D(index)
            end do
        end do
    end subroutine reshape_1D_to_2D

end module array_conversion_mod

