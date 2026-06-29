module calculate_vector_mod
    use iso_fortran_env, only: real64
    use grid_mod, only: grid_type
    implicit none
    private

    integer, parameter :: dp = real64

    contains
        subroutine calculate_coulomb_FFT(rho, V_c, grid)
            implicit none
            class(grid_type), intent(in) :: grid
            real(dp), intent(in) :: rho(:,:,:)
            real(dp), intent(out) :: V_c(:,:,:)
            integer :: i, j, k
            integer :: nx, ny, nz
            nx = size(rho, 1)
            ny = size(rho, 2)
            nz = size(rho, 3)
            
        end subroutine calculate_coulomb_FFT

        subroutine calculate_Yukawa(rho, V_yukawa, grid)
            implicit none
            class(grid_type), intent(in) :: grid
            real(dp), intent(in) :: rho(:,:,:)
            real(dp), intent(out) :: V_yukawa(:,:,:)
            integer :: i, j, k
            integer :: nx, ny, nz
            nx = size(rho, 1)
            ny = size(rho, 2)
            nz = size(rho, 3)
            
        end subroutine calculate_Yukawa


end module calculate_vector_mod