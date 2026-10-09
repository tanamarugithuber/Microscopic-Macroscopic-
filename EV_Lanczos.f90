module EV_Lanczos
    use iso_fortran_env, only: real64
    use common_func_mod
    use constant_mod
    use nucleus_mod
    use micro_constant_mod, only: microscopic_variables
    use micro_potential_mod
    implicit none
    private

    integer, parameter :: dp = real64

    contains

        subroutine lanczos_algorithm(nx,ny,nz,l,max_iter, tol, eigenvalue, eigenvector)
            implicit none
            integer, intent(in) :: max_iter
            real(dp), intent(in) :: tol
            real(dp), intent(out) :: eigenvalue
            real(dp), intent(out) :: eigenvector(:)
            real(dp), allocatable :: b(:,:,:,:),b_new(:,:,:,:), b_old(:,:,:,:)
            real(dp), allocatable :: alpha(:,:,:,:), beta(:,:,:,:)
            integer :: nx, ny, nz, l
            ! Implement the Lanczos algorithm here to compute the largest eigenvalue and corresponding eigenvector of matrix A.
            ! Use the initial vector b as the starting point for the iterations

            ! set up initial vector b and matrix A
            allocate(b(nx,ny,nz,l))
            allocate(b_new(nx,ny,nz,l))
            allocate(b_old(nx,ny,nz,l))
            allocate(alpha(nx,ny,nz,l))
            allocate(beta(nx,ny,nz,l))

            ! Initialize b with some values 
            ! use a random number generator


        end subroutine lanczos_algorithm

end module EV_Lanczos