module CG_method_mod
    use iso_fortran_env, only: real64

    implicit none
    private
    integer, parameter :: dp = real64
    real(dp), allocatable, public :: helmholtz_matrix(:,:)
    real(dp), allocatable, public :: poisson_matrix(:,:)

    !public :: CG_method_1 ! Using the mat A of Ax = b explicitly
    public :: CG_method ! Using the mat A of Ax = b implicitly, i.e., using the function to calculate Ax.
    public :: mat_to_vec
    public :: initialize_helmholtz_matrix
    public :: initialize_poisson_matrix

    ! public :: index_3D_to_1D
    


contains
    subroutine mat_to_vec(A, x) ! reshape the 3D array A to a 1D array x, where A is the discretized laplacian operator in 3D, size is (n_x*n_y*n_z, n_x*n_y*n_z), and x is the vectorized form of A, size is (n_x*n_y*n_z * n_x*n_y*n_z)
        implicit none
        real(dp), intent(in) :: A(:,:,:)
        real(dp), intent(out) :: x(:)
        
        integer :: n_x, n_y, n_z
        n_x = size(A, 1)
        n_y = size(A, 2)
        n_z = size(A, 3)
        if (size(x) /= n_x*n_y*n_z) stop "mat_to_vec: size of x must be equal to n_x*n_y*n_z"
        x = reshape(A, [n_x*n_y*n_z]) 
    end subroutine mat_to_vec

    ! subroutine index_3D_to_1D( n_x, n_y, n_z, index)
    !     implicit none
    !     integer :: i, j, k
    !     integer, intent(in) :: n_x, n_y, n_z
    !     integer, intent(out) :: index(n_x,n_y,n_z) ! index(i,j,k) gives the corresponding index in the 1D array for the 3D grid point (i,j,k)

    !     if (i < 1 .or. i > n_x .or. j < 1 .or. j > n_y .or. k < 1 .or. k > n_z) then
    !         stop "index_3D_to_1D: i, j, k must be within the bounds of the grid"
    !     end if
    !     do k = 1, n_z
    !         do j = 1, n_y
    !             do i = 1, n_x
    !                 index(i,j,k) = (k-1)*n_x*n_y + (j-1)*n_x + i
    !             end do
    !         end do
    !     end do
    ! end subroutine index_3D_to_1D

    subroutine initialize_helmholtz_matrix(n_x, h_x, coeff, A)
        implicit none
        integer, intent(in) :: n_x
        real(dp), intent(in) :: h_x, coeff
        real(dp), intent(out) :: A(n_x**3, n_x**3)
        
        integer :: i, j, k, row
        real(dp) :: h2inv

        print *, "Initializing Helmholtz matrix for the grid..."
        
        h2inv = 1.0_dp / (12.0_dp * h_x * h_x)
        A = 0.0_dp

        ! 3次元格子でループを回す
        do k = 1, n_x
            do j = 1, n_x
                do i = 1, n_x
                    ! 現在のグリッド点の1次元インデックス（行番号）
                    row = (k-1)*n_x*n_x + (j-1)*n_x + i
                    
                    ! --- 対角成分 ---
                    A(row, row) = -90.0_dp * h2inv - 1.0_dp / coeff**2
                    
                    ! --- X方向の隣接点 ---
                    if (i+1 <= n_x) A(row + 1, row)     = 16.0_dp * h2inv
                    if (i-1 >= 1)   A(row - 1, row)     = 16.0_dp * h2inv
                    if (i+2 <= n_x) A(row + 2, row)     = -1.0_dp * h2inv
                    if (i-2 >= 1)   A(row - 2, row)     = -1.0_dp * h2inv
                    
                    ! --- Y方向の隣接点 ---
                    if (j+1 <= n_x) A(row + n_x, row)   = 16.0_dp * h2inv
                    if (j-1 >= 1)   A(row - n_x, row)   = 16.0_dp * h2inv
                    if (j+2 <= n_x) A(row + 2*n_x, row) = -1.0_dp * h2inv
                    if (j-2 >= 1)   A(row - 2*n_x, row) = -1.0_dp * h2inv
                    
                    ! --- Z方向の隣接点 ---
                    if (k+1 <= n_x) A(row + n_x*n_x, row)   = 16.0_dp * h2inv
                    if (k-1 >= 1)   A(row - n_x*n_x, row)   = 16.0_dp * h2inv
                    if (k+2 <= n_x) A(row + 2*n_x*n_x, row) = -1.0_dp * h2inv
                    if (k-2 >= 1)   A(row - 2*n_x*n_x, row) = -1.0_dp * h2inv
                end do
            end do
        end do
        print *, "Helmholtz matrix initialized."
    end subroutine initialize_helmholtz_matrix

    subroutine initialize_poisson_matrix(n_x, h_x, A)
        implicit none
        integer, intent(in) :: n_x
        real(dp), intent(in) :: h_x
        real(dp), intent(out) :: A(n_x**3, n_x**3)
        
        integer :: i, j, k, row
        real(dp) :: h2inv
        print *, "Initializing Poisson matrix for the grid..."
        h2inv = 1.0_dp / (h_x * h_x)
        A = 0.0_dp

        ! 3次元格子でループを回す
        do k = 1, n_x
            do j = 1, n_x
                do i = 1, n_x
                    ! 現在のグリッド点の1次元インデックス（行番号）
                    row = (k-1)*n_x*n_x + (j-1)*n_x + i
                    
                    ! --- 対角成分 ---
                    A(row, row) = -6.0_dp * h2inv
                    
                    ! --- X方向の隣接点 ---
                    if (i+1 <= n_x) A(row + 1, row)     = h2inv
                    if (i-1 >= 1)   A(row - 1, row)     = h2inv
                    
                    ! --- Y方向の隣接点 ---
                    if (j+1 <= n_x) A(row + n_x, row)   = h2inv
                    if (j-1 >= 1)   A(row - n_x, row)   = h2inv
                    
                    ! --- Z方向の隣接点 ---
                    if (k+1 <= n_x) A(row + n_x*n_x, row)   = h2inv
                    if (k-1 >= 1)   A(row - n_x*n_x, row)   = h2inv
                end do
            end do
        end do
        print *, "Poisson matrix initialized."
    end subroutine initialize_poisson_matrix

    
    subroutine CG_method(A,b,x)
        implicit none
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: tol ! tolerance for convergence
        integer, parameter :: max_iter = 10000
        integer  :: iter
        real(dp), allocatable :: temp(:), r(:), p(:)
        real(dp) :: alpha, beta

        integer :: n
        print *, "Starting Conjugate Gradient method..."
        n = size(b)

        tol = 1.0e-20_dp
        allocate(temp(n), r(n), p(n))
        if (size(A,1) /= n .or. size(A,2) /= n) stop "CG_method: size of A must be (n,n)"

        x(:) = 0.0_dp ! Initial guess

        r(:) = b(:) - matmul(A, x(:)) ! Initial residual
        p(:) = r(:) ! Initial search direction

        iter = 1
        do iter = 1, max_iter
            temp(:) = matmul(A, p(:)) ! A*p
            alpha = dot_product(r(:), r(:)) / dot_product(p(:), temp(:)) ! step size
            x(:) = x(:) + alpha * p(:) ! update solution
            r(:) = r(:) - alpha * temp(:) ! update residual

            if (sqrt(dot_product(r(:), r(:))) < tol) then
                print *, "CG converged in ", iter, " iterations."
                exit
            end if

            beta = dot_product(r(:), r(:)) / dot_product(p(:), temp(:)) ! update beta
            p(:) = r(:) + beta * p(:) ! update search direction

            if (mod(iter, 100) == 0) then
                print *, "CG iteration ", iter, ": residual norm = ", sqrt(dot_product(r(:), r(:)))
            end if
        end do

        if (iter > max_iter) then
            print *, "CG did not converge within the maximum number of iterations."
        end if
        print *, "Conjugate Gradient method finished."
        print *, "Final residual norm: ", sqrt(dot_product(r(:), r(:)))
        deallocate(temp, r, p)
    end subroutine CG_method

    
    



end module CG_method_mod