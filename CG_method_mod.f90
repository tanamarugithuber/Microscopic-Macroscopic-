module CG_method_mod
    !$ use omp_lib
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
    public :: CG_method_helmholtz_EQ
    public :: helmholtz_x
    public :: poisson_x
    public :: CG_method_poisson

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
        
        h2inv = 1.0_dp / (5040.0_dp * h_x * h_x)
        A = 0.0_dp

        ! loop over the 3Dgrid points
        ! 9points finite difference formula for the second derivative in 3D
        do k = 1, n_x
            do j = 1, n_x
                do i = 1, n_x
                    ! current grid point's 1D index (row number)
                    row = (k-1)*n_x*n_x + (j-1)*n_x + i
                    
                    ! --- diagonal element ---
                    A(row, row) = -43050.0_dp * h2inv - 1.0_dp / coeff**2
                    
                    ! --- X direction neighbors ---
                    if (i+1 <= n_x) A(row + 1, row)     = 8064.0_dp * h2inv
                    if (i-1 >= 1)   A(row - 1, row)     = 8064.0_dp * h2inv
                    if (i+2 <= n_x) A(row + 2, row)     = -1008.0_dp * h2inv
                    if (i-2 >= 1)   A(row - 2, row)     = -1008.0_dp * h2inv   
                    if (i+3 <= n_x) A(row + 3, row)     = 128.0_dp * h2inv
                    if (i-3 >= 1)   A(row - 3, row)     = 128.0_dp * h2inv
                    if (i+4 <= n_x) A(row + 4, row)     = -9.0_dp * h2inv
                    if (i-4 >= 1)   A(row - 4, row)     = -9.0_dp * h2inv

                    ! --- Y direction neighbors ---
                    if (j+1 <= n_x) A(row + n_x, row)   = 8064.0_dp * h2inv
                    if (j-1 >= 1)   A(row - n_x, row)   = 8064.0_dp * h2inv
                    if (j+2 <= n_x) A(row + 2*n_x, row) = -1008.0_dp * h2inv
                    if (j-2 >= 1)   A(row - 2*n_x, row) = -1008.0_dp * h2inv
                    if (j+3 <= n_x) A(row + 3*n_x, row) = 128.0_dp * h2inv
                    if (j-3 >= 1)   A(row - 3*n_x, row) = 128.0_dp * h2inv
                    if (j+4 <= n_x) A(row + 4*n_x, row) = -9.0_dp * h2inv
                    if (j-4 >= 1)   A(row - 4*n_x, row) = -9.0_dp * h2inv

                    ! --- Z direction neighbors ---
                    if (k+1 <= n_x) A(row + n_x*n_x, row)   = 8064.0_dp * h2inv
                    if (k-1 >= 1)   A(row - n_x*n_x, row)   = 8064.0_dp * h2inv
                    if (k+2 <= n_x) A(row + 2*n_x*n_x, row) = -1008.0_dp * h2inv
                    if (k-2 >= 1)   A(row - 2*n_x*n_x, row) = -1008.0_dp * h2inv
                    if (k+3 <= n_x) A(row + 3*n_x*n_x, row) = 128.0_dp * h2inv
                    if (k-3 >= 1)   A(row - 3*n_x*n_x, row) = 128.0_dp * h2inv
                    if (k+4 <= n_x) A(row + 4*n_x*n_x, row) = -9.0_dp * h2inv
                    if (k-4 >= 1)   A(row - 4*n_x*n_x, row) = -9.0_dp * h2inv
                    
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
        h2inv = 1.0_dp / (5040.0_dp * h_x * h_x)
        A = 0.0_dp

        ! loop over the 3Dgrid points
        do k = 1, n_x
            do j = 1, n_x
                do i = 1, n_x
                    ! current grid point's 1D index (row number)
                    row = (k-1)*n_x*n_x + (j-1)*n_x + i

                    ! --- diagonal element ---
                    A(row, row) = -43050.0_dp * h2inv

                    ! --- x direction neighbors ---
                    if (i+1 <= n_x) A(row + 1, row)     = 8064.0_dp * h2inv
                    if (i-1 >= 1)   A(row - 1, row)     = 8064.0_dp * h2inv
                    if (i+2 <= n_x) A(row + 2, row)     = -1008.0_dp * h2inv
                    if (i-2 >= 1)   A(row - 2, row)     = -1008.0_dp * h2inv
                    if (i+3 <= n_x) A(row + 3, row)     = 128.0_dp * h2inv
                    if (i-3 >= 1)   A(row - 3, row)     = 128.0_dp * h2inv
                    if (i+4 <= n_x) A(row + 4, row)     = -9.0_dp * h2inv
                    if (i-4 >= 1)   A(row - 4, row)     = -9.0_dp * h2inv

                    ! --- y direction neighbors ---
                    if (j+1 <= n_x) A(row + n_x, row)   = 8064.0_dp * h2inv
                    if (j-1 >= 1)   A(row - n_x, row)   = 8064.0_dp * h2inv
                    if (j+2 <= n_x) A(row + 2*n_x, row) = -1008.0_dp * h2inv
                    if (j-2 >= 1)   A(row - 2*n_x, row) = -1008.0_dp * h2inv
                    if (j+3 <= n_x) A(row + 3*n_x, row) = 128.0_dp * h2inv
                    if (j-3 >= 1)   A(row - 3*n_x, row) = 128.0_dp * h2inv
                    if (j+4 <= n_x) A(row + 4*n_x, row) = -9.0_dp * h2inv
                    if (j-4 >= 1)   A(row - 4*n_x, row) = -9.0_dp * h2inv

                    ! --- z direction neighbors ---
                    if (k+1 <= n_x) A(row + n_x*n_x, row)   = 8064.0_dp * h2inv
                    if (k-1 >= 1)   A(row - n_x*n_x, row)   = 8064.0_dp * h2inv
                    if (k+2 <= n_x) A(row + 2*n_x*n_x, row) = -1008.0_dp * h2inv
                    if (k-2 >= 1)   A(row - 2*n_x*n_x, row) = -1008.0_dp * h2inv
                    if (k+3 <= n_x) A(row + 3*n_x*n_x, row) = 128.0_dp * h2inv
                    if (k-3 >= 1)   A(row - 3*n_x*n_x, row) = 128.0_dp * h2inv
                    if (k+4 <= n_x) A(row + 4*n_x*n_x, row) = -9.0_dp * h2inv
                    if (k-4 >= 1)   A(row - 4*n_x*n_x, row) = -9.0_dp * h2inv
                    
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
        real(dp) :: alpha, beta, rr_old, rr_new

        integer :: n
        print *, "Starting Conjugate Gradient method..."
        n = size(b)

        tol = 1.0e-20_dp
        allocate(temp(n), r(n), p(n))
        if (size(A,1) /= n .or. size(A,2) /= n) stop "CG_method: size of A must be (n,n)"

        x(:) = 0.0_dp ! Initial guess

        r(:) = b(:) - matmul(A, x(:)) ! Initial residual
        p(:) = r(:) ! Initial search direction
        rr_old = dot_product(r(:), r(:)) ! Initial residual squared norm
        iter = 1
        do iter = 1, max_iter
            temp(:) = matmul(A, p(:)) ! A*p
            alpha = rr_old / dot_product(p(:), temp(:)) ! step size
            x(:) = x(:) + alpha * p(:) ! update solution
            r(:) = r(:) - alpha * temp(:) ! update residual
            rr_new = dot_product(r(:), r(:)) ! new residual squared norm

            if (sqrt(rr_new) < tol) then
                print *, "CG converged in ", iter, " iterations."
                exit
            end if

            beta = rr_new / rr_old ! update beta
            p(:) = r(:) + beta * p(:) ! update search direction

            if (mod(iter, 100) == 0) then
                print *, "CG iteration ", iter, ": residual norm = ", sqrt(rr_new)
            end if
            rr_old = rr_new ! update for next iteration
        end do

        if (iter > max_iter) then
            print *, "CG did not converge within the maximum number of iterations."
        end if
        print *, "Conjugate Gradient method finished."
        print *, "Final residual norm: ", sqrt(rr_new)
        deallocate(temp, r, p)
    end subroutine CG_method



    subroutine CG_method_helmholtz_EQ(b, x, n_x, h_x, coeff)
        !$ use omp_lib
        implicit none
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(in) :: n_x
        real(dp), intent(in) :: h_x, coeff
        real(dp), allocatable :: Ax(:)
        real(dp) :: tol ! tolerance for convergence
        integer, parameter :: max_iter = 10000
        integer  :: iter
        real(dp), allocatable :: temp(:), r(:), p(:)
        real(dp) :: alpha, beta, rr_old, rr_new, denom
        integer :: n_size

        integer :: n, n_count
        print *, "Starting Conjugate Gradient method..."
        n = size(b)

        tol = 1.0e-20_dp
        allocate(temp(n), r(n), p(n), Ax(n))

        x(:) = 0.0_dp ! Initial guess

        call helmholtz_x(x, n_x, h_x, coeff, Ax)
        r(:) = b(:) - Ax(:) ! Initial residual
        p(:) = r(:) ! Initial search direction
        rr_old = dot_product(r(:), r(:)) ! Initial residual squared norm
        iter = 1
        do iter = 1, max_iter
            call helmholtz_x(p, n_x, h_x, coeff, temp) ! temp = A*p

            n_size = size(p)
            denom = 0.0_dp
            !$omp parallel do default(none) private(n) shared(p,temp,n_size) &
            !$omp reduction(+:denom) schedule(static)
            do n = 1, n_size
                denom = denom + p(n) * temp(n)
            end do
            !$omp end parallel do
            alpha = rr_old / denom

            x(:) = x(:) + alpha * p(:) ! update solution
            
            n_size = size(r)
            rr_new = 0.0_dp
            !$omp parallel do default(none) private(n) shared(r,alpha,temp,n_size) &
            !$omp reduction(+:rr_new) schedule(static)
            do n = 1, n_size
                r(n) = r(n) - alpha * temp(n)
                rr_new = rr_new + r(n) * r(n)
            end do
            !$omp end parallel do

            if (sqrt(rr_new) < tol) then
                print *, "CG converged in ", iter, " iterations."
                exit
            end if

            beta = rr_new / rr_old ! update beta
            p(:) = r(:) + beta * p(:) ! update search direction

            if (mod(iter, 100) == 0) then
                print *, "CG iteration ", iter, ": residual norm = ", sqrt(rr_new)
            end if
            rr_old = rr_new ! update for next iteration
        end do

        if (iter > max_iter) then
            print *, "CG did not converge within the maximum number of iterations."
        end if
        print *, "Conjugate Gradient method finished."
        print *, "Final residual norm: ", sqrt(rr_new)
        deallocate(temp, r, p, Ax)

        print *, "Helmholtz equation solved."
    end subroutine CG_method_helmholtz_EQ

    subroutine helmholtz_x(x, n_x, h_x, coeff, Ax)
        !$ use omp_lib
        implicit none
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: n_x
        real(dp), intent(in) :: h_x
        real(dp), intent(in) :: coeff
        real(dp), intent(out) :: Ax(:)
        
        integer :: i, j, k, row
        real(dp) :: h2inv
        ! print *, "Calculating Ax for the Helmholtz equation..."
        if (size(x) /= n_x**3) stop "helmholtz_x: size of x must be equal to n_x^3"

        h2inv = 1.0_dp / (5040.0_dp * h_x * h_x)
        Ax = 0.0_dp

        ! loop over the 3Dgrid points
        !$omp parallel do collapse(3) default(none) &
        !$omp private(i,j,k,row) &
        !$omp shared(x,Ax,n_x,h2inv,coeff) &
        !$omp schedule(static)
        do k = 1, n_x
            do j = 1, n_x
                do i = 1, n_x
                    ! current grid point's 1D index (row number)
                    row = (k-1)*n_x*n_x + (j-1)*n_x + i

                    ! --- diagonal element ---
                    Ax(row) = -43050.0_dp * h2inv * x(row) - 1.0_dp * x(row) * coeff**2

                    ! --- x direction neighbors ---
                    if (i+1 <= n_x) Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row + 1)
                    if (i-1 >= 1)   Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row - 1)
                    if (i+2 <= n_x) Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row + 2)
                    if (i-2 >= 1)   Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row - 2)
                    if (i+3 <= n_x) Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row + 3)
                    if (i-3 >= 1)   Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row - 3)
                    if (i+4 <= n_x) Ax(row) = Ax(row) -9.0_dp * h2inv * x(row + 4)
                    if (i-4 >= 1)   Ax(row) = Ax(row) -9.0_dp * h2inv * x(row - 4)

                    ! --- y direction neighbors ---
                    if (j+1 <= n_x) Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row + n_x)
                    if (j-1 >= 1)   Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row - n_x)
                    if (j+2 <= n_x) Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row + 2*n_x)
                    if (j-2 >= 1)   Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row - 2*n_x)
                    if (j+3 <= n_x) Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row + 3*n_x)
                    if (j-3 >= 1)   Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row - 3*n_x)
                    if (j+4 <= n_x) Ax(row) = Ax(row) -9.0_dp * h2inv * x(row + 4*n_x)
                    if (j-4 >= 1)   Ax(row) = Ax(row) -9.0_dp * h2inv * x(row - 4*n_x)

                    ! --- z direction neighbors ---
                    if (k+1 <= n_x) Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row + n_x*n_x)
                    if (k-1 >= 1)   Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row - n_x*n_x)
                    if (k+2 <= n_x) Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row + 2*n_x*n_x)
                    if (k-2 >= 1)   Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row - 2*n_x*n_x)
                    if (k+3 <= n_x) Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row + 3*n_x*n_x)
                    if (k-3 >= 1)   Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row - 3*n_x*n_x)
                    if (k+4 <= n_x) Ax(row) = Ax(row) -9.0_dp * h2inv * x(row + 4*n_x*n_x)
                    if (k-4 >= 1)   Ax(row) = Ax(row) -9.0_dp * h2inv * x(row - 4*n_x*n_x)

                end do
            end do
        end do
        !$omp end parallel do
        Ax(:) = -Ax(:) ! because the matrix A in the Helmholtz equation is negative definite, we need to negate Ax to get the correct result for Ax = A*x
        ! print *, "Ax calculated for the Helmholtz equation."
    end subroutine helmholtz_x

    ! CG_method for helmholtz EQ with the dirichlet boundary condition.s, using 8th order finite difference formula 
    ! use CG method to solve the points without the points within four grid points of the boundary
    ! the points within four grid points of the boundary are sloved by gauss law
    subroutine CG_method_poisson(density, solution, n_x, h_x)
        !$ use omp_lib
        implicit none
        real(dp), intent(in) :: density(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(in) :: n_x
        real(dp), intent(in) :: h_x
        real(dp), allocatable :: x(:), b(:)
        real(dp), allocatable :: Ax(:)
        real(dp), allocatable :: temp(:), r(:), p(:)!,dist(:)
        real(dp) :: alpha, beta, rr_old, rr_new, denom
        real(dp) :: tol = 1.0e-20_dp
        integer, parameter :: max_iter = 10000
        integer :: n_sizes, i, j, k, row, iter
        real(dp) :: h2inv
        real(dp) :: total_charge, x_temp, y_temp, z_temp
        integer :: nx_mid, ny_mid, nz_mid
        integer :: n_inner, ii, jj, kk, idx
        integer :: n_count
        print *, "Solving the equation Ax = b using Conjugate Gradient method..."

        ! prepare the x(:) and b(:) vectors for CG method
        ! caution : x(:) and b(:) don't have the points within four grid points of the boundary
        n_inner = n_x - 8
        n_sizes = n_inner**3
        h2inv = 1.0_dp / (5040.0_dp * h_x * h_x)
        allocate(x(n_sizes), b(n_sizes), Ax(n_sizes), temp(n_sizes), r(n_sizes), p(n_sizes))
        if (size(density) /= n_x**3) stop "CG_method_2: size of density must be equal to n_x^3"
        if (size(solution) /= n_x**3) stop "CG_method_2: size of solution must be equal to n_x^3"

        total_charge = sum(density)
        print *, "Total charge:", total_charge

        ! dist(:) is used to store the distance from zero to the points within four grid points of the boundary, which is used to calculate the contribution from these points to b(:) using gauss law
    
        ! fill the solution(:) with the initial guess, which is the contribution from the points within four grid points of the boundary to b(:) using gauss law, and then we will update solution(:) in each iteration of CG method
            solution(:) = 0.0_dp
            nx_mid = (n_x + 1) / 2
            ny_mid = (n_x + 1) / 2
            nz_mid = (n_x + 1) / 2

            !$omp parallel do collapse(3) default(none) &
            !$omp private(i,j,k,row,x_temp,y_temp,z_temp) &
            !$omp shared(solution,n_x,h_x,total_charge,nx_mid,ny_mid,nz_mid) &
            !$omp schedule(static)
            do k = 1, n_x
                do j = 1, n_x
                    do i = 1, n_x
                        if (i < 5 .or. i > n_x-4 .or. j < 5 .or. j > n_x-4 .or. k < 5 .or. k > n_x-4) then
                            row = (k-1)*n_x*n_x + (j-1)*n_x + i
                            x_temp = (i - nx_mid) * h_x
                            y_temp = (j - ny_mid) * h_x
                            z_temp = (k - nz_mid) * h_x
                            solution(row) = total_charge *h_x**3 / (sqrt(x_temp**2 + y_temp**2 + z_temp**2)) ! initial guess for the points within four grid points of the boundary, using the formula for the potential of a point charge, where the total charge is assumed to be concentrated at the center of the grid
                        end if
                    end do
                end do
            end do
            !$omp end parallel do

        ! fill the b(:)
        !$omp parallel do collapse(3) default(none) &
        !$omp private(i,j,k,row,ii,jj,kk,idx) &
        !$omp shared(b,density,n_x,h_x,solution,h2inv,n_inner) &
        !$omp schedule(static)
        do k = 5, n_x-4
            do j = 5, n_x-4
                do i = 5, n_x-4
                    row = (k-1)*n_x*n_x + (j-1)*n_x + i
                    ii = i - 4
                    jj = j - 4
                    kk = k - 4
                    idx = (kk-1)*n_inner*n_inner + (jj-1)*n_inner + ii
                    b(idx) = density(row)  ! because the matrix A in the Poisson equation is negative definite, we need to negate b to get the correct result for Ax = b

                    ! add the contribution from the points within four grid points of the boundary to b(:) using gauss law
                    if (i-4 < 5) then
                        b(idx) = b(idx) &
                        + (-9.0_dp * h2inv *solution(row - 4))
                    end if
                    if (i+4 > n_x - 4) then
                        b(idx) = b(idx) &
                        + (-9.0_dp * h2inv * solution(row + 4))
                    end if
                    if (j-4 < 5) then
                        b(idx) = b(idx) &
                        + (-9.0_dp * h2inv * solution(row - 4*n_x))
                    end if
                    if (j+4 > n_x - 4) then
                        b(idx) = b(idx) &
                        + (-9.0_dp * h2inv * solution(row + 4*n_x))
                    end if
                    if (k-4 < 5) then
                        b(idx) = b(idx) &
                        + (-9.0_dp * h2inv * solution(row - 4*n_x*n_x))
                    end if
                    if (k+4 > n_x - 4) then
                        b(idx) = b(idx) &
                        + (-9.0_dp * h2inv * solution(row + 4*n_x*n_x))
                    end if
                    if (i-3 < 5) then
                        b(idx) = b(idx) &
                        + (128.0_dp * h2inv * solution(row - 3))
                    end if
                    if (i+3 > n_x - 4) then
                        b(idx) = b(idx) &
                        + (128.0_dp * h2inv * solution(row + 3))
                    end if
                    if (j-3 < 5) then
                        b(idx) = b(idx) &
                        + (128.0_dp * h2inv * solution(row - 3*n_x))
                    end if
                    if (j+3 > n_x - 4) then
                        b(idx) = b(idx) &
                        + (128.0_dp * h2inv * solution(row + 3*n_x))
                    end if
                    if (k-3 < 5) then
                        b(idx) = b(idx) &
                        + (128.0_dp * h2inv * solution(row - 3*n_x*n_x))
                    end if
                    if (k+3 > n_x - 4) then
                        b(idx) = b(idx) &
                        + (128.0_dp * h2inv * solution(row + 3*n_x*n_x))
                    end if
                    if (i-2 < 5) then
                        b(idx) = b(idx) &
                        - (1008.0_dp * h2inv * solution(row - 2))
                    end if
                    if (i+2 > n_x - 4) then
                        b(idx) = b(idx) &
                        - (1008.0_dp * h2inv * solution(row + 2))
                    end if
                    if (j-2 < 5) then
                        b(idx) = b(idx) &
                        - (1008.0_dp * h2inv * solution(row - 2*n_x))
                    end if
                    if (j+2 > n_x - 4) then
                        b(idx) = b(idx) &
                        - (1008.0_dp * h2inv * solution(row + 2*n_x))
                    end if
                    if (k-2 < 5) then
                        b(idx) = b(idx) &
                        - (1008.0_dp * h2inv * solution(row - 2*n_x*n_x))
                    end if
                    if (k+2 > n_x - 4) then
                        b(idx) = b(idx) &
                        - (1008.0_dp * h2inv * solution(row + 2*n_x*n_x))
                    end if
                    if (i-1 < 5) then
                        b(idx) = b(idx) &
                        + (8064.0_dp * h2inv * solution(row - 1))
                    end if
                    if (i+1 > n_x - 4) then
                        b(idx) = b(idx) &
                        + (8064.0_dp * h2inv * solution(row + 1))
                    end if
                    if (j-1 < 5) then
                        b(idx) = b(idx) &
                        + (8064.0_dp * h2inv * solution(row - n_x))
                    end if
                    if (j+1 > n_x - 4) then
                        b(idx) = b(idx) &
                        + (8064.0_dp * h2inv * solution(row + n_x))
                    end if
                    if (k-1 < 5) then
                        b(idx) = b(idx) &
                        + (8064.0_dp * h2inv * solution(row - n_x*n_x))
                    end if
                    if (k+1 > n_x - 4) then
                        b(idx) = b(idx) &
                        + (8064.0_dp * h2inv * solution(row + n_x*n_x))
                    end if

                end do
            end do
        end do
        !$omp end parallel do

        ! CG method to solve Ax = b
        x(:) = 0.0_dp ! Initial guess
        call poisson_x(x, n_x-8, h_x, Ax) ! Ax = A*x
        r(:) = b(:) - Ax(:) ! Initial residual
        p(:) = r(:) ! Initial search direction
        rr_old = dot_product(r(:), r(:)) ! Initial residual squared norm
        iter = 1
        do iter = 1, max_iter
            call poisson_x(p, n_x-8, h_x, temp) ! temp = A*p

            n_sizes = size(p)
            denom = 0.0_dp
            !$omp parallel do default(none) private(i) shared(p,temp,n_sizes) &
            !$omp reduction(+:denom) schedule(static)
            do i = 1, n_sizes
                denom = denom + p(i) * temp(i)
            end do
            !$omp end parallel do
            alpha = rr_old / denom

            x(:) = x(:) + alpha * p(:) ! update solution
            
            n_sizes = size(r)
            rr_new = 0.0_dp
            !$omp parallel do default(none) private(i) shared(r,alpha,temp,n_sizes) &
            !$omp reduction(+:rr_new) schedule(static)
            do i = 1, n_sizes
                r(i) = r(i) - alpha * temp(i)
                rr_new = rr_new + r(i) * r(i)
            end do
            !$omp end parallel do

            if (sqrt(rr_new) < tol) then
                print *, "CG converged in ", iter, " iterations."
                exit
            end if

            beta = rr_new / rr_old ! update beta
            p(:) = r(:) + beta * p(:) ! update search direction

            if (mod(iter, 100) == 0) then
                print *, "CG iteration ", iter, ": residual norm = ", sqrt(rr_new)
            end if
            rr_old = rr_new ! update for next iteration
        end do

        if (iter > max_iter) then
            print *, "CG did not converge within the maximum number of iterations."
        end if
        print *, "Conjugate Gradient method finished."
        print *, "Final residual norm: ", sqrt(rr_new)
        deallocate(b, Ax, temp, r, p)

        ! fill the solution(:) with the values from x(:)
        print *, "Filling the solution vector with the values from CG method..."

        n_count = 0
        !$omp parallel do collapse(3) default(none) &
        !$omp private(i,j,k,row,ii,jj,kk,idx) &
        !$omp shared(solution,x,n_x,n_inner) &
        !$omp schedule(static)
        do k = 5, n_x-4
            do j = 5, n_x-4
                do i = 5, n_x-4
                    ! n_count = n_count + 1
                    ! print *, n_count, " / ", n_inner**3
                    row = (k-1)*n_x*n_x + (j-1)*n_x + i

                    ii = i - 4
                    jj = j - 4
                    kk = k - 4
                    idx = (kk-1)*n_inner*n_inner + (jj-1)*n_inner + ii

                    solution(row) = x(idx)

                end do
            end do
        end do
        !$omp end parallel do
        deallocate(x)
        print *, "Solution vector filled with the values from CG method."
    
    end subroutine CG_method_poisson

    subroutine poisson_x(x, n_x_in, h_x, Ax)
        !$ use omp_lib
        implicit none
        real(dp), intent(out) :: x(:)
        integer, intent(in) :: n_x_in
        real(dp), intent(in) :: h_x
        real(dp), intent(out) :: Ax(:)
        integer :: i, j, k, row
        real(dp) :: h2inv
        

        if (size(x) /= n_x_in**3) stop "poisson_x_: size of x must be equal to n_x^3"
        if (size(Ax) /= n_x_in**3) stop "poisson_x_: size of Ax must be equal to n_x^3"

        h2inv = 1.0_dp / (5040.0_dp * h_x * h_x)
        Ax = 0.0_dp
        !$omp parallel do collapse(3) default(none) &
        !$omp private(i,j,k,row) &
        !$omp shared(x,Ax,n_x_in,h2inv) &
        !$omp schedule(static)
        do k = 1, n_x_in
            do j = 1, n_x_in
                do i = 1, n_x_in
                    row = (k-1)*n_x_in*n_x_in + (j-1)*n_x_in + i

                    ! --- diagonal element ---
                    Ax(row) = -43050.0_dp * h2inv * x(row)

                    ! --- x direction neighbors ---
                    if (i+1 <= n_x_in) Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row + 1)
                    if (i-1 >= 1)   Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row - 1)
                    if (i+2 <= n_x_in) Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row + 2)
                    if (i-2 >= 1)   Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row - 2)
                    if (i+3 <= n_x_in) Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row + 3)
                    if (i-3 >= 1)   Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row - 3)
                    if (i+4 <= n_x_in) Ax(row) = Ax(row) -9.0_dp * h2inv * x(row + 4)
                    if (i-4 >= 1)   Ax(row) = Ax(row) -9.0_dp * h2inv * x(row - 4)

                    ! --- y direction neighbors ---
                    if (j+1 <= n_x_in) Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row + n_x_in)
                    if (j-1 >= 1)   Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row - n_x_in)
                    if (j+2 <= n_x_in) Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row + 2*n_x_in)
                    if (j-2 >= 1)   Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row - 2*n_x_in)
                    if (j+3 <= n_x_in) Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row + 3*n_x_in)
                    if (j-3 >= 1)   Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row - 3*n_x_in)
                    if (j+4 <= n_x_in) Ax(row) = Ax(row) -9.0_dp * h2inv * x(row + 4*n_x_in)
                    if (j-4 >= 1)   Ax(row) = Ax(row) -9.0_dp * h2inv * x(row - 4*n_x_in)

                    ! --- z direction neighbors ---
                    if (k+1 <= n_x_in) Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row + n_x_in*n_x_in)
                    if (k-1 >= 1)   Ax(row) = Ax(row) + 8064.0_dp * h2inv * x(row - n_x_in*n_x_in)
                    if (k+2 <= n_x_in) Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row + 2*n_x_in*n_x_in)
                    if (k-2 >= 1)   Ax(row) = Ax(row) -1008.0_dp * h2inv * x(row - 2*n_x_in*n_x_in)
                    if (k+3 <= n_x_in) Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row + 3*n_x_in*n_x_in)
                    if (k-3 >= 1)   Ax(row) = Ax(row) + 128.0_dp * h2inv * x(row - 3*n_x_in*n_x_in)
                    if (k+4 <= n_x_in) Ax(row) = Ax(row) -9.0_dp * h2inv * x(row + 4*n_x_in*n_x_in)
                    if (k-4 >= 1)   Ax(row) = Ax(row) -9.0_dp * h2inv * x(row - 4*n_x_in*n_x_in)    
                end do
            end do
        end do
        !$omp end parallel do
        Ax(:) = -Ax(:) ! because the matrix A in the Poisson equation

    end subroutine poisson_x
    
    



end module CG_method_mod