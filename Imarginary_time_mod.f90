module Imarginary_time_mod
    use iso_fortran_env, only: real64
    use constant_mod
    use calculate_vector_cp_mod
    use three_D_derivative_cp_mod
    use micro_constant_mod
    use, intrinsic :: ieee_arithmetic

    implicit none

    integer, parameter :: dp = real64

    !---------------------------
    ! Variables for Imaginary Time calculations
    !---------------------------


    contains
        subroutine initialize_imaginary_time_variables(psi,dh)
            implicit none
            complex(dp), intent(inout) :: psi(:,:,:,:,:)
            real(dp), intent(in) :: dh
            integer :: n_states, n_x, n_y, n_z, n_middle, n_spin
            integer :: i, j, k, n, l
            real(dp) :: norm
            n_states = size(psi, 5)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)
            n_spin = size(psi, 4)


            n_middle = (n_x + 1) / 2  ! Assuming n_x is odd, this gives the middle index
            
            ! ---------------------------
            ! Initialize the wavefunction psi for imaginary time evolution.
            ! Use Gaussian wave packets or any other suitable initial guess for the wavefunction.
            ! ---------------------------
            do n = 1, n_states
                do l = 1, n_spin
                    do k = 1, n_z
                        do j = 1, n_y
                            do i = 1, n_x
                                psi(i,j,k,l,n) = exp(-((real(i - n_middle, dp)**2 + real(j - n_middle, dp)**2 + real(k - n_middle, dp)**2) / (2.0_dp * 5.0_dp**2)))
                            end do
                        end do
                    end do
                end do
            end do


            !---------------------------
            ! Normalize the initial wavefunction
            !---------------------------
            do n = 1, n_states
                call normalize_wavefunction_cp(psi(:,:,:,:,n), dh)
            end do
            
        end subroutine initialize_imaginary_time_variables

        subroutine imaginary_time_evolution(micro_vars, psi, dh, V_ext, dis)
            implicit none
            type(microscopic_variables), intent(in) :: micro_vars
            complex(dp), intent(inout) :: psi(:,:,:,:,:)
            complex(dp), allocatable :: H_psi(:,:,:,:,:)
            
            character(len=*), intent(in) :: dis
            real(dp), intent(in) :: V_ext(:,:,:)
            real(dp), intent(in) :: dh
            integer :: i, j, k, n, l, m, iter, max_iter
            integer :: n_states, n_x, n_y, n_z, n_spin
            real(dp) :: norm
            complex(dp) :: overlap
            real(dp) :: xi_0
            real(dp), allocatable :: Eigen_energy(:)
            real(dp), allocatable :: Eigen_energy_old(:)
            complex(dp), allocatable :: DPI_psi(:,:,:,:,:)
            complex(dp), allocatable :: b(:,:,:,:,:)
            complex(dp), allocatable :: psi_old(:,:,:,:,:)
            n_states = size(psi, 5)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)
            n_spin = size(psi, 4)

           
            allocate(H_psi(n_x, n_y, n_z, n_spin, n_states))
            allocate(Eigen_energy(n_states))
            allocate(DPI_psi(n_x, n_y, n_z, n_spin, n_states))
            allocate(b(n_x, n_y, n_z, n_spin, n_states))
            allocate(psi_old(n_x, n_y, n_z, n_spin, n_states))
            psi_old(:,:,:,:,:) = 0.0_dp
            max_iter = 100000  ! Set the maximum number of iterations for imaginary time evolution

            do iter = 1, max_iter
                b(:,:,:,:,:) = 0.0_dp
                DPI_psi(:,:,:,:,:) = 0.0_dp
                H_psi(:,:,:,:,:) = 0.0_dp
                ! check convergence by calculating the overlap between the new and old wavefunctions
                


                ! Calculate the Hamiltonian acting on psi
                do n = 1, n_states
                    call Hamiltonian_on_psi(micro_vars, psi(:,:,:,:,n), dh, V_ext, H_psi(:,:,:,:,n), dis)
                    call overlap_calculation_cp(psi(:,:,:,:,n), H_psi(:,:,:,:,n), overlap, dh)
                    Eigen_energy(n) = real(overlap, dp)
                end do

                if (mod(iter - 1, 100) == 0) then
                          ! Store the current wavefunction for comparison
                    do n = 1, n_states
                        call overlap_calculation_cp(psi(:,:,:,:,n), psi_old(:,:,:,:,n), overlap, dh)
                        norm = abs(overlap)
                        if (norm < 1.0e-6_dp .and. all(abs(Eigen_energy - Eigen_energy_old) < 1.0e-6_dp)) then
                            print *, "Converged after ", iter, " iterations."
                            exit
                        end if
                    end do

                end if
                

                ! pretereatment 
                
                do n = 1, n_states
                    b(:,:,:,:,n) = H_psi(:,:,:,:,n) - Eigen_energy(n) * psi(:,:,:,:,n)
                end do

                do n = 1, n_states
                    do l = 1, n_spin
                        call CG_method_3D(b(:,:,:,l,n), DPI_psi(:,:,:,l,n), n_x, dh, Eigen_energy(n))
                    end do
                end do

                ! Imaginary time evolution step
                do n = 1, n_states
                    psi(:,:,:,:,n) = psi(:,:,:,:,n) + DPI_psi(:,:,:,:,n)
                    call normalize_wavefunction_cp(psi(:,:,:,:,n), dh)
                end do

                if (mod(iter, 100) == 0) then
                    
                    Eigen_energy_old(:) = Eigen_energy(:)
                    
                    call gram_schmidt_orthogonalization_cp(psi, dh)
                    ! check convergence by calculating the overlap between the new and old wavefunctions
                    do n = 1, n_states
                        call overlap_calculation_cp(psi(:,:,:,:,n), psi_old(:,:,:,:,n), overlap, dh)
                        norm = abs(overlap)
                        if (abs(norm - 1.0_dp) < 1.0e-6_dp .and. all(abs(Eigen_energy - Eigen_energy_old) < 1.0e-6_dp)) then
                            print *, "Converged after ", iter, " iterations."
                            exit
                        end if
                    end do 
                    psi_old(:,:,:,:,:) = psi(:,:,:,:,:)
                end if


                

            end do
            
           

            
                    
          
            
        end subroutine imaginary_time_evolution

        subroutine V_so(micro_vars, dis, psi, dh, V_so_psi)
            implicit none
            type(microscopic_variables), intent(in) :: micro_vars
            complex(dp), intent(in) :: psi(:,:,:,:)
            character(len=*), intent(in) :: dis
            real(dp), intent(in) :: dh
            integer :: i, j, k, n, l, m
            integer :: n_spin, n_x, n_y, n_z

            ! V_so
            complex(dp) :: paulix(2,2), pauliy(2,2), pauliz(2,2)
            complex(dp), allocatable :: px_psi(:,:,:,:), py_psi(:,:,:,:), pz_psi(:,:,:,:)
            complex(dp), allocatable :: cross_x(:,:,:,:), cross_y(:,:,:,:), cross_z(:,:,:,:)
            complex(dp), intent(out) :: V_so_psi(:,:,:,:)

            n_spin = size(psi, 4)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)

            ! allocate the derivative arrays
            allocate(px_psi(n_x, n_y, n_z, n_spin))
            allocate(py_psi(n_x, n_y, n_z, n_spin))
            allocate(pz_psi(n_x, n_y, n_z, n_spin))
            allocate(cross_x(n_x, n_y, n_z, n_spin))
            allocate(cross_y(n_x, n_y, n_z, n_spin))
            allocate(cross_z(n_x, n_y, n_z, n_spin))


            ! set up Pauli matrices
            paulix = reshape([(0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)], [2, 2])
            pauliy = reshape([(0.0_dp, 0.0_dp), (0.0_dp, -1.0_dp), (0.0_dp, 1.0_dp), (0.0_dp, 0.0_dp)], [2, 2])
            pauliz = reshape([(1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (-1.0_dp, 0.0_dp)], [2, 2])

            ! calculate the derivatives of psi
            do l = 1, n_spin
                call d_x_cp(psi(:,:,:,l), dh, px_psi(:,:,:,l))
                call d_y_cp(psi(:,:,:,l), dh, py_psi(:,:,:,l))
                call d_z_cp(psi(:,:,:,l), dh, pz_psi(:,:,:,l))
            end do

            call cross_product_cp(px_psi, py_psi, pz_psi, px_psi, py_psi, pz_psi, cross_x, cross_y, cross_z)
            call overlap_of_2Dmatrices_and_spinor_cp(cross_x, cross_y, cross_z, paulix, pauliy, pauliz, V_so_psi)

            deallocate(px_psi)
            deallocate(py_psi)
            deallocate(pz_psi)
            deallocate(cross_x)
            deallocate(cross_y)
            deallocate(cross_z)
            if (dis == "proton") then
                V_so_psi = V_so_psi * micro_vars%lambda_p* i_complex * hbar_c**2 / (2.0_dp * m_nucleon )**2

            else if (dis == "neutron") then
                V_so_psi = V_so_psi * micro_vars%lambda_n* i_complex * hbar_c**2 / (2.0_dp * m_nucleon )**2
            else
                print *, "Error: Unknown particle type in Hamiltonian_on_psi."
            end if
 

        end subroutine V_so



 
        subroutine Hamiltonian_on_psi(micro_vars, psi, dh, V_ext, H_psi, dis)
            implicit none
            type(microscopic_variables), intent(in) :: micro_vars
            character(len=*), intent(in) :: dis
            complex(dp), intent(in) :: psi(:,:,:,:)
            real(dp), intent(in) :: dh
            real(dp), intent(in) :: V_ext(:,:,:)
            complex(dp), intent(out) :: H_psi(:,:,:,:)
            complex(dp), allocatable :: laplacian_psi(:,:,:,:)
            complex(dp), allocatable :: V_so_psi(:,:,:,:)
            integer :: n_spin, n_x, n_y, n_z
            integer :: i, j, k, n, l

            n_spin = size(psi, 4)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)

            allocate(laplacian_psi(n_x, n_y, n_z, n_spin))
            allocate(V_so_psi(n_x, n_y, n_z, n_spin))

            call V_so(micro_vars, dis, psi, dh, V_so_psi)
            do l = 1, n_spin
                call laplacian_cp(psi(:,:,:,l), dh, laplacian_psi(:,:,:,l))
            end do

            do l = 1, n_spin
                H_psi(:,:,:,l) = - (hbar_c**2 / (2.0_dp * m_nucleon)) * laplacian_psi(:,:,:,l) + V_ext(:,:,:) * psi(:,:,:,l) + V_so_psi(:,:,:,l)
            end do
         
        end subroutine Hamiltonian_on_psi   
        
        subroutine CG_method_3D(b, x, n_x, h_x, E_n)
            !$ use omp_lib
            implicit none
            complex(dp), intent(in) :: b(:,:,:)
            complex(dp), intent(out) :: x(:,:,:)
            integer, intent(in) :: n_x
            real(dp), intent(in) :: h_x
            complex(dp), allocatable :: Ax(:,:,:)
            real(dp), intent(in) :: E_n
            real(dp) :: tol ! tolerance for convergence
            integer, parameter :: max_iter = 10000
            integer  :: iter
            complex(dp), allocatable :: temp(:,:,:), r(:,:,:), p(:,:,:)
            complex(dp) :: alpha, beta, rr_old, rr_new, denom
            integer :: n_size, i, j, k,l

            integer :: n, n_count
            print *, "Starting Conjugate Gradient method..."
            n = size(b)

            tol = 1.0e-20_dp
            allocate(temp(n_x,n_x,n_x), r(n_x,n_x,n_x), p(n_x,n_x,n_x), Ax(n_x,n_x,n_x))

            x(:,:,:) = 0.0_dp ! Initial guess

            call Ax_3D(x, n_x, h_x, Ax, E_n) ! Calculate Ax for the initial guess
            r(:,:,:) = b(:,:,:) - Ax(:,:,:) ! Initial residual
            p(:,:,:) = r(:,:,:) ! Initial search direction
            ! rr_old = dot_product(r(:,:,:), r(:,:,:)) ! Initial residual squared norm
            rr_old = 0.0_dp
            do i = 1, n_x
                do j = 1, n_x
                    do k = 1, n_x
                        rr_old = rr_old + r(i,j,k) * r(i,j,k)
                    end do
                end do
            end do

            iter = 1
            do iter = 1, max_iter
                call Ax_3D(p, n_x, h_x, temp, E_n) ! temp = A*p

                n_size = size(p)
                denom = 0.0_dp
                !$omp parallel do default(none) private(n) shared(p,temp,n_x) &
                !$omp reduction(+:denom) schedule(static)
                do k = 1, n_x
                    do j = 1, n_x
                        do i = 1, n_x
                            denom = denom + p(i,j,k) * temp(i,j,k)
                        end do
                    end do
                end do
                !$omp end parallel do
                if (.not. ieee_is_finite(abs(denom))) then
                    print *, "ERROR: denom is not finite: ", denom
                    stop
                end if

                if (abs(denom) < 1.0e-300_dp) then
                    print *, "ERROR: denom is too small: ", denom
                    stop
                end if
                alpha = rr_old / denom

                x(:,:,:) = x(:,:,:) + alpha * p(:,:,:) ! update solution
                
                n_size = size(r)
                rr_new = 0.0_dp
                !$omp parallel do default(none) private(n) shared(r,alpha,temp,n_x) &
                !$omp reduction(+:rr_new) schedule(static)
                do k = 1, n_x
                    do j = 1, n_x
                        do i = 1, n_x
                            r(i,j,k) = r(i,j,k) - alpha * temp(i,j,k)
                            rr_new = rr_new + r(i,j,k) * r(i,j,k)
                        end do
                    end do
                end do
                !$omp end parallel do

                if (sqrt(abs(rr_new)) < tol) then
                    print *, "CG converged in ", iter, " iterations."
                    exit
                end if

                beta = rr_new / rr_old ! update beta
                p(:,:,:) = r(:,:,:) + beta * p(:,:,:) ! update search direction

                if (mod(iter, 100) == 0) then
                    print *, "CG iteration ", iter, ": residual norm = ", sqrt(abs(rr_new))
                end if
                rr_old = rr_new ! update for next iteration
            end do

            if (iter > max_iter) then
                print *, "CG did not converge within the maximum number of iterations."
            end if
            print *, "Conjugate Gradient method finished."
            print *, "Final residual norm: ", sqrt(abs(rr_new))
            deallocate(temp, r, p, Ax)

            print *, "Helmholtz equation solved."
        end subroutine CG_method_3D

        subroutine Ax_3D(x, n_x, h_x, Ax, E_n)
            !$ use omp_lib
            implicit none
            complex(dp), intent(in) :: x(:,:,:)
            integer, intent(in) :: n_x
            real(dp), intent(in) :: h_x
            complex(dp), intent(out) :: Ax(:,:,:)
            real(dp), intent(in) :: E_n
            
            integer :: i, j, k, row
            real(dp) :: h2inv
            ! print *, "Calculating Ax for the Helmholtz equation..."
            if (size(x) /= n_x**3) stop "Ax_3D: size of x must be equal to n_x^3"

            h2inv = 1.0_dp / (5040.0_dp * h_x * h_x)
            Ax = 0.0_dp

            ! loop over the 3Dgrid points
            !$omp parallel do collapse(3) default(none) &
            !$omp private(i,j,k,row) &
            !$omp shared(x,Ax,n_x,h2inv) &
            !$omp schedule(static)
            do k = 1, n_x
                do j = 1, n_x
                    do i = 1, n_x
                        ! current grid point's 1D index (row number)
                        row = (k-1)*n_x*n_x + (j-1)*n_x + i

                        ! --- diagonal element ---
                        Ax(i,j,k) = -43050.0_dp * h2inv * x(i,j,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) ) + E_n * x(i,j,k)

                        ! --- x direction neighbors ---
                        if (i+1 <= n_x) Ax(i,j,k) = Ax(i,j,k) + 8064.0_dp * h2inv * x(i+1,j,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (i-1 >= 1)   Ax(i,j,k) = Ax(i,j,k) + 8064.0_dp * h2inv * x(i-1,j,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (i+2 <= n_x) Ax(i,j,k) = Ax(i,j,k) -1008.0_dp * h2inv * x(i+2,j,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (i-2 >= 1)   Ax(i,j,k) = Ax(i,j,k) -1008.0_dp * h2inv * x(i-2,j,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (i+3 <= n_x) Ax(i,j,k) = Ax(i,j,k) + 128.0_dp * h2inv * x(i+3,j,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (i-3 >= 1)   Ax(i,j,k) = Ax(i,j,k) + 128.0_dp * h2inv * x(i-3,j,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (i+4 <= n_x) Ax(i,j,k) = Ax(i,j,k) -9.0_dp * h2inv * x(i+4,j,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (i-4 >= 1)   Ax(i,j,k) = Ax(i,j,k) -9.0_dp * h2inv * x(i-4,j,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )

                        ! --- y direction neighbors ---
                        if (j+1 <= n_x) Ax(i,j,k) = Ax(i,j,k) + 8064.0_dp * h2inv * x(i,j+1,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (j-1 >= 1)   Ax(i,j,k) = Ax(i,j,k) + 8064.0_dp * h2inv * x(i,j-1,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (j+2 <= n_x) Ax(i,j,k) = Ax(i,j,k) -1008.0_dp * h2inv * x(i,j+2,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (j-2 >= 1)   Ax(i,j,k) = Ax(i,j,k) -1008.0_dp * h2inv * x(i,j-2,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (j+3 <= n_x) Ax(i,j,k) = Ax(i,j,k) + 128.0_dp * h2inv * x(i,j+3,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (j-3 >= 1)   Ax(i,j,k) = Ax(i,j,k) + 128.0_dp * h2inv * x(i,j-3,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (j+4 <= n_x) Ax(i,j,k) = Ax(i,j,k) -9.0_dp * h2inv * x(i,j+4,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (j-4 >= 1)   Ax(i,j,k) = Ax(i,j,k) -9.0_dp * h2inv * x(i,j-4,k) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )

                        ! --- z direction neighbors ---
                        if (k+1 <= n_x) Ax(i,j,k) = Ax(i,j,k) + 8064.0_dp * h2inv * x(i,j,k+1) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (k-1 >= 1)   Ax(i,j,k) = Ax(i,j,k) + 8064.0_dp * h2inv * x(i,j,k-1) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (k+2 <= n_x) Ax(i,j,k) = Ax(i,j,k) -1008.0_dp * h2inv * x(i,j,k+2) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (k-2 >= 1)   Ax(i,j,k) = Ax(i,j,k) -1008.0_dp * h2inv * x(i,j,k-2) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (k+3 <= n_x) Ax(i,j,k) = Ax(i,j,k) + 128.0_dp * h2inv * x(i,j,k+3) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (k-3 >= 1)   Ax(i,j,k) = Ax(i,j,k) + 128.0_dp * h2inv * x(i,j,k-3) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (k+4 <= n_x) Ax(i,j,k) = Ax(i,j,k) -9.0_dp * h2inv * x(i,j,k+4) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )
                        if (k-4 >= 1)   Ax(i,j,k) = Ax(i,j,k) -9.0_dp * h2inv * x(i,j,k-4) * ( - hbar_c**2 / (2.0_dp * m_nucleon) )

                    end do
                end do
            end do
            !$omp end parallel do

            Ax = Ax * xi_0
            
        end subroutine Ax_3D

        

end module Imarginary_time_mod