module calculate_vector_cp_mod
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter :: dp = real64

    public :: overlap_calculation_cp, gram_schmidt_orthogonalization_cp, normalize_wavefunction_cp, cross_product_cp

    contains
        subroutine overlap_calculation_cp(psi1, psi2, overlap, dh)
            implicit none
            complex(dp), intent(in) :: psi1(:,:,:,:), psi2(:,:,:,:)
            complex(dp), intent(out) :: overlap
            integer :: i, j, k, n, l
            integer :: n_states, n_x, n_y, n_z, n_spin
            real(dp), intent(in) :: dh
            ! n_states = size(psi1, 4)
            n_spin = size(psi1, 4)
            n_x = size(psi1, 1)
            n_y = size(psi1, 2)
            n_z = size(psi1, 3)

            !---------------------------
            ! This subroutine calculates the overlap between two wavefunctions psi1 and psi2. 
            ! The actual implementation of the function will depend on the specific model being used for the energy calculation.
            !----------------------------
            overlap = 0.0_dp
            
                do l = 1, n_spin
                    do k = 1, n_x
                        do j = 1, n_y
                            do i = 1, n_z
                                overlap = overlap + conjg(psi1(i,j,k,l)) * psi1(i,j,k,l) * dh**3 + conjg(psi2(i,j,k,l)) * psi2(i,j,k,l) * dh**3
                            end do
                        end do
                    end do
                end do
            
            
        end subroutine overlap_calculation_cp

        subroutine gram_schmidt_orthogonalization_cp(psi, dh)
            implicit none
            complex(dp), intent(inout) :: psi(:,:,:,:,:)
            integer :: n_states
            real(dp), intent(in) :: dh
            integer :: i, j, k, m, n, l
            integer :: n_x, n_y, n_z, n_spin
            real(dp) :: norm
            complex(dp) :: overlap

            n_states = size(psi, 5)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)
            n_spin = size(psi, 4)

            do m = 1, n_states
                ! subtract the projections of the m-th state onto all previous states
                do n = 1, m - 1
                    call overlap_calculation_cp(psi(:,:,:,:,m), psi(:,:,:,:,n), overlap, dh)
                    do l = 1, n_spin
                        do k = 1, n_z
                            do j = 1, n_y
                                do i = 1, n_x
                                    psi(i,j,k,l,m) = psi(i,j,k,l,m) - overlap * psi(i,j,k,l,n)
                                end do
                            end do
                        end do
                    end do
                end do

                ! normalize the m-th state
                call Normalize_wavefunction_cp(psi(:,:,:,:,m), dh)
                

            end do

        end subroutine gram_schmidt_orthogonalization_cp

        subroutine normalize_wavefunction_cp(psi, dh)
            implicit none
            complex(dp), intent(inout) :: psi(:,:,:,:)
            ! integer :: n_states
            real(dp), intent(in) :: dh
            integer :: i, j, k, n, l
            integer :: n_x, n_y, n_z, n_spin
            real(dp) :: norm
            complex(dp) :: overlap

            ! n_states = size(psi, 5)
            n_spin = size(psi, 4)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)


                

                call overlap_calculation_cp(psi(:,:,:,:), psi(:,:,:,:), overlap, dh)
                norm = abs(overlap)
                norm = sqrt(norm)
            do l = 1, n_spin
                do k = 1, n_x
                    do j = 1, n_y
                        do i = 1, n_z
                            psi(i,j,k,l) = psi(i,j,k,l) / norm
                        end do
                    end do
                end do
            end do


        end subroutine normalize_wavefunction_cp

        subroutine cross_product_cp(psi1_x, psi1_y, psi1_z, psi2_x, psi2_y, psi2_z, cross_x, cross_y, cross_z)
            implicit none
            complex(dp), intent(in) :: psi1_x(:,:,:,:), psi1_y(:,:,:,:), psi1_z(:,:,:,:)
            complex(dp), intent(in) :: psi2_x(:,:,:,:), psi2_y(:,:,:,:), psi2_z(:,:,:,:)
            complex(dp), intent(out) :: cross_x(:,:,:,:), cross_y(:,:,:,:), cross_z(:,:,:,:)
            integer :: i, j, k, l
            integer :: n_x, n_y, n_z, n_spin
            n_x = size(psi1_x,1)
            n_y = size(psi1_x,2)
            n_z = size(psi1_x,3)
            n_spin = size(psi1_x,4)

            do l = 1, n_spin
                do k = 1, n_z
                    do j = 1, n_y
                        do i = 1, n_x
                            cross_x(i,j,k,l) = psi1_y(i,j,k,l) * psi2_z(i,j,k,l) - psi1_z(i,j,k,l) * psi2_y(i,j,k,l)
                            cross_y(i,j,k,l) = psi1_z(i,j,k,l) * psi2_x(i,j,k,l) - psi1_x(i,j,k,l) * psi2_z(i,j,k,l)
                            cross_z(i,j,k,l) = psi1_x(i,j,k,l) * psi2_y(i,j,k,l) - psi1_y(i,j,k,l) * psi2_x(i,j,k,l)
                        end do
                    end do
                end do
            end do
        end subroutine cross_product_cp
    
        subroutine overlap_of_2Dmatrices_and_spinor_cp(vec_x,vec_y,vec_z,matrix_x,matrix_y,matrix_z,overlap)
            implicit none
            ! spinor components of the wavefunction
            complex(dp), intent(in) :: vec_x(:,:,:,:), vec_y(:,:,:,:), vec_z(:,:,:,:)
             
            complex(dp), intent(in) :: matrix_x(:,:), matrix_y(:,:), matrix_z(:,:)
            complex(dp), intent(out) :: overlap(:,:,:,:)
            integer :: i, j, k, m, n
            integer :: n_x, n_y, n_z, n_spinor
            n_spinor = size(vec_x,4)
            n_x = size(vec_x,1)
            n_y = size(vec_x,2)
            n_z = size(vec_x,3)

            if (size(matrix_x,1) /= n_x .or. size(matrix_x,2) /= n_y .or. size(matrix_y,1)&
             /= n_x .or. size(matrix_y,2) /= n_y .or. size(matrix_z,1) /= n_x .or. size(matrix_z,2) /= n_y) then
                print *, "Error: The dimensions of the matrices do not match the dimensions of the spinor components."
                stop
            end if
            overlap = 0.0_dp
            do n = 1, n_spinor
                do m = 1, n_spinor
                    do k = 1, n_z
                        do j = 1, n_y
                            do i = 1, n_x
                                overlap(i,j,k,m) = overlap(i,j,k,m) + vec_x(i,j,k,n) * matrix_x(m,n) &
                                + vec_y(i,j,k,n) * matrix_y(m,n) + vec_z(i,j,k,n) * matrix_z(m,n)
                            end do
                        end do
                    end do
                end do
            end do

        end subroutine overlap_of_2Dmatrices_and_spinor_cp

end module calculate_vector_cp_mod