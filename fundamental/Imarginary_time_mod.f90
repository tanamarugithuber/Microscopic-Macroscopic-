module Imarginary_time_mod
    use iso_fortran_env, only: real64

    implicit none

    integer, parameter :: dp = real64

    !---------------------------
    ! Variables for Imaginary Time calculations
    !---------------------------


    contains
        subroutine initialize_imaginary_time_variables(psi,dh)
            implicit none
            complex(dp), intent(inout) :: psi(:,:,:,:)
            real(dp), intent(in) :: dh
            integer :: n_states, n_x, n_y, n_z, n_middle
            integer :: i, j, k, n
            real(dp) :: norm
            n_states = size(psi, 4)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)

            n_middle = (n_x + 1) / 2  ! Assuming n_x is odd, this gives the middle index
            
            ! ---------------------------
            ! Initialize the wavefunction psi for imaginary time evolution.
            ! Use Gaussian wave packets or any other suitable initial guess for the wavefunction.
            ! ---------------------------
            do n = 1, n_states
                do k = 1, n_x
                    do j = 1, n_y
                        do i = 1, n_z
                            psi(i,j,k,n) = exp(-((real(i - n_middle, dp)**2 + real(j - n_middle, dp)**2 + real(k - n_middle, dp)**2) / (2.0_dp * 5.0_dp**2)))
                        end do
                    end do
                end do
            end do


            !---------------------------
            ! Normalize the initial wavefunction
            !---------------------------
            norm = 0.0_dp
            do n = 1, n_states
                do k = 1, n_x
                    do j = 1, n_y
                        do i = 1, n_z
                            norm = norm + real(conjg(psi(i,j,k,n))*psi(i,j,k,n))*dh**3
                        end do
                    end do
                end do
                norm = sqrt(norm)
                do k = 1, n_x
                    do j = 1, n_y
                        do i = 1, n_z
                            psi(i,j,k,n) = psi(i,j,k,n) / norm
                        end do
                    end do
                end do
            end do
            
        end subroutine initialize_imaginary_time_variables

        subroutine update_imaginary_time_variables(psi, dh)
            implicit none
            complex(dp), intent(inout) :: psi(:,:,:,:)
            real(dp), intent(in) :: dh
            integer :: i, j, k, n
            integer :: n_states, n_x, n_y, n_z
            real(dp) :: norm
            n_states = size(psi, 4)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)

            
            


            
        end subroutine update_imaginary_time_variables



        subroutine Pretreatment(psi, dh, DPI_psi)
            implicit none
            complex(dp), intent(inout) :: psi(:,:,:,:)
            real(dp), intent(in) :: dh
            complex(dp), intent(out) :: DPI_psi(:,:,:,:)
            integer :: n_states, n_x, n_y, n_z, n_middle
            integer :: i, j, k, n
            n_states = size(psi, 4)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)
            n_middle = (n_x + 1) / 2
            ! This subroutine can be used for any necessary pretreatment before starting the imaginary time evolution, such as initializing wavefunctions, setting up the grid, etc.



                 
            print *, "Pretreatment for imaginary time evolution is done."

        end subroutine Pretreatment



        ! subroutine three_point_folding(psi,d)
        !     implicit none
        !     complex(dp), intent(in) :: psi(:,:,:,:)
        !     complex(dp), intent(out) :: d(:,:,:,:)
        !     integer :: i, j, k, n
        !     integer :: n_states, n_x, n_y, n_z
        !     n_states = size(psi, 4)
        !     n_x = size(psi, 1)
        !     n_y = size(psi, 2)
        !     n_z = size(psi, 3)

        !     !---------------------------
        !     ! This subroutine calculates the three-point folding of the wavefunction psi to get the derivative d. 
        !     !----------------------------
           
        ! end subroutine three_point_folding

        

        

end module Imarginary_time_mod