module Micro_potential_mod
    use iso_fortran_env, only: real64
    use CG_method_mod
    use constant_mod
    use nucleus_mod, only: nucleus_property
    use grid_mod, only: grid_type
    use micro_constant_mod
    
    implicit none
    private

    integer, parameter :: dp = real64

    type :: potential_type
        real(dp), allocatable :: V(:,:,:)


        contains
            procedure :: calculate_potential_energy
            procedure :: H_on_psi
    end type potential_type
contains

    subroutine calculate_potential_energy(self, potential, microscopic_vars, nucleus, grid)
        class(potential_type), intent(inout) :: self
        real(dp), intent(out) :: potential(:,:,:)
        

        integer :: i, j, k
        integer :: nx, ny, nz

        nx = size(x)
        ny = size(y)
        nz = size(z)

        
    end subroutine calculate_potential_energy

    subroutine H_on_psi(self, psi, H_psi)
        class(potential_type), intent(inout) :: self
        complex(dp), intent(inout) :: psi(:,:,:,:)
            complex(dp), intent(out) :: H_psi(:,:,:,:)
            integer :: i, j, k, n
            integer :: n_states, n_x, n_y, n_z
            n_states = size(psi, 4)
            n_x = size(psi, 1)
            n_y = size(psi, 2)
            n_z = size(psi, 3)

            !---------------------------
            ! This subroutine calculates H*psi, where H is the Hamiltonian operator. 
            !----------------------------

            
            
            

        end subroutine H_on_psi
end module Micro_potential_mod