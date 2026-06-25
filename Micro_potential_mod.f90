module Micro_potential_mod
    use iso_fortran_env, only: real64
    use CG_method_mod
    use constant_mod
    use three_D_derivative_mod
    use nucleus_mod, only: nucleus_property
    use grid_mod, only: grid_type
    use micro_constant_mod
    
    implicit none
    private

    integer, parameter :: dp = real64

    type :: potential_type
        real(dp), allocatable :: V_c(:,:,:)
        real(dp), allocatable :: V_1(:,:,:)
        real(dp), allocatable :: V_so(:,:,:,:,:)
        real(dp), allocatable :: nabla_V_1_x(:,:,:)
        real(dp), allocatable :: nabla_V_1_y(:,:,:)
        real(dp), allocatable :: nabla_V_1_z(:,:,:)
        real(dp), allocatable :: psi(:,:,:,:,:)

        contains
            procedure :: calculate_potential_energy
            procedure :: H_on_psi
    end type potential_type
contains

    subroutine calculate_potential_energy(self, V_c, V_1, microscopic_vars, nucleus, grid)
        class(potential_type), intent(in) :: self
        class(microscopic_variables), intent(in) :: microscopic_vars
        class(nucleus_property), intent(in) :: nucleus
        class(grid_type), intent(in) :: grid
        real(dp), intent(out) :: V_c(:,:,:)
        real(dp), intent(out) :: V_1(:,:,:)

        
        
        

        integer :: i, j, k
        integer :: nx, ny, nz

        

        
    end subroutine calculate_potential_energy

    subroutine compute_nabla_V_1(self, grid)
        class(potential_type), intent(inout) :: self
        class(grid_type), intent(in) :: grid

        call d_x(self%V_1, grid%h_x, self%nabla_V_1_x)
        call d_y(self%V_1, grid%h_y, self%nabla_V_1_y)
        call d_z(self%V_1, grid%h_z, self%nabla_V_1_z)

    end subroutine compute_nabla_V_1

    subroutine compute_V_so(self, grid, psi, V_so)
        class(potential_type), intent(inout) :: self
        class(grid_type), intent(in) :: grid
        complex(dp), intent(in) :: psi(:,:,:,:)
        complex(dp), intent(out) :: V_so(:,:,:,:)
        complex(dp), allocatable :: nabla_psi_x(:,:,:,:), nabla_psi_y(:,:,:,:), nabla_psi_z(:,:,:,:)
        integer :: i, j, k, n
        integer :: n_states, n_x, n_y, n_z
        n_states = size(psi, 4)
        n_x = size(psi, 1)
        n_y = size(psi, 2)
        n_z = size(psi, 3)
        allocate(nabla_psi_x(n_x, n_y, n_z, n_states))
        allocate(nabla_psi_y(n_x, n_y, n_z, n_states))
        allocate(nabla_psi_z(n_x, n_y, n_z, n_states))

        !---------------------------
        ! Compute the gradient of psi
        !---------------------------
        call d_x_cp(psi, grid%h_x, nabla_psi_x)
        call d_y_cp(psi, grid%h_y, nabla_psi_y)
        call d_z_cp(psi, grid%h_z, nabla_psi_z)

    end subroutine

    

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

            !---------------------------
            ! compute Coulomb
            !---------------------------


            !---------------------------
            ! compute Yukawa
            !---------------------------

            !---------------------------
            ! compute 

            
            
            

        end subroutine H_on_psi
end module Micro_potential_mod