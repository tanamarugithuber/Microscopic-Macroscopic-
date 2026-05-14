module single_particle_solver_mod
    use iso_fortran_env, only: real64
    use constant_mod, only: pi, e2, hbar_c, m_nucleon, alpha, i_complex
    use nucleus_mod, only: nucleus_property
    use grid_mod, only: grid_type
    use micro_constant_mod
    use CG_method_mod
    use three_D_derivative_mod, only: d_x, d_y, d_z, dd_x, dd_y, dd_z, laplacian
    implicit none
    private

    integer, parameter :: dp = real64

    type :: sp_solver_type
        real(dp), allocatable :: V_1(:,:,:)
        real(dp), allocatable :: V_c(:,:,:)
        real(dp), allocatable :: V_so(:,:,:)
        real(dp), allocatable :: psi(:,:,:,:)
        real(dp), allocatable :: psi_new(:,:,:,:)
        real(dp), allocatable :: energy(:)
        real(dp), allocatable :: energy_new(:)
        integer :: n_states

        contains 
            procedure :: solve_schrodinger_equation
            procedure :: imaginary_time_evolution
            procedure :: compute_V_so
            procedure :: check_orthogonality
            procedure :: Gram_Schmidt_orthogonalization
            

    end type sp_solver_type

    public :: sp_solver_type

    contains
        subroutine solve_schrodinger_equation(this, nucleus, grid, micro_const,number_of_states,energy)
            implicit none
            class(sp_solver_type), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            type(grid_type), intent(in) :: grid
            type(microscopic_variables), intent(in) :: micro_const
            integer, intent(in) :: number_of_states
            real(dp), intent(out) :: energy(number_of_states)
            real(dp) :: psi(grid%n_x_points,grid%n_y_points,grid%n_z_points,number_of_states)
            
            print *, "Solving the Schrödinger equation for single-particle states..."
            
        end subroutine solve_schrodinger_equation

        subroutine imaginary_time_evolution(this, nucleus, grid, micro_const)
            implicit none
            class(sp_solver_type), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            type(grid_type), intent(in) :: grid
            type(microscopic_variables), intent(in) :: micro_const
            print *, "Performing imaginary time evolution to find the ground state wavefunction..."
        end subroutine imaginary_time_evolution


        subroutine compute_V_so(this, nucleus, grid)
            implicit none
            class(sp_solver_type), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            type(grid_type), intent(in) :: grid
            print *, "Computing the spin-orbit potential V_so..."



        end subroutine compute_V_so

        subroutine check_orthogonality(this)
            implicit none
            class(sp_solver_type), intent(inout) :: this
            print *, "Checking the orthogonality of the computed wavefunctions..."
        end subroutine check_orthogonality

        subroutine Gram_Schmidt_orthogonalization(this, wavefunctions)
            implicit none
            class(sp_solver_type), intent(inout) :: this
            real(dp), intent(inout) :: wavefunctions(:,:,:,:)
            print *, "Performing Gram-Schmidt orthogonalization on the wavefunctions..."
        end subroutine Gram_Schmidt_orthogonalization

        subroutine compute_V_1_V_c(this, nucleus, grid, micro_const)
            implicit none
            class(sp_solver_type), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            type(grid_type), intent(in) :: grid
            type(microscopic_variables), intent(in) :: micro_const
            integer :: i, j, k
            
            print *, "Computing the central potential V_1 and Coulomb potential V_c..."
            allocate(this%V_1(grid%n_x_points, grid%n_y_points, grid%n_z_points))
            allocate(this%V_c(grid%n_x_points, grid%n_y_points, grid%n_z_points))
            


        end subroutine compute_V_1_V_c
        

end module single_particle_solver_mod