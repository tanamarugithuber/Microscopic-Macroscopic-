program main
    use iso_fortran_env, only: real64
    use constant_mod
    use nucleus_mod
    use grid_mod
    use CG_method_mod
    use frldm_mod

    implicit none
    integer, parameter :: dp = real64
    type(nucleus_property) :: nucleus
    type(grid_type) :: g_mod
    type(frldm_variables) :: frldm
    


    ! real(dp), allocatable :: A(:,:), b(:), x(:)
    ! Allocate(A(2,2), b(2), x(2))
    ! A(1,1) = 4.0_dp
    ! A(1,2) = 1.0_dp
    ! A(2,1) = 1.0_dp
    ! A(2,2) = 3.0_dp
    ! b(1) = 1.0_dp
    ! b(2) = 2.0_dp
    ! call CG_method(A, b, x)
    ! print *, "Solution x: ", x

    ! Initialize nucleus properties
    nucleus%N = 8
    nucleus%Z = 8
    nucleus%semi1 = nucleus%R0*(nucleus%Z + nucleus%N)**(1.0_dp/3.0_dp)
    call nucleus%calculate_fundamental_properties()

    ! Initialize grid
    g_mod%h_x = 0.4_dp
    g_mod%h_y = 0.4_dp
    g_mod%h_z = 0.4_dp
    g_mod%n_times = 2
    call g_mod%initialize_grid(nucleus)
    allocate(g_mod%density_index(g_mod%n_points))
    call g_mod%inside_outside_nucleus(nucleus, g_mod%density_index)

    ! Calculate FRLDM variables
    call frldm%calculate_exp(g_mod)
    call frldm%calculate_frldm_energy(nucleus, g_mod, frldm%E_frldm)

    print *, "FRLDM Energy: ", frldm%E_frldm

end program main