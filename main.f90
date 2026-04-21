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

    ! Initialize nucleus properties
    nucleus%A = 8
    nucleus%Z = 8
    nucleus%semi1 = nucleus%R0*nucleus%A**(1.0_dp/3.0_dp)
    call nucleus%calculate_fundamental_properties()

    ! Initialize grid
    g_mod%h_x = 0.4_dp
    g_mod%h_y = 0.4_dp
    g_mod%h_z = 0.4_dp
    g_mod%n_times = 3
    call g_mod%initialize_grid(nucleus)
    allocate(g_mod%density_index(g_mod%n_points))
    call g_mod%inside_outside_nucleus(nucleus, g_mod%density_index)

    ! Calculate FRLDM variables
    call frldm%calculate_frldm_energy(nucleus, g_mod, frldm%E_frldm)

    print *, "FRLDM Energy: ", frldm%E_frldm

end program main