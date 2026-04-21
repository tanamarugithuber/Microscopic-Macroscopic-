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
    call calculate_fundamental_properties(nucleus)

    ! Initialize grid
    g_mod%h_x = 0.4_dp
    g_mod%h_y = 0.4_dp
    g_mod%h_z = 0.4_dp
    g_mod%n_times = 3
    call initialize_grid(g_mod, nucleus)
    allocate(g_mod%density_index(g_mod%n_points))
    call inside_outside_nucleus(g_mod, nucleus)

    ! Calculate FRLDM variables
    call calculate_frldm_variables(frldm, nucleus, g_mod)

    print *, "FRLDM Energy: ", frldm%E_frldm

end program main