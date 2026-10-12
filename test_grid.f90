program test_grid
    !---------------------------
    ! Tests for grid_mod with the box determined from (epsilon, gamma)
    !   1. semi-axes: volume conservation, axial symmetry for gamma = 0 and 60 deg
    !   2. box extents (cubic and per-axis) cover the nucleus
    !   3. number of inside points * dV reproduces the nuclear volume
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use nucleus_mod
    use grid_mod
    implicit none

    integer, parameter :: dp = real64
    type(nucleus_property) :: nucleus
    type(grid_type) :: grid
    real(dp) :: s(3), s2(3), vol_grid
    logical :: ok_all
    integer :: icase

    ok_all = .true.

    print *, "--- 1. semi-axes from (eps, gamma) ---"
    s = semi_axes_from_eps_gamma(7.0_dp, 0.3_dp, 25.0_dp * pi / 180.0_dp)
    call check("volume conservation R_x R_y R_z = R^3", abs(product(s) / 7.0_dp**3 - 1.0_dp), 1.0e-14_dp)
    s = semi_axes_from_eps_gamma(7.0_dp, 0.3_dp, 0.0_dp)
    call check("gamma = 0: R_x = R_y               ", abs(s(1) - s(2)), 1.0e-12_dp)
    call check("gamma = 0, eps > 0: prolate along z  ", merge(0.0_dp, 1.0_dp, s(3) > s(1)), 0.0_dp)
    ! (eps, gamma = 60 deg) is an oblate shape with symmetry axis x
    s2 = semi_axes_from_eps_gamma(7.0_dp, 0.3_dp, pi / 3.0_dp)
    call check("gamma = 60 deg: R_y = R_z (oblate)   ", abs(s2(2) - s2(3)), 1.0e-12_dp)
    call check("gamma = 60 deg: R_x is the short axis", merge(0.0_dp, 1.0_dp, s2(1) < s2(2)), 0.0_dp)

    nucleus%Z = 82
    nucleus%N = 126
    nucleus%semi1 = nucleus%R0 * 208.0_dp**(1.0_dp / 3.0_dp)
    nucleus%epsilon2 = 0.3_dp
    nucleus%gamma_deg = 25.0_dp
    call nucleus%calculate_fundamental_properties()

    do icase = 1, 2
        print *, "--- box: cubic = ", icase == 1, " ---"
        grid%h_x = 0.25_dp
        grid%h_y = 0.25_dp
        grid%h_z = 0.25_dp
        grid%n_times = 2
        grid%cubic = (icase == 1)
        call grid%initialize_grid(nucleus)
        allocate(grid%density_index3D(grid%n_x_points, grid%n_y_points, grid%n_z_points))
        allocate(grid%density_index(grid%n_points))
        call grid%inside_outside_nucleus3D(nucleus, grid%density_index3D)
        call grid%inside_outside_nucleus(nucleus, grid%density_index)

        call check("box covers the nucleus (x)           ", &
                   merge(0.0_dp, 1.0_dp, grid%nu_x_max * grid%h_x >= grid%semi(1)), 0.0_dp)
        call check("box covers the nucleus (y)           ", &
                   merge(0.0_dp, 1.0_dp, grid%nu_y_max * grid%h_y >= grid%semi(2)), 0.0_dp)
        call check("box covers the nucleus (z)           ", &
                   merge(0.0_dp, 1.0_dp, grid%nu_z_max * grid%h_z >= grid%semi(3)), 0.0_dp)
        if (grid%cubic) then
            call check("cubic: n_x = n_y = n_z               ", &
                       real(abs(grid%n_x_points - grid%n_y_points) + abs(grid%n_x_points - grid%n_z_points), dp), 0.0_dp)
        else
            call check("per axis: box shorter along x than z ", &
                       merge(0.0_dp, 1.0_dp, grid%n_x_points < grid%n_z_points), 0.0_dp)
        end if
        vol_grid = sum(grid%density_index3D) * grid%dV
        call check("inside points * dV vs volume (rel)   ", abs(vol_grid / nucleus%volume - 1.0_dp), 1.0e-3_dp)
        call check("1D and 3D inside counts agree        ", &
                   abs(sum(grid%density_index) - sum(grid%density_index3D)), 0.0_dp)
        deallocate(grid%density_index3D, grid%density_index)
    end do

    if (ok_all) then
        print *, "ALL TESTS PASSED"
    else
        print *, "SOME TESTS FAILED"
        stop 1
    end if

contains

    subroutine check(name, err, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: err, tol
        if (err <= tol) then
            print '(a, a, es10.2)', "PASS: ", name, err
        else
            print '(a, a, es10.2)', "FAIL: ", name, err
            ok_all = .false.
        end if
    end subroutine check

end program test_grid
