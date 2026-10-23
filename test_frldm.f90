program test_frldm
    !---------------------------
    ! Tests for frldm_mod
    !   1. B_1, B_3 of a sphere by the double surface integrals against the analytic
    !      Eq. (74) of FRDM2012 (A = 16 and 208)
    !   2. deformed shape (eps2 = 0.22, eps4 = -0.05, eps6 = 0.04): B_3 from the surface
    !      integrals against the grid + Hockney value (B_1 on the grid has an O(h) error
    !      from the smeared surface cells and is only printed)
    !   2b. octupole shape (eps2 = 0.15, eps3 = 0.10, eps4 = -0.04): the same comparison; the surface
    !      integrals use only the x, y reflections, the grid uses no symmetry at all
    !   3. macroscopic mass excess of 208Pb, 132Sn, 40Ca against FRDM2012
    !      (M_th^FL - E_mic^FL: -8.72, -64.98, -37.52 MeV)
    !
    ! build:
    !   gfortran -O2 -fopenmp -I/usr/include constant_mod.f90 nucleus_mod.f90 shape_mod.f90 grid_mod.f90 \
    !       calculate_potential_hockney.f90 frldm_mod.f90 test_frldm.f90 -lfftw3 -o test_frldm.exe
    !---------------------------
    ! 【日本語】frldm_mod のテスト: 球の B_1, B_3 を解析式 (74) と比較, 変形核の B_3 を
    ! 表面積分と格子で比較, 208Pb, 132Sn, 40Ca の巨視的質量超過を FRDM2012 と比較する.
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use nucleus_mod, only: nucleus_property
    use grid_mod, only: grid_type
    use shape_mod, only: epsilon_shape, epsilon_shape_from
    use frldm_mod
    implicit none
    integer, parameter :: dp = real64
    logical :: ok_all

    ok_all = .true.
    print *, "--- 1. sphere, surface integrals vs Eq. (74) ---"
    call test_sphere(8, 8)
    call test_sphere(82, 126)
    print *, "--- 2. deformed shape, surface integrals vs grid ---"
    call test_deformed(94, 146, 0.22_dp, 0.0_dp, -0.05_dp, 0.04_dp)
    print *, "--- 2b. octupole shape (eps3 /= 0), surface integrals vs grid ---"
    call test_deformed(88, 136, 0.15_dp, 0.10_dp, -0.04_dp, 0.0_dp)
    print *, "--- 3. macroscopic mass excess vs FRDM2012 ---"
    call test_mac(82, 126, -8.72_dp)
    call test_mac(50, 82, -64.98_dp)
    call test_mac(20, 20, -37.52_dp)

    if (ok_all) then
        print *, "ALL TESTS PASSED"
    else
        print *, "SOME TESTS FAILED"
        stop 1
    end if

contains

    subroutine check(name, e, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: e, tol
        if (e <= tol) then
            print '(a, a, es10.2, a, es9.1, a)', "PASS: ", name, e, "  (tol ", tol, ")"
        else
            print '(a, a, es10.2, a, es9.1, a)', "FAIL: ", name, e, "  (tol ", tol, ")"
            ok_all = .false.
        end if
    end subroutine check

    subroutine setup(Z, N, e2_, e3_, e4_, e6_, nuc, shape)
        integer, intent(in) :: Z, N
        real(dp), intent(in) :: e2_, e3_, e4_, e6_
        type(nucleus_property), intent(out) :: nuc
        type(epsilon_shape), intent(out) :: shape
        nuc%Z = Z
        nuc%N = N
        nuc%A = Z + N
        nuc%semi1 = nuc%R0 * (Z + N)**(1.0_dp / 3.0_dp)
        nuc%epsilon2 = e2_
        nuc%epsilon3 = e3_
        nuc%epsilon4 = e4_
        nuc%epsilon6 = e6_
        nuc%R_a = nuc%R0 * (Z + N)**(1.0_dp / 3.0_dp)
        shape = epsilon_shape_from(nuc%R_a, e2_, 0.0_dp, e3_, e4_, e6_)
    end subroutine setup

    subroutine test_sphere(Z, N)
        integer, intent(in) :: Z, N
        type(nucleus_property) :: nuc
        type(epsilon_shape) :: shape
        type(frldm_variables) :: fr
        real(dp) :: x0, y0, B10, B30
        call setup(Z, N, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, nuc, shape)
        call fr%calculate_b_1_b_3_surface(nuc, shape)
        x0 = nuc%R_a / a_Yukawa
        y0 = nuc%R_a / a_den
        B10 = 1.0_dp - 3.0_dp / x0**2 + (1.0_dp + x0) * (2.0_dp + 3.0_dp / x0 + 3.0_dp / x0**2) * exp(-2.0_dp * x0)
        B30 = 1.0_dp - 5.0_dp / y0**2 * (1.0_dp - 15.0_dp / (8.0_dp * y0) + 21.0_dp / (8.0_dp * y0**3) &
            - 0.75_dp * (1.0_dp + 4.5_dp / y0 + 7.0_dp / y0**2 + 3.5_dp / y0**3) * exp(-2.0_dp * y0))
        print '(a, i4, a, 2f12.7, a, 2f12.7)', " A = ", Z + N, "  B_1, B_3 = ", fr%B_1, fr%B_3, &
              "   Eq. (74): ", B10, B30
        call check("|B_1 - Eq. (74)|                 ", abs(fr%B_1 - B10), 1.0e-4_dp)
        call check("|B_3 - Eq. (74)|                 ", abs(fr%B_3 - B30), 1.0e-8_dp)
    end subroutine test_sphere

    subroutine test_deformed(Z, N, e2_, e3_, e4_, e6_)
        integer, intent(in) :: Z, N
        real(dp), intent(in) :: e2_, e3_, e4_, e6_
        type(nucleus_property) :: nuc
        type(epsilon_shape) :: shape
        type(frldm_variables) :: fr
        type(grid_type) :: g
        real(dp) :: b1s, b3s
        call setup(Z, N, e2_, e3_, e4_, e6_, nuc, shape)
        call fr%calculate_b_1_b_3_surface(nuc, shape)
        b1s = fr%B_1
        b3s = fr%B_3
        call nuc%calculate_fundamental_properties()
        g%h_x = 0.4_dp
        g%h_y = 0.4_dp
        g%h_z = 0.4_dp
        g%n_times = 3
        call g%initialize_grid(nuc)
        allocate(g%density_index3D(g%n_x_points, g%n_y_points, g%n_z_points))
        call g%inside_outside_nucleus3D(nuc, g%density_index3D)
        call fr%calculate_b_1_b_3(nuc, g)
        print '(a, 2f12.7, a, 2f12.7)', " surface: B_1, B_3 = ", b1s, b3s, "   grid (h = 0.4): ", fr%B_1, fr%B_3
        call check("|B_3 surface - grid|             ", abs(b3s - fr%B_3), 2.0e-3_dp)
    end subroutine test_deformed

    subroutine test_mac(Z, N, ref)
        integer, intent(in) :: Z, N
        real(dp), intent(in) :: ref
        type(nucleus_property) :: nuc
        type(epsilon_shape) :: shape
        type(frldm_variables) :: fr
        real(dp) :: E
        call setup(Z, N, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, nuc, shape)
        fr%verbose = .false.
        call fr%calculate_frldm_energy(nuc, shape, E)
        print '(a, i4, a, i4, a, f10.4, a, f10.4)', " Z = ", Z, "  N = ", N, "  E_mac = ", E, "   FRDM2012: ", ref
        call check("|E_mac - FRDM2012| MeV           ", abs(E - ref), 0.1_dp)
    end subroutine test_mac

end program test_frldm
