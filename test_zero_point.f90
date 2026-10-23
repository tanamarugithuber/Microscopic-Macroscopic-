program test_zero_point
    !---------------------------
    ! Tests for zero_point_mod
    !   1. stiffness_fit returns 2a for any parabola a x^2 + b x + c, and the least-squares
    !      curvature (not the local one) when a quartic term is added
    !   2. inertias: B_eps2^irr(0) = 2/15 M0 R0^2, B_gamma^irr -> 2/15 eps2^2 M0 R0^2 for small eps2
    !   3. 208Pb: E_pot(eps2 = 0, +-0.05, +-0.10, +-0.15) of main.exe gives E_0,eps2 = 0.687 MeV
    !      (E_zp = 2 E_0,eps2 = 1.37 MeV; FRDM2012 table E_mic^FL - E_s+p = 1.34 MeV)
    !
    ! build:
    !   gfortran -O2 constant_mod.f90 nucleus_mod.f90 micro_constant_mod.f90 zero_point_mod.f90 \
    !       test_zero_point.f90 -o test_zero_point.exe
    !---------------------------
    ! 【日本語】zero_point_mod のテスト: 最小二乗フィットのかたさ C, 慣性 B^irr の極限,
    ! 208Pb の零点エネルギーが FRDM2012 の表から読み取れる値と合うこと.
    use iso_fortran_env, only: real64
    use nucleus_mod, only: nucleus_property
    use zero_point_mod
    implicit none
    integer, parameter :: dp = real64
    logical :: ok_all
    type(nucleus_property) :: nuc
    type(zero_eps2) :: zp
    real(dp) :: E(7), x, e0, c_ls
    integer :: k
    real(dp), parameter :: h = 0.05_dp
    ! E_pot of 208Pb at eps2 = -0.15 .. 0.15 (main.exe, develop1 05e03e3)
    real(dp), parameter :: E_pb(7) = [-14.0765_dp, -17.2285_dp, -21.1505_dp, -22.6370_dp, &
                                      -21.0528_dp, -16.0063_dp, -8.9127_dp]

    ok_all = .true.
    print *, "--- 1. least-squares stiffness ---"
    do k = -3, 3
        x = k * h
        E(k + 4) = 300.0_dp * x**2 - 7.0_dp * x + 2.0_dp
    end do
    call check("parabola: C - 2a                 ", abs(stiffness_fit(E, h) - 600.0_dp), 1.0e-9_dp)
    do k = -3, 3
        x = k * h
        E(k + 4) = x**4
    end do
    ! sum k^6 = 2 (1 + 64 + 729) = 1588, sum k^4 = 196: a = (1588 - 4 * 196) h^2 / 84
    c_ls = 2.0_dp * (1588.0_dp - 4.0_dp * 196.0_dp) * h**2 / 84.0_dp
    call check("x^4: C - least-squares value     ", abs(stiffness_fit(E, h) - c_ls), 1.0e-12_dp)

    print *, "--- 2. inertias ---"
    nuc%Z = 82
    nuc%N = 126
    nuc%semi1 = nuc%R0 * 208.0_dp**(1.0_dp / 3.0_dp)
    call nuc%calculate_fundamental_properties()
    call check("B_eps2(0) / (2/15 M0 R0^2) - 1   ", &
               abs(B_eps2_irr(nuc, 0.0_dp) / (2.0_dp / 15.0_dp * nuc%mass * nuc%R_a**2) - 1.0_dp), 1.0e-14_dp)
    call check("B_gamma(0)                       ", abs(B_gamma_irr(nuc, 0.0_dp)), 1.0e-14_dp)
    call check("B_gamma(1e-4) / (2/15 e^2 M0 R0^2)", &
               abs(B_gamma_irr(nuc, 1.0e-4_dp) / (2.0_dp / 15.0_dp * 1.0e-8_dp * nuc%mass * nuc%R_a**2) - 1.0_dp), &
               1.0e-3_dp)

    print *, "--- 3. 208Pb ---"
    call zp%cal_energy(nuc, 0.0_dp, h, E_pb, e0)
    print '(a, f10.3, a, f8.4, a)', " C = ", stiffness_fit(E_pb, h), " MeV,  E_0,eps2 = ", e0, " MeV"
    call check("E_0,eps2 - 0.687 MeV             ", abs(e0 - 0.687_dp), 1.0e-3_dp)
    call check("2 E_0,eps2 vs FRDM2012 1.34 MeV  ", abs(2.0_dp * e0 - 1.34_dp), 0.05_dp)

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

end program test_zero_point
