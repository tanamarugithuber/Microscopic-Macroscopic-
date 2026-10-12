program test_shell_bcs
    !---------------------------
    ! Tests for shell_bcs_mod
    !   1. Strutinsky on a spherical oscillator spectrum (e = n + 3/2, (n+1)(n+2) states):
    !      smooth particle number, plateau in gamma_S, and the shell correction at the
    !      magic number 168 against its Thomas-Fermi value
    !   2. BCS on equidistant levels: the microscopic gap must reproduce Delta_G and the
    !      pairing correction must nearly vanish (this is how G is defined)
    !   3. 208Pb (spherical) and 240Pu (eps = 0.2): full shell + pairing correction
    !
    ! build:
    !   gfortran -O2 -fopenmp constant_mod.f90 nucleus_mod.f90 grid_mod.f90 micro_constant_mod.f90 \
    !       ho_basis_mod.f90 potential_gh_mod.f90 shell_bcs_mod.f90 test_shell_bcs.f90 -llapack -o test_shell_bcs.exe
    !---------------------------
    use iso_fortran_env, only: real64
    use nucleus_mod, only: nucleus_property
    use micro_constant_mod, only: microscopic_variables
    use shell_bcs_mod
    implicit none
    integer, parameter :: dp = real64
    logical :: ok_all

    ok_all = .true.
    call test_oscillator(ok_all)
    call test_equidistant(ok_all)
    call test_nucleus(82, 126, 0.0_dp, ok_all)
    call test_nucleus(94, 146, 0.2_dp, ok_all)

    if (ok_all) then
        print *, "ALL TESTS PASSED"
    else
        print *, "SOME TESTS FAILED"
        stop 1
    end if

contains

    subroutine check(name, err, tol, ok_all)
        character(*), intent(in) :: name
        real(dp), intent(in) :: err, tol
        logical, intent(inout) :: ok_all
        if (err <= tol) then
            print '(a, a, es10.2, a, es9.1, a)', "PASS: ", name, err, "  (tol ", tol, ")"
        else
            print '(a, a, es10.2, a, es9.1, a)', "FAIL: ", name, err, "  (tol ", tol, ")"
            ok_all = .false.
        end if
    end subroutine check

    subroutine test_oscillator(ok_all)
        !---------------------------
        ! hw = 1, shells n = 0 .. 30 with (n+1)(n+2) states (spin included), e_n = n + 3/2.
        ! Smooth part of the spectrum (Poisson summation, k = 0 term) is g(e) = e^2 - 1/4;
        ! the Euler-Maclaurin terms at the lower edge e = -1/2 (where the weight vanishes) give
        !   N(lambda) = lambda^3/3 - lambda/4
        !   E(lambda) = lambda^4/4 - lambda^2/8 + 1/64 - 1/24 + 1/120
        ! which the Strutinsky smoothing (p = 8, gamma >~ hw) must reproduce.
        !---------------------------
        logical, intent(inout) :: ok_all
        integer, parameter :: nsh = 30
        real(dp), allocatable :: e(:)
        real(dp) :: e_sum, e_smooth, lam, g, de(5), gam(5), lam_tf, e_tf, n_tf
        integer :: n, k, i, nF, n_part

        print *, "--- 1. Strutinsky, spherical oscillator ---"
        allocate(e(0))
        do n = 0, nsh
            e = [e, (n + 1.5_dp, k = 1, (n + 1) * (n + 2))]
        end do
        nF = 6
        n_part = (nF + 1) * (nF + 2) * (nF + 3) / 3      ! 168
        gam = [1.4_dp, 1.6_dp, 1.8_dp, 2.0_dp, 2.2_dp]
        do i = 1, 5
            call strutinsky_correction(e, n_part, gam(i), 8, e_sum, e_smooth, lam, g)
            de(i) = e_sum - e_smooth
            print '(a, f5.2, a, f12.6, a, f10.6, a, f10.4)', " gamma = ", gam(i), "  dE = ", de(i), &
                  "  lambda = ", lam, "  g = ", g
        end do
        call check("plateau: spread of dE over gamma ", maxval(de) - minval(de), 1.0e-3_dp, ok_all)

        n_tf = real(n_part, dp)
        lam_tf = solve_tf(n_tf)
        e_tf = lam_tf**4 / 4.0_dp - lam_tf**2 / 8.0_dp + 1.0_dp / 64.0_dp - 1.0_dp / 24.0_dp + 1.0_dp / 120.0_dp
        print '(a, f12.6, a, f12.6)', " Thomas-Fermi lambda = ", lam_tf, "   smooth energy = ", e_tf
        call strutinsky_correction(e, n_part, 1.5_dp, 8, e_sum, e_smooth, lam, g)
        call check("smooth lambda vs Thomas-Fermi    ", abs(lam - lam_tf), 1.0e-3_dp, ok_all)
        call check("smooth energy vs Thomas-Fermi    ", abs(e_smooth - e_tf), 1.0e-2_dp, ok_all)
        call check("smooth density vs lambda^2 - 1/4 ", abs(g - (lam**2 - 0.25_dp)), 1.0e-3_dp, ok_all)
    end subroutine test_oscillator

    function solve_tf(n) result(lam)
        ! lambda^3/3 - lambda/4 = n
        real(dp), intent(in) :: n
        real(dp) :: lam, f, c0
        integer :: it
        c0 = 0.0_dp
        lam = (3.0_dp * n)**(1.0_dp / 3.0_dp)
        do it = 1, 50
            f = lam**3 / 3.0_dp - lam / 4.0_dp - c0 - n
            lam = lam - f / (lam**2 - 0.25_dp)
        end do
    end function solve_tf

    subroutine test_equidistant(ok_all)
        logical, intent(inout) :: ok_all
        real(dp), parameter :: d = 0.25_dp, delta_G = 1.2_dp
        integer, parameter :: nlev = 400
        real(dp) :: e(2 * nlev), G, delta, lambda, e_pc, e_pc_ave
        integer :: k, N1, N2, n_part, iodd

        print *, "--- 2. BCS, equidistant doubly degenerate levels (spacing 0.25 MeV) ---"
        do k = 1, nlev
            e(2 * k - 1) = -60.0_dp + k * d
            e(2 * k) = e(2 * k - 1)
        end do
        do iodd = 0, 1
            n_part = 300 + iodd
            call bcs_correction(e, n_part, 1.0_dp / d, delta_G, G, delta, lambda, N1, N2, e_pc, e_pc_ave)
            print '(a, i4, a, f8.4, a, f8.4, a, 2i5, a, 2f10.4)', " N = ", n_part, "  G = ", G, &
                  "  Delta = ", delta, "  window ", N1, N2, "  E_pc, average: ", e_pc, e_pc_ave
            call check("Delta vs Delta_G (relative)      ", abs(delta / delta_G - 1.0_dp), 2.0e-2_dp, ok_all)
            call check("pairing correction |E_pair| MeV  ", abs(e_pc - e_pc_ave), 0.15_dp, ok_all)
        end do
    end subroutine test_equidistant

    subroutine test_nucleus(Z, N, eps, ok_all)
        integer, intent(in) :: Z, N
        real(dp), intent(in) :: eps
        logical, intent(inout) :: ok_all
        type(nucleus_property) :: nuc
        type(microscopic_variables) :: mv
        type(shell_bcs_result) :: res
        integer :: t0, t1, rate

        print '(a, i3, a, i3, a, f5.2, a)', " --- 3. Z = ", Z, ", N = ", N, ", eps = ", eps, " ---"
        call system_clock(t0, rate)
        nuc%Z = Z
        nuc%N = N
        nuc%A = Z + N
        nuc%epsilon2 = eps
        nuc%gamma_deg = 0.0_dp
        call shell_bcs_calculation(nuc, mv, res)
        call system_clock(t1)
        print '(a, 3f10.4)', " R_pot [fm], B_s, gamma_S [MeV]  : ", mv%R_pot, res%B_s, res%gamma_S
        call report("neutrons", res%neutron)
        call report("protons ", res%proton)
        print '(a, 3f10.4)', " E_shell, E_pair, E_micro [MeV]  : ", res%e_shell, res%e_pair, res%e_micro
        print '(a, f8.2, a)', " time: ", real(t1 - t0, dp) / rate, " s"
        if (eps == 0.0_dp) then
            ! doubly magic: large negative shell correction and no BCS solution at N = 126, Z = 82
            call check("208Pb E_shell in [-18, -10] MeV  ", &
                       max(0.0_dp, res%e_shell + 10.0_dp, -18.0_dp - res%e_shell), 0.0_dp, ok_all)
        end if
    end subroutine test_nucleus

    subroutine report(name, nc)
        character(*), intent(in) :: name
        type(nucleon_correction), intent(in) :: nc
        print '(1x, a, a, f9.3, a, f9.3, a, f8.4, a, f7.3, a, f7.3, a, f8.3, a, f8.3)', name, &
              ": E_shell ", nc%e_shell, "  lambda~ ", nc%lambda_smooth, "  G ", nc%G, &
              "  Delta_G ", nc%delta_G, "  Delta ", nc%delta, "  lambda ", nc%lambda, "  E_pair ", nc%e_pair
    end subroutine report

end program test_shell_bcs
