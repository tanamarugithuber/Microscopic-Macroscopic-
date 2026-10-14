program test_hockney
    !---------------------------
    ! Tests for hockney_mod against analytic potentials of Gaussian charges
    !   Coulomb : Q erf(d / (sqrt(2) s)) / d
    !   Yukawa  : Q/(2d) exp(k^2 s^2/2) [exp(-k d) erfc((k s^2 - d)/(sqrt(2) s))
    !                                   - exp(k d) erfc((k s^2 + d)/(sqrt(2) s))]
    ! three off-centre Gaussians (no symmetry), cubic and non-cubic boxes, smooth = .true.
    ! and a uniform cube made of whole cells (smooth = .false.), whose Coulomb potential
    ! at the centre is 2.3800773... L^2. The exp kernel is checked against -d/dkappa of the Yukawa result.
    ! The box is cut tight (charge reaches the edge at ~1e-4) to show that the
    ! boundary does not matter.
    !---------------------------
    ! 【日本語】hockney_mod のテスト: ガウス型電荷のクーロン・湯川ポテンシャルの解析解,
    ! 一様な立方体の中心でのクーロンポテンシャル (2.3800773 L^2) と比べる.
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use hockney_mod
    implicit none

    integer, parameter :: dp = real64
    logical :: ok_all, exp_kernel = .false.
    real(dp) :: err_max, err_w

    ok_all = .true.

    print *, "--- 1. Coulomb, 51^3, h = 0.3 fm ---"
    call run(51, 51, 51, 0.3_dp, 0.0_dp, .true., err_max, err_w)
    call check("max |V - V_exact| / max V          ", err_max, 5.0e-5_dp)
    call check("int rho (V - V_exact) / int rho V   ", err_w, 2.0e-5_dp)

    print *, "--- 2. Coulomb, non-cubic 41 x 51 x 61, h = 0.3 fm ---"
    call run(41, 51, 61, 0.3_dp, 0.0_dp, .true., err_max, err_w)
    call check("max |V - V_exact| / max V          ", err_max, 5.0e-5_dp)
    call check("int rho (V - V_exact) / int rho V   ", err_w, 2.0e-5_dp)

    print *, "--- 3. Yukawa a = 0.68 fm, 51^3, h = 0.3 fm ---"
    call run(51, 51, 51, 0.3_dp, 1.0_dp / 0.68_dp, .true., err_max, err_w)
    call check("max |V - V_exact| / max V          ", err_max, 5.0e-4_dp)
    call check("int rho (V - V_exact) / int rho V   ", err_w, 2.0e-4_dp)

    print *, "--- 4. Yukawa a = 0.70 fm, non-cubic 81 x 101 x 121, h = 0.15 fm ---"
    call run(81, 101, 121, 0.15_dp, 1.0_dp / 0.70_dp, .true., err_max, err_w)
    call check("max |V - V_exact| / max V          ", err_max, 1.0e-4_dp)
    call check("int rho (V - V_exact) / int rho V   ", err_w, 2.0e-5_dp)

    print *, "--- 5. exp(-r/a) kernel, a = 0.70 fm, 51^3, h = 0.3 fm ---"
    exp_kernel = .true.
    call run(51, 51, 51, 0.3_dp, 1.0_dp / 0.70_dp, .true., err_max, err_w)
    exp_kernel = .false.
    call check("max |W - W_exact| / max W          ", err_max, 1.0e-4_dp)
    call check("int rho (W - W_exact) / int rho W   ", err_w, 2.0e-5_dp)

    print *, "--- 6. uniform cube of 21^3 cells (cell-constant density), Coulomb ---"
    call cube_test(err_max)
    call check("|V(centre) - 2.3800773 L^2| / V     ", err_max, 1.0e-6_dp)

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
            print '(a, a, es10.2, a, es9.1, a)', "PASS: ", name, err, "  (tol ", tol, ")"
        else
            print '(a, a, es10.2, a, es9.1, a)', "FAIL: ", name, err, "  (tol ", tol, ")"
            ok_all = .false.
        end if
    end subroutine check

    subroutine run(nx, ny, nz, h, kappa, smooth, err_max, err_w)
        integer, intent(in) :: nx, ny, nz
        real(dp), intent(in) :: h, kappa
        logical, intent(in) :: smooth
        real(dp), intent(out) :: err_max, err_w
        real(dp), parameter :: c(3, 3) = reshape([1.0_dp, 0.5_dp, -0.8_dp, &
                                                 -1.5_dp, 0.3_dp, 1.0_dp, &
                                                  0.2_dp, -1.2_dp, 0.4_dp], [3, 3])
        real(dp), parameter :: sig(3) = [1.2_dp, 1.0_dp, 1.5_dp]
        real(dp), parameter :: Q(3) = [1.0_dp, 0.7_dp, 0.5_dp]
        real(dp), allocatable :: rho(:,:,:), V(:,:,:), Vex(:,:,:)
        real(dp) :: p(3), d
        integer :: i, j, k, g

        allocate(rho(nx, ny, nz), V(nx, ny, nz), Vex(nx, ny, nz))
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    p = [(i - 0.5_dp * (nx + 1)) * h, (j - 0.5_dp * (ny + 1)) * h, (k - 0.5_dp * (nz + 1)) * h]
                    rho(i, j, k) = 0.0_dp
                    Vex(i, j, k) = 0.0_dp
                    do g = 1, 3
                        d = norm2(p - c(:, g))
                        rho(i, j, k) = rho(i, j, k) + Q(g) / ((2.0_dp * pi)**1.5_dp * sig(g)**3) &
                                     * exp(-0.5_dp * d * d / sig(g)**2)
                        if (exp_kernel) then
                            ! int rho exp(-kappa r) = - d/dkappa int rho exp(-kappa r)/r
                            Vex(i, j, k) = Vex(i, j, k) - Q(g) * (gauss_potential(d, sig(g), kappa + 1.0e-4_dp) &
                                         - gauss_potential(d, sig(g), kappa - 1.0e-4_dp)) / 2.0e-4_dp
                        else
                            Vex(i, j, k) = Vex(i, j, k) + Q(g) * gauss_potential(d, sig(g), kappa)
                        end if
                    end do
                end do
            end do
        end do
        if (exp_kernel) then
            call hockney_exp(rho, h, 1.0_dp / kappa, V, smooth)
        else if (kappa == 0.0_dp) then
            call hockney_coulomb(rho, h, V, smooth)
        else
            call hockney_yukawa(rho, h, 1.0_dp / kappa, V, smooth)
        end if
        err_max = maxval(abs(V - Vex)) / maxval(Vex)
        err_w = abs(sum(rho * (V - Vex)) / sum(rho * Vex))
        print '(a, f10.6, a, f10.6)', " V at the box centre: ", V((nx+1)/2, (ny+1)/2, (nz+1)/2), &
              "   exact: ", Vex((nx+1)/2, (ny+1)/2, (nz+1)/2)
        print '(a, f10.6, a, f10.6)', " V at a box corner  : ", V(1, 1, 1), "   exact: ", Vex(1, 1, 1)
    end subroutine run

    subroutine cube_test(err)
        real(dp), intent(out) :: err
        integer, parameter :: n = 41, m = 10
        real(dp), parameter :: h = 0.3_dp
        real(dp) :: rho(n, n, n), V(n, n, n), Vex
        integer :: c
        c = (n + 1) / 2
        rho = 0.0_dp
        rho(c - m:c + m, c - m:c + m, c - m:c + m) = 1.0_dp
        call hockney_coulomb(rho, h, V)
        Vex = 2.3800773639795_dp * ((2 * m + 1) * h)**2
        print '(a, f12.8, a, f12.8)', " V at the centre: ", V(c, c, c), "   exact: ", Vex
        err = abs(V(c, c, c) - Vex) / Vex
    end subroutine cube_test

    function gauss_potential(d, s, kappa) result(f)
        ! potential of a unit Gaussian charge of width s at distance d
        real(dp), intent(in) :: d, s, kappa
        real(dp) :: f, t
        t = sqrt(2.0_dp) * s
        if (kappa == 0.0_dp) then
            if (d < 1.0e-8_dp) then
                f = sqrt(2.0_dp / pi) / s
            else
                f = erf(d / t) / d
            end if
        else
            if (d < 1.0e-8_dp) then
                ! limit d -> 0
                f = sqrt(2.0_dp / pi) / s - kappa * exp(0.5_dp * (kappa * s)**2) * erfc(kappa * s / sqrt(2.0_dp))
            else
                f = 0.5_dp / d * (exp(0.5_dp * (kappa * s)**2 - kappa * d) * erfc((kappa * s * s - d) / t) &
                                - exp(0.5_dp * (kappa * s)**2 + kappa * d) * erfc((kappa * s * s + d) / t))
            end if
        end if
    end function gauss_potential

end program test_hockney
