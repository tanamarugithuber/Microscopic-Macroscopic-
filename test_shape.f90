program test_shape
    !---------------------------
    ! Tests for shape_mod
    !   1. epsilon shape with eps4 = eps6 = 0 is the ellipsoid of ellipsoid_from_eps_gamma
    !   2. volume conservation, inside() consistent with eval(), N = dr/du x dr/dv
    !   3. beta_l of the FRDM2012 ground state of 240Pu (eps2 = 0.22, eps4 = -0.05, eps6 = 0.04)
    !      against the tabulated beta2 = 0.237, beta4 = 0.086, beta6 = -0.024
    !
    ! build:
    !   gfortran -O2 constant_mod.f90 shape_mod.f90 test_shape.f90 -o test_shape.exe
    !---------------------------
    ! 【日本語】shape_mod のテスト: epsilon 形状が楕円体と一致すること, 体積保存,
    ! inside と eval の整合, N = dr/du x dr/dv, 240Pu の beta_l が FRDM2012 の表と合うこと.
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use shape_mod
    implicit none
    integer, parameter :: dp = real64
    logical :: ok_all
    type(epsilon_shape) :: es
    type(ellipsoid_shape) :: el
    real(dp) :: R, err, r1(3), N1(3), rp(3), rm(3), dru(3), drv(3), Nfd(3), beta(3)
    real(dp), parameter :: hd = 1.0e-6_dp
    integer :: i, j

    ok_all = .true.
    R = 7.0_dp

    print *, "--- 1. eps4 = eps6 = 0 is the ellipsoid ---"
    es = epsilon_shape_from(R, 0.3_dp, 20.0_dp * pi / 180.0_dp, 0.0_dp, 0.0_dp)
    el = ellipsoid_from_eps_gamma(R, 0.3_dp, 20.0_dp * pi / 180.0_dp)
    call check("extent vs semi-axes (fm)         ", maxval(abs(es%extent() - el%semi)), 1.0e-3_dp)
    call check("relative area vs ellipsoid       ", abs(es%relative_area() - el%relative_area()), 1.0e-10_dp)
    err = 0.0_dp
    do j = 0, 12
        do i = 1, 11
            call es%eval(pi * i / 12.0_dp, 2.0_dp * pi * j / 12.0_dp, r1, N1)
            err = max(err, abs(sum((r1 / el%semi)**2) - 1.0_dp))
        end do
    end do
    call check("epsilon surface on the ellipsoid ", err, 1.0e-12_dp)

    print *, "--- 2. eps2 = 0.25, gamma = 15 deg, eps4 = -0.06, eps6 = 0.03 ---"
    es = epsilon_shape_from(R, 0.25_dp, 15.0_dp * pi / 180.0_dp, -0.06_dp, 0.03_dp)
    call check("volume / (4/3 pi R^3) - 1        ", abs(es%volume() / (4.0_dp * pi * R**3 / 3.0_dp) - 1.0_dp), &
               1.0e-10_dp)
    err = 0.0_dp
    do j = 0, 8
        do i = 1, 7
            call es%eval(pi * i / 8.0_dp, 2.0_dp * pi * j / 8.0_dp + 0.1_dp, r1, N1)
            if (.not. es%inside(r1 * (1.0_dp - 1.0e-7_dp)) .or. es%inside(r1 * (1.0_dp + 1.0e-7_dp))) err = 1.0_dp
            call es%eval(pi * i / 8.0_dp + hd, 2.0_dp * pi * j / 8.0_dp + 0.1_dp, rp, N1)
            call es%eval(pi * i / 8.0_dp - hd, 2.0_dp * pi * j / 8.0_dp + 0.1_dp, rm, N1)
            dru = (rp - rm) / (2.0_dp * hd)
            call es%eval(pi * i / 8.0_dp, 2.0_dp * pi * j / 8.0_dp + 0.1_dp + hd, rp, N1)
            call es%eval(pi * i / 8.0_dp, 2.0_dp * pi * j / 8.0_dp + 0.1_dp - hd, rm, N1)
            drv = (rp - rm) / (2.0_dp * hd)
            Nfd = [dru(2) * drv(3) - dru(3) * drv(2), dru(3) * drv(1) - dru(1) * drv(3), dru(1) * drv(2) - dru(2) * drv(1)]
            call es%eval(pi * i / 8.0_dp, 2.0_dp * pi * j / 8.0_dp + 0.1_dp, r1, N1)
            err = max(err, norm2(N1 - Nfd) / norm2(Nfd))
        end do
    end do
    call check("inside() and N vs finite diff    ", err, 1.0e-6_dp)

    print *, "--- 3. 240Pu ground state, eps2 = 0.22, eps4 = -0.05, eps6 = 0.04 ---"
    es = epsilon_shape_from(R, 0.22_dp, 0.0_dp, -0.05_dp, 0.04_dp)
    beta = betas(es)
    print '(a, 3f9.4, a)', " beta2, beta4, beta6 = ", beta, "   (FRDM2012: 0.237 0.086 -0.024)"
    call check("|beta - FRDM2012| max            ", maxval(abs(beta - [0.237_dp, 0.086_dp, -0.024_dp])), 1.0e-2_dp)

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

    function betas(s) result(b)
        ! r(theta) = R0 [1 + sum beta_l Y_l0],  beta_l = 4 pi int r Y_l0 dOmega / int r dOmega
        ! (axial shape; radius along each lab direction by bisection on inside())
        type(epsilon_shape), intent(in) :: s
        real(dp) :: b(3), tq(64), wq(64), th, rr, lo, hi, r0s, Y(3), x(3)
        integer :: i, it, l
        call gauss_legendre(64, tq, wq)
        b = 0.0_dp
        r0s = 0.0_dp
        do i = 1, 64
            th = 0.5_dp * pi * (tq(i) + 1.0_dp)
            x = [sin(th), 0.0_dp, cos(th)]
            lo = 0.0_dp
            hi = 20.0_dp
            do it = 1, 60
                rr = 0.5_dp * (lo + hi)
                if (s%inside(rr * x)) then
                    lo = rr
                else
                    hi = rr
                end if
            end do
            do l = 1, 3
                Y(l) = legendre(2 * l, cos(th)) * sqrt((4.0_dp * l + 1.0_dp) / (4.0_dp * pi))
            end do
            b = b + wq(i) * sin(th) * rr * Y
            r0s = r0s + wq(i) * sin(th) * rr
        end do
        b = 4.0_dp * pi * b / r0s
    end function betas

    pure function legendre(l, x) result(p)
        integer, intent(in) :: l
        real(dp), intent(in) :: x
        real(dp) :: p, p0, p1
        integer :: k
        p0 = 1.0_dp
        p1 = x
        if (l == 0) then
            p = 1.0_dp
            return
        end if
        do k = 2, l
            p = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            p0 = p1
            p1 = p
        end do
        p = p1
    end function legendre

end program test_shape
