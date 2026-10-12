program test_coulomb_fft
    !---------------------------
    ! Tests for calculate_coulomb_FFT_mod
    !   1. three off-centre Gaussian charges (non-zero dipole and quadrupole moments)
    !      against the analytic potential  Q erf(d / (sqrt(2) sigma)) / d
    !   2. independence of the reference density (R_pot, a_ref)
    !   3. non-cubic box
    !   4. D2h-symmetric charge (no odd multipoles), the case met for nuclear shapes
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use calculate_coulomb_FFT_mod
    implicit none

    integer, parameter :: dp = real64
    logical :: ok_all
    real(dp) :: err_max, err_w, err2

    ok_all = .true.

    print *, "--- 1. Gaussian charges, cubic box 71^3, h = 0.3 fm ---"
    call run(71, 71, 71, 0.3_dp, 4.0_dp, 0.5_dp, .false., err_max, err_w)
    call check("max |V - V_exact| / max V          ", err_max, 1.0e-4_dp)
    call check("int rho (V - V_exact) / int rho V   ", err_w, 3.0e-6_dp)

    print *, "--- 2. other reference density (R = 2.5 fm, a = 0.8 fm) ---"
    call run(71, 71, 71, 0.3_dp, 2.5_dp, 0.8_dp, .false., err2, err_w)
    call check("max |V - V_exact| / max V          ", err2, 1.0e-4_dp)

    print *, "--- 3. non-cubic box 61 x 71 x 81, h = 0.3 fm ---"
    call run(61, 71, 81, 0.3_dp, 4.0_dp, 0.5_dp, .false., err_max, err_w)
    call check("max |V - V_exact| / max V          ", err_max, 1.0e-4_dp)
    call check("int rho (V - V_exact) / int rho V   ", err_w, 3.0e-6_dp)

    print *, "--- 4. D2h-symmetric charge, 61 x 71 x 81, h = 0.3 fm ---"
    call run(61, 71, 81, 0.3_dp, 4.0_dp, 0.5_dp, .true., err_max, err_w)
    call check("max |V - V_exact| / max V          ", err_max, 1.0e-5_dp)
    call check("int rho (V - V_exact) / int rho V   ", err_w, 3.0e-6_dp)

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

    subroutine run(nx, ny, nz, h, R_ref, a_ref, symmetric, err_max, err_w)
        integer, intent(in) :: nx, ny, nz
        real(dp), intent(in) :: h, R_ref, a_ref
        logical, intent(in) :: symmetric
        real(dp), intent(out) :: err_max, err_w
        real(dp), parameter :: c_asym(3, 3) = reshape([1.0_dp, 0.5_dp, -0.8_dp, &
                                                 -1.5_dp, 0.3_dp, 1.0_dp, &
                                                  0.2_dp, -1.2_dp, 0.4_dp], [3, 3])
        real(dp), parameter :: sig_asym(3) = [1.2_dp, 1.0_dp, 1.5_dp]
        real(dp), parameter :: Q_asym(3) = [1.0_dp, 0.7_dp, 0.5_dp]
        ! D2h: pairs at (+-1.5, 0, 0), (0, +-1, 0), (0, 0, +-2) and one at the centre
        real(dp), parameter :: c_sym(3, 7) = reshape([1.5_dp, 0.0_dp, 0.0_dp, -1.5_dp, 0.0_dp, 0.0_dp, &
                                                      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, &
                                                      0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, -2.0_dp, &
                                                      0.0_dp, 0.0_dp, 0.0_dp], [3, 7])
        real(dp), parameter :: sig_sym(7) = [1.2_dp, 1.2_dp, 1.0_dp, 1.0_dp, 1.3_dp, 1.3_dp, 1.5_dp]
        real(dp), parameter :: Q_sym(7) = [0.5_dp, 0.5_dp, 0.3_dp, 0.3_dp, 0.4_dp, 0.4_dp, 1.0_dp]
        real(dp), allocatable :: c(:,:), sig(:), Q(:)
        real(dp), allocatable :: rho(:,:,:), V(:,:,:), Vex(:,:,:)
        real(dp) :: p(3), d, num, den
        integer :: i, j, k, g

        if (symmetric) then
            c = c_sym; sig = sig_sym; Q = Q_sym
        else
            c = c_asym; sig = sig_asym; Q = Q_asym
        end if
        allocate(rho(nx, ny, nz), V(nx, ny, nz), Vex(nx, ny, nz))
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    p = [(i - 0.5_dp * (nx + 1)) * h, (j - 0.5_dp * (ny + 1)) * h, (k - 0.5_dp * (nz + 1)) * h]
                    rho(i, j, k) = 0.0_dp
                    Vex(i, j, k) = 0.0_dp
                    do g = 1, size(Q)
                        d = norm2(p - c(:, g))
                        rho(i, j, k) = rho(i, j, k) + Q(g) / ((2.0_dp * pi)**1.5_dp * sig(g)**3) &
                                     * exp(-0.5_dp * d * d / sig(g)**2)
                        if (d < 1.0e-12_dp) then
                            Vex(i, j, k) = Vex(i, j, k) + Q(g) * sqrt(2.0_dp / pi) / sig(g)
                        else
                            Vex(i, j, k) = Vex(i, j, k) + Q(g) * erf(d / (sqrt(2.0_dp) * sig(g))) / d
                        end if
                    end do
                end do
            end do
        end do
        call calculate_coulomb(rho, R_ref, h, V, a_ref)
        err_max = maxval(abs(V - Vex)) / maxval(Vex)
        print '(a, es10.2)', " max rel. error inside |r| < 8 fm: ", inner_error(V, Vex, nx, ny, nz, h, 8.0_dp)
        num = sum(rho * (V - Vex))
        den = sum(rho * Vex)
        err_w = abs(num / den)
        print '(a, f10.6, a, f10.6)', " V at the box centre: ", V((nx+1)/2, (ny+1)/2, (nz+1)/2), &
              "   exact: ", Vex((nx+1)/2, (ny+1)/2, (nz+1)/2)
        print '(a, f10.6, a, f10.6)', " V at a box corner  : ", V(1, 1, 1), "   exact: ", Vex(1, 1, 1)
    end subroutine run

    function inner_error(V, Vex, nx, ny, nz, h, rcut) result(e)
        real(dp), intent(in) :: V(:,:,:), Vex(:,:,:), h, rcut
        integer, intent(in) :: nx, ny, nz
        real(dp) :: e, r
        integer :: i, j, k
        e = 0.0_dp
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    r = h * sqrt((i - 0.5_dp * (nx + 1))**2 + (j - 0.5_dp * (ny + 1))**2 + (k - 0.5_dp * (nz + 1))**2)
                    if (r < rcut) e = max(e, abs(V(i, j, k) - Vex(i, j, k)))
                end do
            end do
        end do
        e = e / maxval(Vex)
    end function inner_error

end program test_coulomb_fft
