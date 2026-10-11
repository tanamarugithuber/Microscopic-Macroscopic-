program test_potential_gh
    !---------------------------
    ! Tests for potential_gh_mod
    !   1. sphere: Y, dY/dr, C against the analytic folded-Yukawa and Coulomb potentials
    !   2. triaxial ellipsoid: C inside against the analytic ellipsoid potential,
    !      grad Y against finite differences of Y, enclosed volume
    !   3. on_grid: octant (D2h) filling agrees with the full calculation
    !   4. 208Pb (spherical, FRDM2012 constants) in the HO basis: N=126 and Z=82 gaps
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi, e2
    use potential_gh_mod
    use ho_basis_mod
    implicit none

    integer, parameter :: dp = real64
    logical :: ok_all

    ok_all = .true.
    call test_sphere(ok_all)
    call test_ellipsoid(ok_all)
    call test_grid_symmetry(ok_all)
    call test_pb208(ok_all)

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

    subroutine test_sphere(ok_all)
        logical, intent(inout) :: ok_all
        type(ellipsoid_shape) :: sph
        type(surface_quadrature) :: q
        real(dp), parameter :: R = 6.5_dp, a = 0.8_dp
        real(dp) :: rr(14), dir(3), Y, gY(3), C, Ya, dYa, Ca, eY, edY, eC, rho, K
        integer :: i

        print *, "--- 1. sphere R = 6.5 fm, a = 0.8 fm ---"
        sph%semi = R
        call q%initialize(sph, a, 8, 16, 8)
        rr = [0.0_dp, 1.0_dp, 3.0_dp, 5.0_dp, 6.0_dp, 6.4_dp, 6.49_dp, 6.5_dp, 6.51_dp, &
              6.6_dp, 7.0_dp, 8.0_dp, 10.0_dp, 15.0_dp]
        dir = [0.3_dp, -0.5_dp, 0.81_dp]
        dir = dir / norm2(dir)
        eY = 0.0_dp
        edY = 0.0_dp
        eC = 0.0_dp
        K = (R / a) * cosh(R / a) - sinh(R / a)
        do i = 1, size(rr)
            call q%at_point(rr(i) * dir, Y, gY, C)
            rho = rr(i) / a
            if (rr(i) <= R) then
                if (rho < 1.0e-8_dp) then
                    Ya = 1.0_dp - (1.0_dp + R / a) * exp(-R / a)
                    dYa = 0.0_dp
                else
                    Ya = 1.0_dp - (1.0_dp + R / a) * exp(-R / a) * sinh(rho) / rho
                    dYa = -(1.0_dp + R / a) * exp(-R / a) * (cosh(rho) / rho - sinh(rho) / rho**2) / a
                end if
                Ca = 2.0_dp * pi * (R**2 - rr(i)**2 / 3.0_dp)
            else
                Ya = K * exp(-rho) / rho
                dYa = -K * exp(-rho) * (1.0_dp / rho + 1.0_dp / rho**2) / a
                Ca = 4.0_dp * pi * R**3 / (3.0_dp * rr(i))
            end if
            eY = max(eY, abs(Y - Ya))
            edY = max(edY, norm2(gY - dYa * dir))
            eC = max(eC, abs(C - Ca) / abs(Ca))
        end do
        call check("sphere Y, max abs error        ", eY, 1.0e-8_dp, ok_all)
        call check("sphere grad Y, max abs err 1/fm", edY, 1.0e-6_dp, ok_all)
        call check("sphere C, max rel error        ", eC, 1.0e-8_dp, ok_all)
        call check("sphere volume, rel error       ", abs(q%volume() / (4.0_dp * pi * R**3 / 3.0_dp) - 1.0_dp), &
                   1.0e-12_dp, ok_all)
    end subroutine test_sphere

    subroutine test_ellipsoid(ok_all)
        logical, intent(inout) :: ok_all
        type(ellipsoid_shape) :: ell
        type(surface_quadrature) :: q
        real(dp), parameter :: a = 0.8_dp, h = 1.0e-4_dp
        real(dp) :: pts(3, 6), Y, gY(3), C, Yp, Ym, gd(3), dum(3), Cd, eC, edY, ax2(3)
        real(dp) :: tq(200), wq(200), t, u, du, s, L
        integer :: ip, m, ic
        real(dp) :: e(3)

        print *, "--- 2. triaxial ellipsoid (eps = 0.3, gamma = 25 deg, R = 6.5 fm) ---"
        ell = ellipsoid_from_eps_gamma(6.5_dp, 0.3_dp, 25.0_dp * pi / 180.0_dp)
        print '(a, 3f9.4)', " semi-axes [fm]: ", ell%semi
        call q%initialize(ell, a, 8, 16, 8)
        call check("ellipsoid volume, rel error    ", &
                   abs(q%volume() / (4.0_dp * pi * product(ell%semi) / 3.0_dp) - 1.0_dp), 1.0e-12_dp, ok_all)

        ! points inside, near and on the surface, and outside
        pts(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        pts(:, 2) = [1.5_dp, -2.0_dp, 2.5_dp]
        pts(:, 3) = 0.98_dp * [ell%semi(1) * 0.6_dp, ell%semi(2) * 0.0_dp, ell%semi(3) * 0.8_dp]
        pts(:, 4) = 1.00_dp * [ell%semi(1) * 0.0_dp, ell%semi(2) * 0.6_dp, ell%semi(3) * 0.8_dp]
        pts(:, 5) = 1.03_dp * [ell%semi(1) * 0.48_dp, ell%semi(2) * 0.6_dp, ell%semi(3) * 0.64_dp]
        pts(:, 6) = [9.0_dp, 4.0_dp, -3.0_dp]

        ! analytic C inside a homogeneous ellipsoid
        !   C(x) = pi abc int_0^inf [1 - sum x_i^2/(a_i^2+u)] du / sqrt(prod(a_i^2+u))
        ! with u = L (1/(1-t)^2 - 1)
        call gauss_legendre(200, tq, wq)
        ax2 = ell%semi**2
        L = maxval(ax2)
        eC = 0.0_dp
        do ip = 1, 3
            Cd = 0.0_dp
            do m = 1, 200
                t = 0.5_dp * (tq(m) + 1.0_dp)
                u = L * (1.0_dp / (1.0_dp - t)**2 - 1.0_dp)
                du = 2.0_dp * L / (1.0_dp - t)**3 * 0.5_dp * wq(m)
                s = 1.0_dp - sum(pts(:, ip)**2 / (ax2 + u))
                Cd = Cd + s / sqrt(product(ax2 + u)) * du
            end do
            Cd = pi * product(ell%semi) * Cd
            call q%at_point(pts(:, ip), Y, gY, C)
            eC = max(eC, abs(C - Cd) / Cd)
        end do
        call check("ellipsoid C inside, rel error  ", eC, 1.0e-8_dp, ok_all)

        ! grad Y against central differences of Y
        edY = 0.0_dp
        do ip = 1, 6
            call q%at_point(pts(:, ip), Y, gY, C)
            do ic = 1, 3
                e = 0.0_dp
                e(ic) = h
                call q%at_point(pts(:, ip) + e, Yp, dum, C)
                call q%at_point(pts(:, ip) - e, Ym, dum, C)
                gd(ic) = (Yp - Ym) / (2.0_dp * h)
            end do
            if (ip /= 4) edY = max(edY, norm2(gY - gd))
        end do
        call check("ellipsoid grad Y vs finite diff", edY, 1.0e-6_dp, ok_all)

        ! point 4 lies on the surface, where d2Y/dn2 jumps by 1/a^2 and the central
        ! difference has an O(h) error. Check instead that grad Y is continuous there.
        edY = 0.0_dp
        call q%at_point(pts(:, 4), Y, gY, C)
        call q%at_point(pts(:, 4) * (1.0_dp + 1.0e-7_dp), Yp, gd, C)
        edY = max(edY, norm2(gY - gd))
        call q%at_point(pts(:, 4) * (1.0_dp - 1.0e-7_dp), Ym, gd, C)
        edY = max(edY, norm2(gY - gd))
        call check("grad Y continuity on surface   ", edY, 1.0e-6_dp, ok_all)
    end subroutine test_ellipsoid

    subroutine test_grid_symmetry(ok_all)
        logical, intent(inout) :: ok_all
        type(ellipsoid_shape) :: ell
        type(surface_quadrature) :: q
        real(dp), allocatable :: xi(:), W(:), xg(:), yg(:), zg(:)
        real(dp), allocatable :: Y1(:,:,:), Y2(:,:,:), G1(:,:,:,:), G2(:,:,:,:), C1(:,:,:), C2(:,:,:)
        integer, parameter :: n = 9
        real(dp) :: err

        print *, "--- 3. on_grid: D2h octant filling ---"
        ell = ellipsoid_from_eps_gamma(6.5_dp, 0.3_dp, 25.0_dp * pi / 180.0_dp)
        call q%initialize(ell, 0.8_dp, 8, 16, 8)
        allocate(xi(n), W(n))
        call gauss_hermite(n, xi, W)
        xg = 2.3_dp * xi
        yg = 2.0_dp * xi
        zg = 2.6_dp * xi
        allocate(Y1(n, n, n), Y2(n, n, n), C1(n, n, n), C2(n, n, n), G1(n, n, n, 3), G2(n, n, n, 3))
        call q%on_grid(xg, yg, zg, Y1, G1(:,:,:,1), G1(:,:,:,2), G1(:,:,:,3), C1, .true.)
        call q%on_grid(xg, yg, zg, Y2, G2(:,:,:,1), G2(:,:,:,2), G2(:,:,:,3), C2, .false.)
        err = max(maxval(abs(Y1 - Y2)), maxval(abs(G1 - G2)), maxval(abs(C1 - C2)) / maxval(abs(C2)))
        call check("octant fill vs full grid       ", err, 1.0e-10_dp, ok_all)
    end subroutine test_grid_symmetry

    subroutine test_pb208(ok_all)
        !---------------------------
        ! FRDM2012 single-particle potential for spherical 208Pb (Eqs. (81)-(92), (86), (87))
        !---------------------------
        logical, intent(inout) :: ok_all
        real(dp), parameter :: r0 = 1.16_dp, a_pot = 0.8_dp
        real(dp), parameter :: V_s = 52.5_dp, V_a = 48.7_dp, A_den = 0.82_dp, B_den = 0.56_dp
        real(dp), parameter :: C_cur = 41.0_dp, a_2 = 22.0_dp, J = 35.0_dp, L = 99.0_dp
        real(dp), parameter :: Q_stiff = 25.0_dp, K = 300.0_dp
        integer, parameter :: N_bas = 12, n_gh = 30
        real(dp) :: Z, N, A, I, c1, delta, epsb, R_den, R_pot, e_rho_c, Vp, Vn, lam_p, lam_n, gap
        type(ellipsoid_shape) :: sph
        type(surface_quadrature) :: q
        type(ho_basis_type) :: basis
        real(dp), allocatable :: xg(:), Y(:,:,:), dYx(:,:,:), dYy(:,:,:), dYz(:,:,:), C(:,:,:)
        real(dp), allocatable :: T(:,:), Vm(:,:), Axy(:,:), Ayz(:,:), Azx(:,:), e(:)
        complex(dp), allocatable :: H(:,:), U(:,:)
        integer :: ns, iq, nocc, t0, t1, rate

        print *, "--- 4. 208Pb, spherical, FRDM2012 single-particle potential ---"
        call system_clock(t0, rate)
        Z = 82.0_dp
        N = 126.0_dp
        A = Z + N
        I = (N - Z) / A
        c1 = 3.0_dp * e2 / (5.0_dp * r0)
        delta = (I + 3.0_dp / 8.0_dp * c1 / Q_stiff * Z**2 / A**(5.0_dp / 3.0_dp)) &
              / (1.0_dp + 9.0_dp / 4.0_dp * J / Q_stiff / A**(1.0_dp / 3.0_dp))
        epsb = (-2.0_dp * a_2 / A**(1.0_dp / 3.0_dp) + L * delta**2 + c1 * Z**2 / A**(4.0_dp / 3.0_dp)) / K
        R_den = r0 * A**(1.0_dp / 3.0_dp) * (1.0_dp + epsb)
        R_pot = R_den + A_den - B_den / R_den
        e_rho_c = Z * e2 / (4.0_dp * pi * R_pot**3 / 3.0_dp)
        Vp = V_s + V_a * delta
        Vn = V_s - V_a * delta
        lam_p = 0.025_dp * A + 28.0_dp
        lam_n = 0.01875_dp * A + 31.5_dp
        print '(a, 2f10.5)', " delta_bar, eps_bar : ", delta, epsb
        print '(a, 2f10.5)', " R_den, R_pot [fm]  : ", R_den, R_pot
        print '(a, 2f10.4)', " V_p, V_n [MeV]     : ", Vp, Vn

        call basis%initialize([1.0_dp, 1.0_dp, 1.0_dp] * C_cur / A**(1.0_dp / 3.0_dp), N_bas, n_gh)
        ns = basis%n_states
        sph%semi = R_pot
        call q%initialize(sph, a_pot, 8, 16, 8)

        allocate(Y(n_gh, n_gh, n_gh), dYx(n_gh, n_gh, n_gh), dYy(n_gh, n_gh, n_gh), &
                 dYz(n_gh, n_gh, n_gh), C(n_gh, n_gh, n_gh))
        ! spherical basis: b_x = b_y = b_z
        xg = basis%b(1) * basis%xi
        call q%on_grid(xg, xg, xg, Y, dYx, dYy, dYz, C, .true.)
        call system_clock(t1)
        print '(a, f8.2, a)', " potentials on the GH grid: ", real(t1 - t0, dp) / rate, " s"

        allocate(T(ns, ns), Vm(ns, ns), Axy(ns, ns), Ayz(ns, ns), Azx(ns, ns))
        allocate(H(2 * ns, 2 * ns), U(2 * ns, 2 * ns), e(2 * ns))
        call basis%kinetic_matrix(T)

        do iq = 1, 2
            if (iq == 1) then
                ! neutrons
                call basis%potential_matrix(-Vn * Y, Vm)
                call basis%spin_orbit_matrices(-Vn * dYx, -Vn * dYy, -Vn * dYz, Axy, Ayz, Azx)
                call basis%build_hamiltonian(T, Vm, Axy, Ayz, Azx, so_strength(lam_n), H)
                nocc = nint(N)
            else
                ! protons
                call basis%potential_matrix(-Vp * Y + e_rho_c * C, Vm)
                call basis%spin_orbit_matrices(-Vp * dYx, -Vp * dYy, -Vp * dYz, Axy, Ayz, Azx)
                call basis%build_hamiltonian(T, Vm, Axy, Ayz, Azx, so_strength(lam_p), H)
                nocc = nint(Z)
            end if
            call diagonalize_hamiltonian(H, e, U)
            gap = e(nocc + 1) - e(nocc)
            if (iq == 1) then
                print '(a, 2f10.4)', " neutrons: e(126), e(127) [MeV] = ", e(nocc), e(nocc + 1)
                call print_levels(e, nocc)
                call check("N = 126 gap shortfall from 2.5 ", max(0.0_dp, 2.5_dp - gap), 0.0_dp, ok_all)
            else
                print '(a, 2f10.4)', " protons : e(82),  e(83)  [MeV] = ", e(nocc), e(nocc + 1)
                call print_levels(e, nocc)
                call check("Z = 82 gap shortfall from 2.5  ", max(0.0_dp, 2.5_dp - gap), 0.0_dp, ok_all)
            end if
            print '(a, f8.4, a)', " gap = ", gap, " MeV"
        end do
        call system_clock(t1)
        print '(a, f8.2, a)', " total time: ", real(t1 - t0, dp) / rate, " s"
    end subroutine test_pb208

    subroutine print_levels(e, nocc)
        ! levels around the Fermi surface, grouped into (nearly) degenerate multiplets.
        ! The spread within a multiplet comes from the Gauss-Hermite quadrature of
        ! grad V1, which has a kink at the surface (see the README).
        real(dp), intent(in) :: e(:)
        integer, intent(in) :: nocc
        integer :: i, i0, nfill
        i = 1
        nfill = 0
        do while (i <= size(e))
            i0 = i
            do while (i < size(e))
                if (e(i + 1) - e(i0) > 0.1_dp) exit
                i = i + 1
            end do
            nfill = nfill + (i - i0 + 1)
            if (nfill > nocc - 30 .and. nfill <= nocc + 20) &
                print '(a, i3, a, f10.4, a, f7.4, a, i4)', "   degeneracy ", i - i0 + 1, "  e = ", &
                    sum(e(i0:i)) / (i - i0 + 1), "  spread = ", e(i) - e(i0), "   filled up to ", nfill
            i = i + 1
        end do
    end subroutine print_levels

end program test_potential_gh
