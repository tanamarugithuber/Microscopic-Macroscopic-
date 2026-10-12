module shell_bcs_mod
    !---------------------------
    ! Microscopic corrections: Strutinsky shell correction + BCS pairing correction
    !   P. Moller et al., ADNDT 109-110 (2016) 1  (FRDM2012), Secs. 2.7-2.10
    !   P. Moller and J. R. Nix, Nucl. Phys. A 536 (1992) 20  (thesis/Nuclear_Pairing_Model.pdf)
    !
    ! 1. single-particle levels
    !    V = V1 + V_so + V_C (Eqs. 80-92) for the (eps, gamma) ellipsoid of radius R_pot,
    !    V1, grad V1, V_C on the Gauss-Hermite points by surface integrals (potential_gh_mod),
    !    diagonalised in the deformed HO basis with hw0 = C_cur / A^(1/3) (ho_basis_mod).
    ! 2. Strutinsky shell correction (Eq. 78)
    !    E_shell = sum_{i <= N} e_i - E_smooth,
    !    smoothing function: Gaussian times the correction polynomial of order p (Hermite form)
    !      f(x) = exp(-x^2)/sqrt(pi) * sum_{m=0}^{p/2} c_m H_2m(x),  c_m = (-1)^m / (4^m m!)
    !    range gamma_S = C_S hw0 B_s (Eq. 111), smooth Fermi level from N_smooth(lambda) = N.
    ! 3. BCS pairing correction (Eq. 79)
    !    E_pair = E_pc - E_pc_average,
    !    E_pc (Moller-Nix Eq. 5) from the BCS equations between levels N1 and N2,
    !    G from the effective-interaction gap Delta_G = r_mic B_s / N^(1/3) (Eqs. 102-107),
    !    E_pc_average (Moller-Nix Eq. 27) for the same G, Delta_G and level density.
    !    Pairing window (Moller-Nix Sec. 2.5): N2 = last bound level or last level within
    !    5 MeV above the Fermi level, whichever is higher; the same number of levels below.
    !    Odd particle numbers: no blocking, quasi-particle energy of the odd level added.
    !
    ! Each eigenvalue of the 2 n_states x 2 n_states Hamiltonian holds one nucleon;
    ! Kramers pairs (time-reversal partners) form the doubly degenerate levels of BCS.
    ! Constants are taken from micro_constant_mod.
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi, e2
    use nucleus_mod, only: nucleus_property
    use ho_basis_mod, only: ho_basis_type, frequencies_from_eps_gamma, so_strength, diagonalize_hamiltonian
    use potential_gh_mod, only: ellipsoid_shape, ellipsoid_from_eps_gamma, surface_quadrature, gauss_legendre
    use micro_constant_mod, only: microscopic_variables, C_cur, V_s, V_a, mi_A_den, mi_B_den, a_pot, &
                                  k_p, l_p, k_n, l_n, N_bas, mi_p, C_s, r_mic, a_2, mi_J, mi_L, mi_Q, mi_K
    implicit none
    private

    integer, parameter :: dp = real64

    !---------------------------
    ! results for one kind of nucleon
    !---------------------------
    type, public :: nucleon_correction
        real(dp) :: e_sum = 0.0_dp          ! MeV, sum of the occupied single-particle energies
        real(dp) :: e_smooth = 0.0_dp       ! MeV, Strutinsky smooth energy
        real(dp) :: lambda_smooth = 0.0_dp  ! MeV, Fermi level of the smooth distribution
        real(dp) :: g_smooth = 0.0_dp       ! 1/MeV, smooth level density at lambda_smooth (both spins)
        real(dp) :: e_shell = 0.0_dp        ! MeV, shell correction
        real(dp) :: delta_G = 0.0_dp        ! MeV, effective-interaction pairing gap
        real(dp) :: G = 0.0_dp              ! MeV, pairing strength
        real(dp) :: delta = 0.0_dp          ! MeV, BCS pairing gap (0: no non-trivial solution)
        real(dp) :: lambda = 0.0_dp         ! MeV, BCS Fermi level
        integer :: N1 = 0, N2 = 0           ! pairing window (doubly degenerate levels)
        real(dp) :: e_pc = 0.0_dp           ! MeV, pairing correlation + quasi-particle energy
        real(dp) :: e_pc_ave = 0.0_dp       ! MeV, its average
        real(dp) :: e_pair = 0.0_dp         ! MeV, pairing correction
        real(dp), allocatable :: e(:)       ! MeV, single-particle energies (ascending, one nucleon each)
    end type nucleon_correction

    type, public :: shell_bcs_result
        type(nucleon_correction) :: neutron, proton
        real(dp) :: B_s = 1.0_dp            ! relative surface area of the (eps, gamma) shape
        real(dp) :: gamma_S = 0.0_dp        ! MeV, Strutinsky range
        real(dp) :: e_shell = 0.0_dp        ! MeV, neutron + proton shell correction
        real(dp) :: e_pair = 0.0_dp         ! MeV, neutron + proton pairing correction
        real(dp) :: e_micro = 0.0_dp        ! MeV, e_shell + e_pair
    end type shell_bcs_result

    public :: shell_bcs_calculation
    public :: strutinsky_correction
    public :: bcs_correction
    public :: ellipsoid_relative_area

    contains

        subroutine shell_bcs_calculation(nucleus, mv, res, n_gh)
            !---------------------------
            ! shell + pairing correction of the nucleus at (nucleus%epsilon2, nucleus%gamma_deg).
            ! mv receives the potential quantities (delta_ave, epsilon_ave, R_den, R_pot,
            ! e_rho_c, hbar_omega0, lambda_p, lambda_n).
            ! n_gh: Gauss-Hermite points per axis (default 30)
            !---------------------------
            implicit none
            type(nucleus_property), intent(in) :: nucleus
            type(microscopic_variables), intent(inout) :: mv
            type(shell_bcs_result), intent(out) :: res
            integer, intent(in), optional :: n_gh

            real(dp) :: Z, N, A, I, c1, eps, gam, Vp, Vn
            integer :: ngh, ns
            type(ho_basis_type) :: basis
            type(ellipsoid_shape) :: shape
            type(surface_quadrature) :: q
            real(dp), allocatable :: xg(:), yg(:), zg(:)
            real(dp), allocatable :: Y(:,:,:), dYx(:,:,:), dYy(:,:,:), dYz(:,:,:), C(:,:,:)
            real(dp), allocatable :: T(:,:), Vm(:,:), Axy(:,:), Ayz(:,:), Azx(:,:)
            complex(dp), allocatable :: H(:,:), U(:,:)

            ngh = 30
            if (present(n_gh)) ngh = n_gh

            Z = real(nucleus%Z, dp)
            N = real(nucleus%N, dp)
            A = real(nucleus%A, dp)
            eps = nucleus%epsilon2
            gam = nucleus%gamma_deg * pi / 180.0_dp

            !---------------------------
            ! potential quantities, Eqs. (82)-(87), (89), (90), (92), (93)
            !---------------------------
            I = (N - Z) / A
            c1 = 3.0_dp * e2 / (5.0_dp * nucleus%R0)
            mv%delta_ave = (I + 3.0_dp / 8.0_dp * c1 / mi_Q * Z**2 / A**(5.0_dp / 3.0_dp)) &
                         / (1.0_dp + 9.0_dp / 4.0_dp * mi_J / mi_Q / A**(1.0_dp / 3.0_dp))
            mv%epsilon_ave = (-2.0_dp * a_2 / A**(1.0_dp / 3.0_dp) + mi_L * mv%delta_ave**2 &
                            + c1 * Z**2 / A**(4.0_dp / 3.0_dp)) / mi_K
            mv%R_den = nucleus%R0 * A**(1.0_dp / 3.0_dp) * (1.0_dp + mv%epsilon_ave)
            mv%R_pot = mv%R_den + mi_A_den - mi_B_den / mv%R_den
            mv%e_rho_c = Z * e2 / (4.0_dp * pi * mv%R_pot**3 / 3.0_dp)
            mv%hbar_omega0 = C_cur / A**(1.0_dp / 3.0_dp)
            mv%lambda_p = k_p * A + l_p
            mv%lambda_n = k_n * A + l_n
            Vp = V_s + V_a * mv%delta_ave
            Vn = V_s - V_a * mv%delta_ave

            !---------------------------
            ! basis, shape and potentials on the Gauss-Hermite points
            !---------------------------
            call basis%initialize(frequencies_from_eps_gamma(mv%hbar_omega0, eps, gam), nint(N_bas), ngh)
            ns = basis%n_states
            shape = ellipsoid_from_eps_gamma(mv%R_pot, eps, gam)
            call q%initialize(shape, a_pot, 8, 16, 8)
            xg = basis%b(1) * basis%xi
            yg = basis%b(2) * basis%xi
            zg = basis%b(3) * basis%xi
            allocate(Y(ngh, ngh, ngh), dYx(ngh, ngh, ngh), dYy(ngh, ngh, ngh), dYz(ngh, ngh, ngh), C(ngh, ngh, ngh))
            call q%on_grid(xg, yg, zg, Y, dYx, dYy, dYz, C, .true.)

            !---------------------------
            ! single-particle spectra
            !---------------------------
            allocate(T(ns, ns), Vm(ns, ns), Axy(ns, ns), Ayz(ns, ns), Azx(ns, ns))
            allocate(H(2 * ns, 2 * ns), U(2 * ns, 2 * ns))
            allocate(res%neutron%e(2 * ns), res%proton%e(2 * ns))
            call basis%kinetic_matrix(T)

            call basis%potential_matrix(-Vn * Y, Vm)
            call basis%spin_orbit_matrices(-Vn * dYx, -Vn * dYy, -Vn * dYz, Axy, Ayz, Azx)
            call basis%build_hamiltonian(T, Vm, Axy, Ayz, Azx, so_strength(mv%lambda_n), H)
            call diagonalize_hamiltonian(H, res%neutron%e, U)

            call basis%potential_matrix(-Vp * Y + mv%e_rho_c * C, Vm)
            call basis%spin_orbit_matrices(-Vp * dYx, -Vp * dYy, -Vp * dYz, Axy, Ayz, Azx)
            call basis%build_hamiltonian(T, Vm, Axy, Ayz, Azx, so_strength(mv%lambda_p), H)
            call diagonalize_hamiltonian(H, res%proton%e, U)

            !---------------------------
            ! Strutinsky and BCS
            !---------------------------
            res%B_s = ellipsoid_relative_area(ellipsoid_from_eps_gamma(1.0_dp, eps, gam))
            res%gamma_S = C_s * mv%hbar_omega0 * res%B_s
            call correct_one(res%neutron, nucleus%N)
            call correct_one(res%proton, nucleus%Z)
            res%e_shell = res%neutron%e_shell + res%proton%e_shell
            res%e_pair = res%neutron%e_pair + res%proton%e_pair
            res%e_micro = res%e_shell + res%e_pair

        contains

            subroutine correct_one(nc, n_part)
                type(nucleon_correction), intent(inout) :: nc
                integer, intent(in) :: n_part
                call strutinsky_correction(nc%e, n_part, res%gamma_S, nint(mi_p), &
                                           nc%e_sum, nc%e_smooth, nc%lambda_smooth, nc%g_smooth)
                nc%e_shell = nc%e_sum - nc%e_smooth
                nc%delta_G = r_mic * res%B_s / real(n_part, dp)**(1.0_dp / 3.0_dp)
                call bcs_correction(nc%e, n_part, 0.5_dp * nc%g_smooth, nc%delta_G, &
                                    nc%G, nc%delta, nc%lambda, nc%N1, nc%N2, nc%e_pc, nc%e_pc_ave)
                nc%e_pair = nc%e_pc - nc%e_pc_ave
            end subroutine correct_one

        end subroutine shell_bcs_calculation

        !==========================================================================
        ! Strutinsky shell correction
        !==========================================================================

        subroutine strutinsky_correction(e, n_part, gam, p, e_sum, e_smooth, lambda, g_lambda)
            !---------------------------
            ! e(:)     : single-particle energies, ascending, one nucleon each
            ! n_part   : number of nucleons
            ! gam      : smoothing range gamma_S (MeV),  p : order of the correction polynomial (even)
            ! e_sum    : sum of the lowest n_part energies
            ! e_smooth : smooth energy  int^lambda e g(e) de,  with  int^lambda g(e) de = n_part
            ! g_lambda : smooth level density at lambda
            !---------------------------
            implicit none
            real(dp), intent(in) :: e(:), gam
            integer, intent(in) :: n_part, p
            real(dp), intent(out) :: e_sum, e_smooth, lambda, g_lambda
            real(dp) :: lo, hi, n_lo
            integer :: it

            if (n_part > size(e)) stop "strutinsky_correction: more nucleons than levels"
            e_sum = sum(e(1:n_part))

            ! bracket and bisect N_smooth(lambda) = n_part
            lo = e(1) - 10.0_dp * gam
            hi = e(n_part) + 10.0_dp * gam
            n_lo = smooth_number(e, lo, gam, p) - n_part
            if (n_lo > 0.0_dp .or. smooth_number(e, hi, gam, p) - n_part < 0.0_dp) &
                stop "strutinsky_correction: smooth Fermi level not bracketed"
            do it = 1, 200
                lambda = 0.5_dp * (lo + hi)
                if ((smooth_number(e, lambda, gam, p) - n_part) * n_lo > 0.0_dp) then
                    lo = lambda
                else
                    hi = lambda
                end if
                if (hi - lo < 1.0e-12_dp) exit
            end do
            lambda = 0.5_dp * (lo + hi)
            e_smooth = smooth_energy(e, lambda, gam, p)
            g_lambda = smooth_density(e, lambda, gam, p)
        end subroutine strutinsky_correction

        pure subroutine hermite(x, nmax, Hn)
            ! physicists' Hermite polynomials H_0 .. H_nmax at x
            implicit none
            real(dp), intent(in) :: x
            integer, intent(in) :: nmax
            real(dp), intent(out) :: Hn(0:nmax)
            integer :: n
            Hn(0) = 1.0_dp
            if (nmax >= 1) Hn(1) = 2.0_dp * x
            do n = 1, nmax - 1
                Hn(n + 1) = 2.0_dp * x * Hn(n) - 2.0_dp * n * Hn(n - 1)
            end do
        end subroutine hermite

        pure function corr_coef(m) result(c)
            ! c_m = (-1)^m / (4^m m!)
            implicit none
            integer, intent(in) :: m
            real(dp) :: c
            integer :: k
            c = 1.0_dp
            do k = 1, m
                c = -c / (4.0_dp * k)
            end do
        end function corr_coef

        pure function smooth_density(e, lambda, gam, p) result(g)
            ! g(lambda) = (1/gam) sum_i f((lambda - e_i)/gam)
            implicit none
            real(dp), intent(in) :: e(:), lambda, gam
            integer, intent(in) :: p
            real(dp) :: g, x, s, Hn(0:p)
            integer :: i, m
            g = 0.0_dp
            do i = 1, size(e)
                x = (lambda - e(i)) / gam
                if (abs(x) > 12.0_dp) cycle
                call hermite(x, p, Hn)
                s = 0.0_dp
                do m = 0, p / 2
                    s = s + corr_coef(m) * Hn(2 * m)
                end do
                g = g + exp(-x * x) * s
            end do
            g = g / (sqrt(pi) * gam)
        end function smooth_density

        pure function smooth_number(e, lambda, gam, p) result(nn)
            ! N(lambda) = sum_i [ erfc(-x)/2 - exp(-x^2)/sqrt(pi) sum_{m>=1} c_m H_{2m-1}(x) ]
            ! (int_{-inf}^x exp(-t^2) H_n(t) dt = -exp(-x^2) H_{n-1}(x), n >= 1)
            implicit none
            real(dp), intent(in) :: e(:), lambda, gam
            integer, intent(in) :: p
            real(dp) :: nn, x, s, Hn(0:p)
            integer :: i, m
            nn = 0.0_dp
            do i = 1, size(e)
                x = (lambda - e(i)) / gam
                if (x > 12.0_dp) then
                    nn = nn + 1.0_dp
                    cycle
                else if (x < -12.0_dp) then
                    cycle
                end if
                call hermite(x, p, Hn)
                s = 0.0_dp
                do m = 1, p / 2
                    s = s + corr_coef(m) * Hn(2 * m - 1)
                end do
                nn = nn + 0.5_dp * erfc(-x) - exp(-x * x) * s / sqrt(pi)
            end do
        end function smooth_number

        pure function smooth_energy(e, lambda, gam, p) result(es)
            ! E(lambda) = sum_i [ e_i n_i(x) + gam int_{-inf}^x t f(t) dt ],
            ! int_{-inf}^x t f(t) dt = -exp(-x^2)/sqrt(pi) sum_m c_m [H_2m(x)/2 + 2m H_{2m-2}(x)]
            implicit none
            real(dp), intent(in) :: e(:), lambda, gam
            integer, intent(in) :: p
            real(dp) :: es, x, s1, s2, Hn(0:p)
            integer :: i, m
            es = 0.0_dp
            do i = 1, size(e)
                x = (lambda - e(i)) / gam
                if (x > 12.0_dp) then
                    es = es + e(i)
                    cycle
                else if (x < -12.0_dp) then
                    cycle
                end if
                call hermite(x, p, Hn)
                s1 = 0.0_dp
                s2 = 0.5_dp * Hn(0)
                do m = 1, p / 2
                    s1 = s1 + corr_coef(m) * Hn(2 * m - 1)
                    s2 = s2 + corr_coef(m) * (0.5_dp * Hn(2 * m) + 2.0_dp * m * Hn(2 * m - 2))
                end do
                es = es + e(i) * (0.5_dp * erfc(-x) - exp(-x * x) * s1 / sqrt(pi)) &
                        - gam * exp(-x * x) * s2 / sqrt(pi)
            end do
        end function smooth_energy

        !==========================================================================
        ! BCS pairing correction
        !==========================================================================

        subroutine bcs_correction(e, n_part, rho, delta_G, G, delta, lambda, N1, N2, e_pc, e_pc_ave)
            !---------------------------
            ! e(:)    : single-particle energies, ascending, one nucleon each (Kramers pairs)
            ! rho     : average density of doubly degenerate levels = g_smooth(lambda_smooth)/2
            ! delta_G : effective-interaction pairing gap
            ! out: G, BCS gap delta and Fermi level lambda, window N1..N2 (doubly degenerate
            !      levels numbered from the bottom), e_pc (Moller-Nix Eq. 5), e_pc_ave (Eq. 27)
            !---------------------------
            implicit none
            real(dp), intent(in) :: e(:), rho, delta_G
            integer, intent(in) :: n_part
            real(dp), intent(out) :: G, delta, lambda, e_pc, e_pc_ave
            integer, intent(out) :: N1, N2
            real(dp), allocatable :: ek(:), v2(:), nk(:)
            real(dp) :: y1, y2, n_int, lo, hi, f_lo, Ei
            integer :: nlev, kF, k, it
            logical :: odd

            nlev = size(e) / 2
            allocate(ek(nlev))
            do k = 1, nlev
                ek(k) = 0.5_dp * (e(2 * k - 1) + e(2 * k))
            end do
            odd = mod(n_part, 2) == 1
            kF = (n_part + 1) / 2          ! last (partly) occupied level

            ! pairing window
            N2 = kF
            do k = kF, nlev
                if (ek(k) < 0.0_dp .or. ek(k) <= ek(kF) + 5.0_dp) N2 = k
            end do
            N2 = max(N2, kF + 1)
            N1 = max(1, kF - (N2 - kF) + 1)
            if (N2 > nlev) stop "bcs_correction: pairing window exceeds the basis"

            ! sharp occupation (2, 1 or 0 nucleons per level)
            allocate(nk(N1:N2), v2(N1:N2))
            do k = N1, N2
                if (2 * k <= n_part) then
                    nk(k) = 2.0_dp
                else if (2 * k - 1 == n_part) then
                    nk(k) = 1.0_dp
                else
                    nk(k) = 0.0_dp
                end if
            end do
            n_int = real(n_part - 2 * (N1 - 1), dp)   ! nucleons in the window, Eq. (94)

            !---------------------------
            ! G from the average gap equation (Eq. 107)
            !   1/G = rho/2 [asinh(y2/Delta_G) - asinh(y1/Delta_G)]
            !---------------------------
            y1 = (-0.5_dp * n_part + N1 - 1) / rho
            y2 = (-0.5_dp * n_part + N2) / rho
            G = 2.0_dp / (rho * (asinh(y2 / delta_G) - asinh(y1 / delta_G)))

            ! average pairing correlation + quasi-particle energy (Moller-Nix Eq. 27)
            e_pc_ave = 0.5_dp * rho * ((y2 - G) * (y2 - sqrt(y2**2 + delta_G**2)) &
                                     + (y1 - G) * (y1 + sqrt(y1**2 + delta_G**2))) &
                     + 0.25_dp * G * rho * delta_G * (atan(y2 / delta_G) - atan(y1 / delta_G))
            if (odd) e_pc_ave = e_pc_ave + delta_G

            !---------------------------
            ! BCS equations: gap equation solved by bisection in delta,
            ! lambda(delta) from the particle number
            !---------------------------
            lo = 1.0e-6_dp
            hi = 30.0_dp
            f_lo = gap_function(lo)
            if (f_lo >= 0.0_dp) then
                ! no non-trivial solution: sharp distribution, no correlation energy
                delta = 0.0_dp
                lambda = 0.5_dp * (ek(kF) + ek(min(kF + 1, nlev)))
                e_pc = 0.0_dp
                return
            end if
            if (gap_function(hi) < 0.0_dp) stop "bcs_correction: gap not bracketed"
            do it = 1, 200
                delta = 0.5_dp * (lo + hi)
                if (gap_function(delta) < 0.0_dp) then
                    lo = delta
                else
                    hi = delta
                end if
                if (hi - lo < 1.0e-12_dp) exit
            end do
            delta = 0.5_dp * (lo + hi)
            lambda = fermi_level(delta)
            v2 = 0.5_dp * (1.0_dp - (ek(N1:N2) - lambda) / sqrt((ek(N1:N2) - lambda)**2 + delta**2))

            ! pairing correlation + quasi-particle energy (Moller-Nix Eq. 5)
            e_pc = sum((2.0_dp * v2 - nk) * ek(N1:N2)) - delta**2 / G &
                 - 0.5_dp * G * sum(2.0_dp * v2**2 - nk)
            if (odd) then
                Ei = sqrt((ek(kF) - lambda)**2 + delta**2)
                e_pc = e_pc + Ei
            end if

        contains

            function fermi_level(d) result(lam)
                ! 2 sum v_k^2 = n_int for fixed gap d (monotone in lambda)
                real(dp), intent(in) :: d
                real(dp) :: lam, a, b
                integer :: j
                a = ek(N1) - 50.0_dp
                b = ek(N2) + 50.0_dp
                do j = 1, 200
                    lam = 0.5_dp * (a + b)
                    if (sum(1.0_dp - (ek(N1:N2) - lam) / sqrt((ek(N1:N2) - lam)**2 + d * d)) < n_int) then
                        a = lam
                    else
                        b = lam
                    end if
                    if (b - a < 1.0e-13_dp) exit
                end do
                lam = 0.5_dp * (a + b)
            end function fermi_level

            function gap_function(d) result(f)
                ! 2/G - sum 1/sqrt((e_k - lambda)^2 + d^2); increases with d, zero at the solution
                real(dp), intent(in) :: d
                real(dp) :: f, lam
                lam = fermi_level(d)
                f = 2.0_dp / G - sum(1.0_dp / sqrt((ek(N1:N2) - lam)**2 + d * d))
            end function gap_function

        end subroutine bcs_correction

        !==========================================================================
        ! relative surface area of an ellipsoid
        !==========================================================================

        function ellipsoid_relative_area(shape) result(B_s)
            ! surface area / (4 pi R^2), R^3 = R_x R_y R_z
            implicit none
            type(ellipsoid_shape), intent(in) :: shape
            real(dp) :: B_s
            integer, parameter :: nq = 64
            real(dp) :: tq(nq), wq(nq), th, ph, a, b, c, st, s
            integer :: i, j

            a = shape%semi(1)
            b = shape%semi(2)
            c = shape%semi(3)
            call gauss_legendre(nq, tq, wq)
            s = 0.0_dp
            ! theta in [0, pi/2], phi in [0, pi/2]: one octant, times 8
            do j = 1, nq
                ph = 0.25_dp * pi * (tq(j) + 1.0_dp)
                do i = 1, nq
                    th = 0.25_dp * pi * (tq(i) + 1.0_dp)
                    st = sin(th)
                    s = s + wq(i) * wq(j) * st * sqrt((b * c * st * cos(ph))**2 + (a * c * st * sin(ph))**2 &
                                                      + (a * b * cos(th))**2)
                end do
            end do
            s = 8.0_dp * s * (0.25_dp * pi)**2
            B_s = s / (4.0_dp * pi * (a * b * c)**(2.0_dp / 3.0_dp))
        end function ellipsoid_relative_area

end module shell_bcs_mod
