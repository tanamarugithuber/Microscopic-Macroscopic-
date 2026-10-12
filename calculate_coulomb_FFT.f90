module calculate_coulomb_FFT_mod
    !---------------------------
    ! Coulomb potential of a density on a 3D grid with the FALR method
    !   G. Lauritsch and P.-G. Reinhard, Int. J. Mod. Phys. C 5 (1994) 65  (thesis/FFTcoulomb.pdf)
    !
    ! Solves  Laplacian V = -4 pi rho,  V -> 0 for |r| -> infinity, i.e.
    !   V(r) = int rho(r') / |r - r'| d^3r'
    ! (the same equation as CG_method_mod::CG_method_poisson3D, but with the exact
    !  boundary condition at infinity instead of a Dirichlet condition on the box).
    !
    ! Steps (Sec. 2.3 of the paper; l_cut = 4 here, the paper uses 2):
    !   1. multipole moments q_lm = sum rho r^l S_lm(e_r) dV,  l <= l_cut       (Eq. 5)
    !   2. long-range density rho_long = sum a_lm rho_l(r) S_lm(e_r),
    !      a_lm chosen so that rho_long has the same moments q_lm on the grid    (Eqs. 6-8)
    !   3. V_long from the radial equation for each l, solved by its Green's
    !      function on a fine radial mesh                                         (Eqs. 9, 10)
    !   4. V_short from rho_short = rho - rho_long with FFT:
    !      V_short(k) = 4 pi rho_short(k) / k^2                                   (Eqs. 11-14)
    !   5. V = V_long + V_short                                                   (Eq. 15)
    !
    ! S_lm are real spherical harmonics. r^l S_lm(e_r) are polynomials in x, y, z,
    ! so no angles (acos, atan2) are needed and r = 0 is regular.
    ! Reference radial densities (Eq. 16):
    !   rho_l(r) = 2 r^l / (R^l + r^l) * rho_0(r),  rho_0(r) = 1 / (1 + exp((r - R)/a))
    !
    ! l_cut = 4 because deformed nuclei have large hexadecapole moments; with l_cut = 2
    ! their periodic images give ~1e-5 relative error in a 20 fm box (1e-6 with l_cut = 4).
    !
    ! k = 0 component: rho_short has no l <= l_cut moments, so rho_short(k) ~ -k^2/6 int r^2 rho_short
    ! and V_short(k->0) -> -(2 pi / 3) int r^2 rho_short. This fixes the average of V_short
    ! over the box, which a plain periodic FFT solution leaves undetermined.
    !
    ! Grid: rho(i,j,k) at x = (i - (nx+1)/2) h, y = (j - (ny+1)/2) h, z = (k - (nz+1)/2) h
    !       (the centre of the box is the origin; nx, ny, nz may differ).
    ! Requires FFTW3.
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
    private

    integer, parameter :: dp = real64
    integer, parameter :: l_cut = 4
    integer, parameter :: n_lm = (l_cut + 1)**2   ! number of (l, m) pairs

    public :: calculate_coulomb

    contains

        subroutine calculate_coulomb(rho, R_pot, hx, V, a_ref)
            !---------------------------
            ! rho   : density on the grid (e.g. grid%density_index3D)
            ! R_pot : radius R of the reference densities rho_l (fm), e.g. nucleus%R_a
            ! hx    : grid spacing (fm), equal on all axes
            ! V     : int rho(r') / |r - r'| d^3r'  on the grid (fm^2 x [rho])
            ! a_ref : diffuseness a of rho_0 (fm), default 0.5 fm. The result does not
            !         depend on R_pot and a_ref except through small aliasing errors.
            !---------------------------
            implicit none
            real(dp), intent(in) :: rho(:,:,:)
            real(dp), intent(in) :: R_pot
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: V(:,:,:)
            real(dp), intent(in), optional :: a_ref

            real(dp), allocatable :: rho_long(:,:,:), rho_short(:,:,:)
            real(dp), allocatable :: V_long(:,:,:), V_short(:,:,:)
            real(dp), allocatable :: r_fine(:), rho_l(:,:), phi_l(:,:)
            real(dp) :: q_lm(n_lm), a_lm(n_lm), M(n_lm, n_lm), sh(n_lm)
            real(dp) :: x, y, z, r, a0, r_max, dr, dV, M2
            integer :: nx, ny, nz, i, j, k, ilm, jlm, n_fine
            integer :: l_of(n_lm)

            nx = size(rho, 1)
            ny = size(rho, 2)
            nz = size(rho, 3)
            if (any(shape(V) /= shape(rho))) stop "calculate_coulomb: V and rho must have the same shape"
            a0 = 0.5_dp
            if (present(a_ref)) a0 = a_ref
            dV = hx**3
            do ilm = 1, n_lm
                l_of(ilm) = int(sqrt(real(ilm - 1, dp)) + 1.0e-9_dp)
            end do

            !---------------------------
            ! reference radial densities rho_l(r)/r^l and radial potentials phi_l(r)/r^l
            ! on a fine radial mesh
            !---------------------------
            r_max = 0.5_dp * hx * sqrt(real(nx**2 + ny**2 + nz**2, dp)) + 2.0_dp * hx
            dr = hx / 64.0_dp
            n_fine = ceiling(max(r_max, R_pot + 40.0_dp * a0) / dr) + 1
            allocate(r_fine(0:n_fine), rho_l(0:n_fine, 0:l_cut), phi_l(0:n_fine, 0:l_cut))
            call radial_solutions(R_pot, a0, dr, n_fine, r_fine, rho_l, phi_l)

            !---------------------------
            ! 1. multipole moments q_lm and
            ! 2. the matrix M(lm, l'm') = sum rho_l'(r) S_l'm' r^l S_lm dV
            !    so that rho_long reproduces q_lm exactly on the grid: M a = q
            !---------------------------
            q_lm = 0.0_dp
            M = 0.0_dp
            !$omp parallel do collapse(3) default(none) private(i, j, k, x, y, z, r, sh, ilm, jlm) &
            !$omp shared(nx, ny, nz, hx, rho, rho_l, dr, n_fine, l_of) reduction(+:q_lm, M)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        call position(i, j, k, nx, ny, nz, hx, x, y, z)
                        r = sqrt(x * x + y * y + z * z)
                        call solid_harmonics(x, y, z, sh)
                        q_lm = q_lm + rho(i, j, k) * sh
                        do jlm = 1, n_lm
                            do ilm = 1, n_lm
                                M(ilm, jlm) = M(ilm, jlm) + sh(ilm) * sh(jlm) &
                                            * interpolate(rho_l(:, l_of(jlm)), dr, n_fine, r)
                            end do
                        end do
                    end do
                end do
            end do
            !$omp end parallel do
            q_lm = q_lm * dV
            M = M * dV
            call solve_small(M, q_lm, a_lm)

            !---------------------------
            ! rho_long, rho_short and V_long
            !   rho_long = sum a_lm [rho_l(r)/r^l] [r^l S_lm]
            !   V_long   = sum a_lm [phi_l(r)/r^l] [r^l S_lm]
            !---------------------------
            allocate(rho_long(nx, ny, nz), rho_short(nx, ny, nz), V_long(nx, ny, nz), V_short(nx, ny, nz))
            M2 = 0.0_dp
            !$omp parallel do collapse(3) default(none) private(i, j, k, x, y, z, r, sh, ilm) &
            !$omp shared(nx, ny, nz, hx, rho, rho_l, phi_l, dr, n_fine, l_of, a_lm, rho_long, rho_short, V_long) &
            !$omp reduction(+:M2)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        call position(i, j, k, nx, ny, nz, hx, x, y, z)
                        r = sqrt(x * x + y * y + z * z)
                        call solid_harmonics(x, y, z, sh)
                        rho_long(i, j, k) = 0.0_dp
                        V_long(i, j, k) = 0.0_dp
                        do ilm = 1, n_lm
                            rho_long(i, j, k) = rho_long(i, j, k) &
                                + a_lm(ilm) * interpolate(rho_l(:, l_of(ilm)), dr, n_fine, r) * sh(ilm)
                            V_long(i, j, k) = V_long(i, j, k) &
                                + a_lm(ilm) * interpolate(phi_l(:, l_of(ilm)), dr, n_fine, r) * sh(ilm)
                        end do
                        rho_short(i, j, k) = rho(i, j, k) - rho_long(i, j, k)
                        M2 = M2 + r * r * rho_short(i, j, k)
                    end do
                end do
            end do
            !$omp end parallel do
            M2 = M2 * dV

            !---------------------------
            ! 3. V_short by FFT
            !---------------------------
            call poisson_fft(rho_short, hx, M2, V_short)

            !---------------------------
            ! 4. V = V_long + V_short
            !---------------------------
            V = V_long + V_short

            deallocate(rho_long, rho_short, V_long, V_short, r_fine, rho_l, phi_l)
        end subroutine calculate_coulomb

        pure subroutine position(i, j, k, nx, ny, nz, hx, x, y, z)
            ! coordinates of grid point (i,j,k); the centre of the box is the origin
            implicit none
            integer, intent(in) :: i, j, k, nx, ny, nz
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: x, y, z
            x = (real(i, dp) - 0.5_dp * real(nx + 1, dp)) * hx
            y = (real(j, dp) - 0.5_dp * real(ny + 1, dp)) * hx
            z = (real(k, dp) - 0.5_dp * real(nz + 1, dp)) * hx
        end subroutine position

        pure subroutine solid_harmonics(x, y, z, sh)
            !---------------------------
            ! sh(lm) = r^l S_lm(e_r), real orthonormal spherical harmonics S_lm, l <= 4
            ! order: l = 0, 1, 2, 3, 4 and m = -l, ..., l within each l, i.e. lm = l**2 + l + m + 1
            !---------------------------
            implicit none
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: sh(n_lm)
            real(dp), parameter :: c0 = 0.28209479177387814_dp   ! 1/(2 sqrt(pi))
            real(dp), parameter :: c1 = 0.48860251190291992_dp   ! sqrt(3/(4 pi))
            real(dp), parameter :: c2 = 1.0925484305920792_dp    ! sqrt(15/pi)/2
            real(dp), parameter :: c20 = 0.31539156525252005_dp  ! sqrt(5/pi)/4
            real(dp), parameter :: c22 = 0.54627421529603959_dp  ! sqrt(15/pi)/4
            real(dp), parameter :: c33 = 0.5900435899266435_dp   ! sqrt(35/(2 pi))/4
            real(dp), parameter :: c32 = 2.890611442640554_dp    ! sqrt(105/pi)/2
            real(dp), parameter :: c31 = 0.4570457994644658_dp   ! sqrt(21/(2 pi))/4
            real(dp), parameter :: c30 = 0.3731763325901154_dp   ! sqrt(7/pi)/4
            real(dp), parameter :: c32b = 1.445305721320277_dp   ! sqrt(105/pi)/4
            real(dp), parameter :: c44 = 2.5033429417967046_dp   ! 3 sqrt(35/pi)/4
            real(dp), parameter :: c43 = 1.7701307697799304_dp   ! 3 sqrt(35/(2 pi))/4
            real(dp), parameter :: c42 = 0.9461746957575601_dp   ! 3 sqrt(5/pi)/4
            real(dp), parameter :: c41 = 0.6690465435572892_dp   ! 3 sqrt(5/(2 pi))/4
            real(dp), parameter :: c40 = 0.10578554691520431_dp  ! 3/(16 sqrt(pi))
            real(dp), parameter :: c42b = 0.47308734787878004_dp ! 3 sqrt(5/pi)/8
            real(dp), parameter :: c44b = 0.6258357354491761_dp  ! 3 sqrt(35/pi)/16
            real(dp) :: x2, y2, z2, r2

            x2 = x * x
            y2 = y * y
            z2 = z * z
            r2 = x2 + y2 + z2

            sh(1) = c0

            sh(2) = c1 * y
            sh(3) = c1 * z
            sh(4) = c1 * x

            sh(5) = c2 * x * y
            sh(6) = c2 * y * z
            sh(7) = c20 * (2.0_dp * z2 - x2 - y2)
            sh(8) = c2 * x * z
            sh(9) = c22 * (x2 - y2)

            sh(10) = c33 * y * (3.0_dp * x2 - y2)
            sh(11) = c32 * x * y * z
            sh(12) = c31 * y * (4.0_dp * z2 - x2 - y2)
            sh(13) = c30 * z * (2.0_dp * z2 - 3.0_dp * x2 - 3.0_dp * y2)
            sh(14) = c31 * x * (4.0_dp * z2 - x2 - y2)
            sh(15) = c32b * z * (x2 - y2)
            sh(16) = c33 * x * (x2 - 3.0_dp * y2)

            sh(17) = c44 * x * y * (x2 - y2)
            sh(18) = c43 * y * z * (3.0_dp * x2 - y2)
            sh(19) = c42 * x * y * (7.0_dp * z2 - r2)
            sh(20) = c41 * y * z * (7.0_dp * z2 - 3.0_dp * r2)
            sh(21) = c40 * (35.0_dp * z2 * z2 - 30.0_dp * z2 * r2 + 3.0_dp * r2 * r2)
            sh(22) = c41 * x * z * (7.0_dp * z2 - 3.0_dp * r2)
            sh(23) = c42b * (x2 - y2) * (7.0_dp * z2 - r2)
            sh(24) = c43 * x * z * (x2 - 3.0_dp * y2)
            sh(25) = c44b * (x2 * (x2 - 3.0_dp * y2) - y2 * (3.0_dp * x2 - y2))
        end subroutine solid_harmonics

        subroutine radial_solutions(R_ref, a_ref, dr, n, r, rho_l, phi_l)
            !---------------------------
            ! On r(i) = i dr:
            !   rho_l(i, l) = rho_l(r)/r^l = 2 / (R_ref^l + r^l) * rho_0(r)       (l = 0: rho_0)
            !   phi_l(i, l) = Phi_l(r)/r^l, Phi_l solving the radial equation (9) with 4 pi rho_l:
            !     Phi_l(r) = 4 pi/(2l+1) [ r^-(l+1) int_0^r s^(l+2) rho_l ds + r^l int_r^inf s^(1-l) rho_l ds ]
            !   so  phi_l = 4 pi/(2l+1) [ r^-(2l+1) A_l(r) + B_l(r) ]
            ! The integrals are accumulated with the trapezoidal rule.
            !---------------------------
            implicit none
            real(dp), intent(in) :: R_ref, a_ref, dr
            integer, intent(in) :: n
            real(dp), intent(out) :: r(0:n), rho_l(0:n, 0:l_cut), phi_l(0:n, 0:l_cut)
            real(dp), allocatable :: f(:), A_int(:), B_int(:)
            real(dp) :: rho0
            integer :: i, l

            allocate(f(0:n), A_int(0:n), B_int(0:n))
            do i = 0, n
                r(i) = i * dr
                rho0 = 1.0_dp / (1.0_dp + exp(min((r(i) - R_ref) / a_ref, 700.0_dp)))
                do l = 0, l_cut
                    if (l == 0) then
                        rho_l(i, l) = rho0
                    else
                        rho_l(i, l) = 2.0_dp / (R_ref**l + r(i)**l) * rho0
                    end if
                end do
            end do

            do l = 0, l_cut
                ! A_l(r) = int_0^r s^(2l+2) [rho_l(s)/s^l] ds
                f = r**(2 * l + 2) * rho_l(:, l)
                A_int(0) = 0.0_dp
                do i = 1, n
                    A_int(i) = A_int(i - 1) + 0.5_dp * dr * (f(i - 1) + f(i))
                end do
                ! B_l(r) = int_r^inf s [rho_l(s)/s^l] ds
                f = r * rho_l(:, l)
                B_int(n) = 0.0_dp
                do i = n - 1, 0, -1
                    B_int(i) = B_int(i + 1) + 0.5_dp * dr * (f(i) + f(i + 1))
                end do
                phi_l(0, l) = 4.0_dp * pi / (2 * l + 1) * B_int(0)
                do i = 1, n
                    phi_l(i, l) = 4.0_dp * pi / (2 * l + 1) * (A_int(i) / r(i)**(2 * l + 1) + B_int(i))
                end do
            end do
            deallocate(f, A_int, B_int)
        end subroutine radial_solutions

        pure function interpolate(f, dr, n, r) result(val)
            ! cubic Lagrange interpolation of f(0:n) tabulated on i*dr
            implicit none
            real(dp), intent(in) :: f(0:), dr, r
            integer, intent(in) :: n
            real(dp) :: val, t
            integer :: i0

            i0 = min(max(int(r / dr) - 1, 0), n - 3)
            t = r / dr - i0
            val = - f(i0)     * (t - 1.0_dp) * (t - 2.0_dp) * (t - 3.0_dp) / 6.0_dp &
                  + f(i0 + 1) * t * (t - 2.0_dp) * (t - 3.0_dp) / 2.0_dp &
                  - f(i0 + 2) * t * (t - 1.0_dp) * (t - 3.0_dp) / 2.0_dp &
                  + f(i0 + 3) * t * (t - 1.0_dp) * (t - 2.0_dp) / 6.0_dp
        end function interpolate

        subroutine solve_small(M, b, x)
            ! x = M^-1 b by Gaussian elimination with partial pivoting (n_lm x n_lm)
            implicit none
            real(dp), intent(in) :: M(n_lm, n_lm), b(n_lm)
            real(dp), intent(out) :: x(n_lm)
            real(dp) :: A(n_lm, n_lm + 1), tmp(n_lm + 1), f
            integer :: i, j, p

            A(:, 1:n_lm) = M
            A(:, n_lm + 1) = b
            do i = 1, n_lm
                p = i - 1 + maxloc(abs(A(i:n_lm, i)), 1)
                if (abs(A(p, i)) < tiny(1.0_dp)) stop "calculate_coulomb: singular moment matrix"
                if (p /= i) then
                    tmp = A(i, :)
                    A(i, :) = A(p, :)
                    A(p, :) = tmp
                end if
                do j = i + 1, n_lm
                    f = A(j, i) / A(i, i)
                    A(j, :) = A(j, :) - f * A(i, :)
                end do
            end do
            do i = n_lm, 1, -1
                x(i) = (A(i, n_lm + 1) - dot_product(A(i, i + 1:n_lm), x(i + 1:n_lm))) / A(i, i)
            end do
        end subroutine solve_small

        subroutine poisson_fft(rho_s, hx, M2, V_s)
            !---------------------------
            ! V_s(k) = 4 pi rho_s(k) / k^2  for k /= 0
            ! k = 0 : mean of V_s over the box = -(2 pi / 3) M2 / L^3,  M2 = int r^2 rho_s
            ! FFTW is row-major: a Fortran array (nx, ny, nz) is passed with dims (nz, ny, nx),
            ! and the r2c output has shape (nx/2+1, ny, nz).
            !---------------------------
            implicit none
            real(dp), intent(in) :: rho_s(:,:,:), hx, M2
            real(dp), intent(out) :: V_s(:,:,:)
            real(C_DOUBLE), allocatable :: in(:,:,:)
            complex(C_DOUBLE_COMPLEX), allocatable :: out(:,:,:)
            type(C_PTR) :: plan_f, plan_b
            real(dp) :: kx, ky, kz, k2, vol
            integer :: nx, ny, nz, i, j, k

            nx = size(rho_s, 1)
            ny = size(rho_s, 2)
            nz = size(rho_s, 3)
            allocate(in(nx, ny, nz), out(nx / 2 + 1, ny, nz))
            plan_f = fftw_plan_dft_r2c_3d(nz, ny, nx, in, out, FFTW_ESTIMATE)
            plan_b = fftw_plan_dft_c2r_3d(nz, ny, nx, out, in, FFTW_ESTIMATE)

            in = rho_s
            call fftw_execute_dft_r2c(plan_f, in, out)

            vol = real(nx, dp) * real(ny, dp) * real(nz, dp) * hx**3
            do k = 1, nz
                kz = wave_number(k, nz, hx)
                do j = 1, ny
                    ky = wave_number(j, ny, hx)
                    do i = 1, nx / 2 + 1
                        kx = 2.0_dp * pi * real(i - 1, dp) / (real(nx, dp) * hx)
                        k2 = kx * kx + ky * ky + kz * kz
                        if (i == 1 .and. j == 1 .and. k == 1) then
                            ! box average of V_s, times the number of points (unnormalised FFT)
                            out(i, j, k) = cmplx(-2.0_dp * pi / 3.0_dp * M2 / vol &
                                                 * real(nx, dp) * real(ny, dp) * real(nz, dp), 0.0_dp, kind=dp)
                        else
                            out(i, j, k) = out(i, j, k) * (4.0_dp * pi / k2)
                        end if
                    end do
                end do
            end do

            call fftw_execute_dft_c2r(plan_b, out, in)
            V_s = in / (real(nx, dp) * real(ny, dp) * real(nz, dp))

            call fftw_destroy_plan(plan_f)
            call fftw_destroy_plan(plan_b)
            deallocate(in, out)
        end subroutine poisson_fft

        pure function wave_number(i, n, hx) result(kv)
            ! FFT frequency of index i (1-based) for n points with spacing hx
            implicit none
            integer, intent(in) :: i, n
            real(dp), intent(in) :: hx
            real(dp) :: kv
            if (i - 1 <= n / 2) then
                kv = 2.0_dp * pi * real(i - 1, dp) / (real(n, dp) * hx)
            else
                kv = 2.0_dp * pi * real(i - 1 - n, dp) / (real(n, dp) * hx)
            end if
        end function wave_number

end module calculate_coulomb_FFT_mod
