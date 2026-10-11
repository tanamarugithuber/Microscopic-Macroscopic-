module potential_gh_mod
    !---------------------------
    ! Folded-Yukawa potential V1, its gradient and the Coulomb potential VC
    ! for a sharp generating shape of volume V (Moller et al., ADNDT 109-110 (2016), Eqs. (81), (91))
    !
    !   Y(r) = 1/(4 pi a^3) int_V exp(-s/a) / (s/a) d^3r'    ->  V1 = -V0 * Y,  grad V1 = -V0 * grad Y
    !   C(r) = int_V d^3r' / s                               ->  VC = e^2 rho_c * C,  rho_c = Z / V
    !   s = |r' - r|,  a = a_pot
    !
    ! The volume integrals are turned into surface integrals with the divergence theorem
    ! (radial vector fields centred at r;  d = r' - r,  n dS' outward):
    !   Y(r)      =  1/(4 pi)     oint [1 - (1 + s/a) exp(-s/a)] (d . n) / s^3  dS'
    !   grad Y(r) = -1/(4 pi a^2) oint exp(-s/a) / s  n dS'
    !   C(r)      =  1/2          oint (d . n) / s  dS'
    ! No Poisson equation is solved and the boundary condition at infinity is exact.
    ! The 1/s singularity of grad Y on the surface is removed with
    !   oint n / s dS' = - oint (d . n) d / s^3 dS'
    ! which gives
    !   grad Y(r) = -1/(4 pi a^2) oint [ (exp(-s/a) - 1) / s  n - (d . n) d / s^3 ] dS'
    ! All integrands are then bounded on the surface. Panels close to the evaluation
    ! point are still subdivided adaptively, because the integrands vary rapidly there.
    !
    ! The surface is parametrised as r(u, v), u in [0, pi], v in [0, 2 pi], with
    !   N(u, v) = dr/du x dr/dv   (outward normal times area element)
    ! New shapes are added by extending surface_shape.
    !
    ! Units:
    !   Length : fm
    !   Y      : dimensionless (1 deep inside the nucleus)
    !   grad Y : 1/fm
    !   C      : fm^2
    !---------------------------
    !$ use omp_lib
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    implicit none
    private

    integer, parameter :: dp = real64

    !---------------------------
    ! Surface shapes
    !---------------------------
    type, abstract, public :: surface_shape
        contains
            procedure(surface_eval), deferred :: eval
    end type surface_shape

    abstract interface
        pure subroutine surface_eval(this, u, v, r, N)
            import :: surface_shape, dp
            class(surface_shape), intent(in) :: this
            real(dp), intent(in) :: u, v
            real(dp), intent(out) :: r(3), N(3)
        end subroutine surface_eval
    end interface

    type, extends(surface_shape), public :: ellipsoid_shape
        real(dp) :: semi(3) = 0.0_dp    ! fm, semi-axes along x, y, z
        contains
            procedure :: eval => ellipsoid_eval
    end type ellipsoid_shape

    !---------------------------
    ! Surface quadrature
    !   the (u, v) rectangle is split into n_panel_u x n_panel_v panels,
    !   each integrated with order x order Gauss-Legendre points.
    !   A panel is subdivided (up to max_depth times) while
    !   |r - panel centre| < eta * panel size.
    !---------------------------
    type, public :: surface_quadrature
        real(dp) :: a = 0.0_dp          ! fm, range of the Yukawa function
        integer :: n_panel_u = 0
        integer :: n_panel_v = 0
        integer :: order = 0
        integer :: max_depth = 8
        real(dp) :: eta = 3.0_dp
        class(surface_shape), allocatable :: shape

        real(dp), allocatable :: gl_x(:), gl_w(:)   ! Gauss-Legendre on [0, 1]
        ! precomputed nodes of the undivided panels
        real(dp), allocatable :: rn(:,:,:)          ! rn(3, order**2, panel)
        real(dp), allocatable :: Nw(:,:,:)          ! Nw(3, order**2, panel) = N * weights
        real(dp), allocatable :: pc(:,:)            ! pc(3, panel), panel centre
        real(dp), allocatable :: psize(:)           ! panel size (fm)
        real(dp), allocatable :: pbox(:,:)          ! pbox(4, panel) = u0, u1, v0, v1

        contains
            procedure :: initialize => initialize_surface_quadrature
            procedure :: at_point
            procedure :: on_grid
            procedure :: volume
    end type surface_quadrature

    public :: ellipsoid_from_eps_gamma
    public :: gauss_legendre

    contains
        pure subroutine ellipsoid_eval(this, u, v, r, N)
            !   r = (a sin u cos v, b sin u sin v, c cos u)
            !   N = (b c sin^2 u cos v, a c sin^2 u sin v, a b sin u cos u)
            implicit none
            class(ellipsoid_shape), intent(in) :: this
            real(dp), intent(in) :: u, v
            real(dp), intent(out) :: r(3), N(3)
            real(dp) :: su, cu, sv, cv

            su = sin(u)
            cu = cos(u)
            sv = sin(v)
            cv = cos(v)
            r = [this%semi(1) * su * cv, this%semi(2) * su * sv, this%semi(3) * cu]
            N = [this%semi(2) * this%semi(3) * su * su * cv, &
                 this%semi(1) * this%semi(3) * su * su * sv, &
                 this%semi(1) * this%semi(2) * su * cu]
        end subroutine ellipsoid_eval

        function ellipsoid_from_eps_gamma(R, eps, gamma) result(shape)
            !---------------------------
            ! Ellipsoid whose semi-axes are inversely proportional to the oscillator
            ! frequencies of ho_basis_mod::frequencies_from_eps_gamma
            !   R_x : R_y : R_z = 1/omega_x : 1/omega_y : 1/omega_z,   R_x R_y R_z = R^3
            ! gamma in radian.
            !---------------------------
            implicit none
            real(dp), intent(in) :: R, eps, gamma
            type(ellipsoid_shape) :: shape
            real(dp) :: f(3)

            f(1) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma + 2.0_dp * pi / 3.0_dp)
            f(2) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma - 2.0_dp * pi / 3.0_dp)
            f(3) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma)
            if (any(f <= 0.0_dp)) stop "ellipsoid_from_eps_gamma: non-positive frequency"
            shape%semi = R / f * product(f)**(1.0_dp / 3.0_dp)
        end function ellipsoid_from_eps_gamma

        subroutine gauss_legendre(n, x, w)
            ! Gauss-Legendre nodes and weights on [-1, 1] (Newton iteration on P_n)
            implicit none
            integer, intent(in) :: n
            real(dp), intent(out) :: x(n), w(n)
            real(dp) :: z, z1, p0, p1, p2, dp_n
            integer :: i, k, iter

            do i = 1, (n + 1) / 2
                z = cos(pi * (i - 0.25_dp) / (n + 0.5_dp))
                do iter = 1, 100
                    p0 = 1.0_dp
                    p1 = z
                    do k = 2, n
                        p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k
                        p0 = p1
                        p1 = p2
                    end do
                    if (n == 1) p0 = 1.0_dp
                    dp_n = n * (z * p1 - p0) / (z * z - 1.0_dp)
                    z1 = z
                    z = z1 - p1 / dp_n
                    if (abs(z - z1) < 1.0e-15_dp) exit
                end do
                x(i) = -z
                x(n + 1 - i) = z
                w(i) = 2.0_dp / ((1.0_dp - z * z) * dp_n * dp_n)
                w(n + 1 - i) = w(i)
            end do
        end subroutine gauss_legendre

        subroutine initialize_surface_quadrature(this, shape, a, n_panel_u, n_panel_v, order)
            implicit none
            class(surface_quadrature), intent(inout) :: this
            class(surface_shape), intent(in) :: shape
            real(dp), intent(in) :: a
            integer, intent(in) :: n_panel_u, n_panel_v, order
            integer :: iu, iv, ip, np

            this%a = a
            this%n_panel_u = n_panel_u
            this%n_panel_v = n_panel_v
            this%order = order
            if (allocated(this%shape)) deallocate(this%shape)
            allocate(this%shape, source=shape)

            if (allocated(this%gl_x)) deallocate(this%gl_x, this%gl_w)
            allocate(this%gl_x(order), this%gl_w(order))
            call gauss_legendre(order, this%gl_x, this%gl_w)
            this%gl_x = 0.5_dp * (this%gl_x + 1.0_dp)
            this%gl_w = 0.5_dp * this%gl_w

            np = n_panel_u * n_panel_v
            if (allocated(this%rn)) deallocate(this%rn, this%Nw, this%pc, this%psize, this%pbox)
            allocate(this%rn(3, order**2, np), this%Nw(3, order**2, np))
            allocate(this%pc(3, np), this%psize(np), this%pbox(4, np))

            do iv = 1, n_panel_v
                do iu = 1, n_panel_u
                    ip = iu + (iv - 1) * n_panel_u
                    this%pbox(:, ip) = [pi * (iu - 1) / n_panel_u, pi * iu / n_panel_u, &
                                        2.0_dp * pi * (iv - 1) / n_panel_v, 2.0_dp * pi * iv / n_panel_v]
                    call panel_nodes(this, this%pbox(:, ip), this%rn(:, :, ip), this%Nw(:, :, ip))
                    call panel_geometry(this, this%pbox(:, ip), this%pc(:, ip), this%psize(ip))
                end do
            end do
        end subroutine initialize_surface_quadrature

        pure subroutine panel_nodes(this, box, rn, Nw)
            ! Gauss points of the panel box = (u0, u1, v0, v1): positions and N * weight
            implicit none
            class(surface_quadrature), intent(in) :: this
            real(dp), intent(in) :: box(4)
            real(dp), intent(out) :: rn(:,:), Nw(:,:)
            real(dp) :: du, dv, u, v, N(3)
            integer :: i, j, m

            du = box(2) - box(1)
            dv = box(4) - box(3)
            m = 0
            do j = 1, this%order
                v = box(3) + dv * this%gl_x(j)
                do i = 1, this%order
                    u = box(1) + du * this%gl_x(i)
                    m = m + 1
                    call this%shape%eval(u, v, rn(:, m), N)
                    Nw(:, m) = N * (du * dv * this%gl_w(i) * this%gl_w(j))
                end do
            end do
        end subroutine panel_nodes

        pure subroutine panel_geometry(this, box, centre, psize)
            ! centre and size (max distance from the centre to corners and edge midpoints)
            implicit none
            class(surface_quadrature), intent(in) :: this
            real(dp), intent(in) :: box(4)
            real(dp), intent(out) :: centre(3), psize
            real(dp) :: r(3), N(3), uu(3), vv(3)
            integer :: i, j

            uu = [box(1), 0.5_dp * (box(1) + box(2)), box(2)]
            vv = [box(3), 0.5_dp * (box(3) + box(4)), box(4)]
            call this%shape%eval(uu(2), vv(2), centre, N)
            psize = 0.0_dp
            do j = 1, 3
                do i = 1, 3
                    call this%shape%eval(uu(i), vv(j), r, N)
                    psize = max(psize, norm2(r - centre))
                end do
            end do
        end subroutine panel_geometry

        pure subroutine add_nodes(x, a, rn, Nw, acc)
            !---------------------------
            ! acc(1)   += [1 - (1 + t) e^-t] (d . Nw) / s^3,   t = s/a
            ! acc(2:4) += (e^-t - 1) / s * Nw - (d . Nw) d / s^3
            ! acc(5)   += (d . Nw) / s
            !---------------------------
            implicit none
            real(dp), intent(in) :: x(3), a, rn(:,:), Nw(:,:)
            real(dp), intent(inout) :: acc(5)
            real(dp) :: d(3), s2, s, t, e, g, em1, dn
            integer :: m

            do m = 1, size(rn, 2)
                d = rn(:, m) - x
                s2 = d(1) * d(1) + d(2) * d(2) + d(3) * d(3)
                if (s2 <= 1.0e-28_dp) cycle    ! node coincides with x (measure zero)
                s = sqrt(s2)
                t = s / a
                e = exp(-t)
                if (t < 0.05_dp) then
                    ! 1 - (1+t) e^-t = t^2/2 - t^3/3 + t^4/8 - t^5/30 + t^6/144 - t^7/840 + ...
                    g = t * t * (0.5_dp + t * (-1.0_dp / 3.0_dp + t * (0.125_dp + t * (-1.0_dp / 30.0_dp &
                        + t * (1.0_dp / 144.0_dp - t / 840.0_dp)))))
                    ! e^-t - 1 = -t (1 - t/2 + t^2/6 - t^3/24 + t^4/120 - t^5/720 + ...)
                    em1 = - t * (1.0_dp + t * (-0.5_dp + t * (1.0_dp / 6.0_dp + t * (-1.0_dp / 24.0_dp &
                        + t * (1.0_dp / 120.0_dp - t / 720.0_dp)))))
                else
                    g = 1.0_dp - (1.0_dp + t) * e
                    em1 = e - 1.0_dp
                end if
                dn = d(1) * Nw(1, m) + d(2) * Nw(2, m) + d(3) * Nw(3, m)
                acc(1) = acc(1) + g * dn / (s2 * s)
                acc(2:4) = acc(2:4) + (em1 / s) * Nw(:, m) - (dn / (s2 * s)) * d
                acc(5) = acc(5) + dn / s
            end do
        end subroutine add_nodes

        pure recursive subroutine refine_panel(this, x, box, depth, acc)
            implicit none
            class(surface_quadrature), intent(in) :: this
            real(dp), intent(in) :: x(3), box(4)
            integer, intent(in) :: depth
            real(dp), intent(inout) :: acc(5)
            real(dp) :: centre(3), psize, um, vm
            real(dp) :: rn(3, this%order**2), Nw(3, this%order**2)

            call panel_geometry(this, box, centre, psize)
            if (depth < this%max_depth .and. norm2(x - centre) < this%eta * psize) then
                um = 0.5_dp * (box(1) + box(2))
                vm = 0.5_dp * (box(3) + box(4))
                call refine_panel(this, x, [box(1), um, box(3), vm], depth + 1, acc)
                call refine_panel(this, x, [um, box(2), box(3), vm], depth + 1, acc)
                call refine_panel(this, x, [box(1), um, vm, box(4)], depth + 1, acc)
                call refine_panel(this, x, [um, box(2), vm, box(4)], depth + 1, acc)
            else
                call panel_nodes(this, box, rn, Nw)
                call add_nodes(x, this%a, rn, Nw, acc)
            end if
        end subroutine refine_panel

        pure subroutine at_point(this, x, Y, gradY, C)
            ! Y, grad Y and C at the point x (fm)
            implicit none
            class(surface_quadrature), intent(in) :: this
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: Y, gradY(3), C
            real(dp) :: acc(5)
            integer :: ip

            acc = 0.0_dp
            do ip = 1, size(this%psize)
                if (norm2(x - this%pc(:, ip)) < this%eta * this%psize(ip)) then
                    call refine_panel(this, x, this%pbox(:, ip), 0, acc)
                else
                    call add_nodes(x, this%a, this%rn(:, :, ip), this%Nw(:, :, ip), acc)
                end if
            end do
            Y = acc(1) / (4.0_dp * pi)
            gradY = - acc(2:4) / (4.0_dp * pi * this%a**2)
            C = 0.5_dp * acc(5)
        end subroutine at_point

        function volume(this) result(vol)
            ! enclosed volume (1/3) oint r . n dS, to check the shape and the quadrature
            implicit none
            class(surface_quadrature), intent(in) :: this
            real(dp) :: vol
            integer :: ip, m

            vol = 0.0_dp
            do ip = 1, size(this%psize)
                do m = 1, size(this%rn, 2)
                    vol = vol + dot_product(this%rn(:, m, ip), this%Nw(:, m, ip))
                end do
            end do
            vol = vol / 3.0_dp
        end function volume

        subroutine on_grid(this, xg, yg, zg, Y, dYx, dYy, dYz, C, d2h)
            !---------------------------
            ! Y, grad Y and C on the tensor grid (xg(i), yg(j), zg(k)),
            ! e.g. the Gauss-Hermite points  xg = b_x * xi  of ho_basis_mod.
            ! d2h = .true.: the shape is reflection symmetric in x, y and z and the grid
            !   satisfies xg(n+1-i) = -xg(i); only the octant x, y, z >= 0 is computed
            !   (Y, C even; dY/dx odd in x and even in y, z, ...).
            !---------------------------
            implicit none
            class(surface_quadrature), intent(in) :: this
            real(dp), intent(in) :: xg(:), yg(:), zg(:)
            real(dp), intent(out) :: Y(size(xg), size(yg), size(zg))
            real(dp), intent(out) :: dYx(size(xg), size(yg), size(zg))
            real(dp), intent(out) :: dYy(size(xg), size(yg), size(zg))
            real(dp), intent(out) :: dYz(size(xg), size(yg), size(zg))
            real(dp), intent(out) :: C(size(xg), size(yg), size(zg))
            logical, intent(in) :: d2h
            integer :: nx, ny, nz, i0, j0, k0, i, j, k, n, ip, mi, mj, mk
            real(dp) :: g(3)

            nx = size(xg)
            ny = size(yg)
            nz = size(zg)
            if (d2h) then
                if (.not. (mirrored(xg) .and. mirrored(yg) .and. mirrored(zg))) &
                    stop "on_grid: d2h requires grids symmetric about zero"
                i0 = nx / 2 + 1
                j0 = ny / 2 + 1
                k0 = nz / 2 + 1
            else
                i0 = 1
                j0 = 1
                k0 = 1
            end if

            n = (nx - i0 + 1) * (ny - j0 + 1) * (nz - k0 + 1)
            !$omp parallel do schedule(dynamic, 16) default(none) private(ip, i, j, k, g) &
            !$omp shared(this, xg, yg, zg, Y, dYx, dYy, dYz, C, n, i0, j0, k0, nx, ny)
            do ip = 0, n - 1
                i = i0 + mod(ip, nx - i0 + 1)
                j = j0 + mod(ip / (nx - i0 + 1), ny - j0 + 1)
                k = k0 + ip / ((nx - i0 + 1) * (ny - j0 + 1))
                call this%at_point([xg(i), yg(j), zg(k)], Y(i, j, k), g, C(i, j, k))
                dYx(i, j, k) = g(1)
                dYy(i, j, k) = g(2)
                dYz(i, j, k) = g(3)
            end do
            !$omp end parallel do

            if (.not. d2h) return
            ! fill the other octants: index i <-> nx+1-i
            do k = 1, nz
                mk = merge(k, nz + 1 - k, k >= k0)
                do j = 1, ny
                    mj = merge(j, ny + 1 - j, j >= j0)
                    do i = 1, nx
                        mi = merge(i, nx + 1 - i, i >= i0)
                        if (i == mi .and. j == mj .and. k == mk) cycle
                        Y(i, j, k) = Y(mi, mj, mk)
                        C(i, j, k) = C(mi, mj, mk)
                        dYx(i, j, k) = merge(1.0_dp, -1.0_dp, i == mi) * dYx(mi, mj, mk)
                        dYy(i, j, k) = merge(1.0_dp, -1.0_dp, j == mj) * dYy(mi, mj, mk)
                        dYz(i, j, k) = merge(1.0_dp, -1.0_dp, k == mk) * dYz(mi, mj, mk)
                    end do
                end do
            end do

        contains
            pure logical function mirrored(g)
                real(dp), intent(in) :: g(:)
                integer :: m
                mirrored = .true.
                do m = 1, size(g)
                    if (abs(g(m) + g(size(g) + 1 - m)) > 1.0e-10_dp * maxval(abs(g))) mirrored = .false.
                end do
            end function mirrored
        end subroutine on_grid

end module potential_gh_mod
