module shape_mod
    !---------------------------
    ! Sharp generating shapes of the nucleus
    !
    !   ellipsoid_shape : triaxial ellipsoid with semi-axes R_x, R_y, R_z
    !   epsilon_shape   : Nilsson perturbed-spheroid (epsilon) parameterization with
    !                     eps2, gamma, eps4, eps6 (Moller et al., ADNDT 109-110 (2016), Eqs. (17)-(25))
    !
    ! epsilon_shape is defined as an equipotential of the Nilsson potential (Eq. 21 without
    ! l.s and l^2) in the stretched coordinates
    !   xi_k = sqrt(f_k) x_k,  f = [1 - 2/3 eps2 cos(gamma + 2pi/3), 1 - 2/3 eps2 cos(gamma - 2pi/3),
    !                               1 - 2/3 eps2 cos(gamma)]
    ! (f_k = omega_k / omega0 of ho_basis_mod::frequencies_from_eps_gamma before volume rescaling):
    !   rho_t^2 F(theta_t, phi_t) = C^2,
    !   F = f_x n_x^2 + f_y n_y^2 + f_z n_z^2 + 2 eps4 V4(u, v) + 2 eps6 P6(u),
    ! where n = (sin theta_t cos phi_t, sin theta_t sin phi_t, cos theta_t), u = cos theta_t,
    ! v = cos 2 phi_t. The first term equals 1 - 2/3 eps2 cos(gamma) P2(u)
    ! + eps2 sin(gamma) sqrt(8 pi/45) (Y22 + Y2-2), as in Eq. (21). C is fixed by the volume
    ! 4/3 pi R^3 (Eq. 25). With eps4 = eps6 = 0 the shape is the ellipsoid with semi-axes
    ! proportional to 1/f_k, i.e. the same as ellipsoid_from_eps_gamma.
    !   V4(u, v) = a40 P4(u) + sqrt(4 pi/9) [a42 SY42 + a44 SY44]                  (Eqs. 22, 23)
    !   a40 = (5 cos^2 gamma + 1)/6, a42 = -sqrt(30) sin(2 gamma)/12, a44 = sqrt(70) sin^2(gamma)/12
    !   SY42 = sqrt(45/(32 pi)) (1 - u^2)(7 u^2 - 1) v,  SY44 = sqrt(315/(128 pi)) (1 - u^2)^2 (2 v^2 - 1)
    ! eps4 = eps6 = 0 and the D2h symmetry (reflections in x, y, z) hold for all these shapes.
    !
    ! Every shape provides
    !   eval(u, v, r, N) : surface point r(u, v), u in [0, pi], v in [0, 2 pi], and
    !                      N = dr/du x dr/dv (outward normal times area element)
    !   inside(x)        : .true. if the point x is inside the surface
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    implicit none
    private

    integer, parameter :: dp = real64

    type, abstract, public :: surface_shape
        contains
            procedure(surface_eval), deferred :: eval
            procedure(surface_inside), deferred :: inside
            procedure :: extent => shape_extent
            procedure :: relative_area => shape_relative_area
            procedure :: volume => shape_volume
    end type surface_shape

    abstract interface
        pure subroutine surface_eval(this, u, v, r, N)
            import :: surface_shape, dp
            class(surface_shape), intent(in) :: this
            real(dp), intent(in) :: u, v
            real(dp), intent(out) :: r(3), N(3)
        end subroutine surface_eval

        pure function surface_inside(this, x) result(yes)
            import :: surface_shape, dp
            class(surface_shape), intent(in) :: this
            real(dp), intent(in) :: x(3)
            logical :: yes
        end function surface_inside
    end interface

    type, extends(surface_shape), public :: ellipsoid_shape
        real(dp) :: semi(3) = 0.0_dp    ! fm, semi-axes along x, y, z
        contains
            procedure :: eval => ellipsoid_eval
            procedure :: inside => ellipsoid_inside
    end type ellipsoid_shape

    type, extends(surface_shape), public :: epsilon_shape
        real(dp) :: eps2 = 0.0_dp
        real(dp) :: gamma = 0.0_dp      ! radian
        real(dp) :: eps4 = 0.0_dp
        real(dp) :: eps6 = 0.0_dp
        real(dp) :: f(3) = 1.0_dp       ! stretching factors f_x, f_y, f_z
        real(dp) :: C = 1.0_dp          ! fm, size of the equipotential in stretched coordinates
        contains
            procedure :: eval => epsilon_eval
            procedure :: inside => epsilon_inside
            procedure :: F_value => epsilon_F
    end type epsilon_shape

    public :: ellipsoid_from_eps_gamma
    public :: epsilon_shape_from
    public :: gauss_legendre

    contains

        !==========================================================================
        ! ellipsoid
        !==========================================================================

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

            f = stretch_factors(eps, gamma)
            shape%semi = R / f * product(f)**(1.0_dp / 3.0_dp)
        end function ellipsoid_from_eps_gamma

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

        pure function ellipsoid_inside(this, x) result(yes)
            implicit none
            class(ellipsoid_shape), intent(in) :: this
            real(dp), intent(in) :: x(3)
            logical :: yes
            yes = sum((x / this%semi)**2) <= 1.0_dp
        end function ellipsoid_inside

        !==========================================================================
        ! epsilon parameterization
        !==========================================================================

        pure function stretch_factors(eps, gamma) result(f)
            implicit none
            real(dp), intent(in) :: eps, gamma
            real(dp) :: f(3)
            f(1) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma + 2.0_dp * pi / 3.0_dp)
            f(2) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma - 2.0_dp * pi / 3.0_dp)
            f(3) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma)
            if (any(f <= 0.0_dp)) error stop "shape_mod: eps2 too large (non-positive stretching factor)"
        end function stretch_factors

        function epsilon_shape_from(R, eps2, gamma, eps4, eps6) result(shape)
            !---------------------------
            ! epsilon shape of volume 4/3 pi R^3; gamma in radian
            !   V = det(D)/3 int F^(-3/2) dOmega C^3,  D = diag(f^(-1/2))      (Eq. 25)
            !---------------------------
            implicit none
            real(dp), intent(in) :: R, eps2, gamma, eps4, eps6
            type(epsilon_shape) :: shape
            integer, parameter :: nq = 64
            real(dp) :: tq(nq), wq(nq), th, ph, s, Fv
            integer :: i, j

            shape%eps2 = eps2
            shape%gamma = gamma
            shape%eps4 = eps4
            shape%eps6 = eps6
            shape%f = stretch_factors(eps2, gamma)
            shape%C = 1.0_dp
            ! one octant (D2h), theta and phi in [0, pi/2]
            call gauss_legendre(nq, tq, wq)
            s = 0.0_dp
            do j = 1, nq
                ph = 0.25_dp * pi * (tq(j) + 1.0_dp)
                do i = 1, nq
                    th = 0.25_dp * pi * (tq(i) + 1.0_dp)
                    Fv = shape%F_value(th, ph)
                    if (Fv <= 0.0_dp) error stop "epsilon_shape_from: F <= 0, deformation too large"
                    s = s + wq(i) * wq(j) * sin(th) * Fv**(-1.5_dp)
                end do
            end do
            s = 8.0_dp * s * (0.25_dp * pi)**2
            shape%C = (4.0_dp * pi * R**3 * sqrt(product(shape%f)) / s)**(1.0_dp / 3.0_dp)
        end function epsilon_shape_from

        pure function epsilon_F(this, th, ph) result(Fv)
            ! F(theta_t, phi_t) of the module header
            implicit none
            class(epsilon_shape), intent(in) :: this
            real(dp), intent(in) :: th, ph
            real(dp) :: Fv, n(3), u, v, u2, w, P4, P6, a40, a42, a44, SY42, SY44, V4

            n = [sin(th) * cos(ph), sin(th) * sin(ph), cos(th)]
            Fv = sum(this%f * n * n)
            if (this%eps4 /= 0.0_dp .or. this%eps6 /= 0.0_dp) then
                u = n(3)
                u2 = u * u
                v = cos(2.0_dp * ph)
                w = 1.0_dp - u2
                P4 = (35.0_dp * u2 * u2 - 30.0_dp * u2 + 3.0_dp) / 8.0_dp
                P6 = (231.0_dp * u2**3 - 315.0_dp * u2 * u2 + 105.0_dp * u2 - 5.0_dp) / 16.0_dp
                a40 = (5.0_dp * cos(this%gamma)**2 + 1.0_dp) / 6.0_dp
                a42 = -sqrt(30.0_dp) * sin(2.0_dp * this%gamma) / 12.0_dp
                a44 = sqrt(70.0_dp) * sin(this%gamma)**2 / 12.0_dp
                SY42 = sqrt(45.0_dp / (32.0_dp * pi)) * w * (7.0_dp * u2 - 1.0_dp) * v
                SY44 = sqrt(315.0_dp / (128.0_dp * pi)) * w * w * (2.0_dp * v * v - 1.0_dp)
                V4 = a40 * P4 + sqrt(4.0_dp * pi / 9.0_dp) * (a42 * SY42 + a44 * SY44)
                Fv = Fv + 2.0_dp * this%eps4 * V4 + 2.0_dp * this%eps6 * P6
            end if
        end function epsilon_F

        pure subroutine epsilon_eval(this, u, v, r, N)
            !---------------------------
            ! u = theta_t, v = phi_t (stretched angles).  r = D rho e,  rho = C F^(-1/2), D = diag(f^(-1/2))
            ! N = dr/du x dr/dv = det(D) D^(-1) [rho^2 sin u e - rho rho_u sin u e_theta - rho rho_v e_phi]
            ! (derivatives of F by central differences)
            !---------------------------
            implicit none
            class(epsilon_shape), intent(in) :: this
            real(dp), intent(in) :: u, v
            real(dp), intent(out) :: r(3), N(3)
            real(dp), parameter :: hd = 1.0e-5_dp
            real(dp) :: Fv, Fu, Fw, rho, rho_u, rho_v, e(3), et(3), ep(3), su, cu, sv, cv, dinv(3)

            su = sin(u)
            cu = cos(u)
            sv = sin(v)
            cv = cos(v)
            e = [su * cv, su * sv, cu]
            et = [cu * cv, cu * sv, -su]
            ep = [-sv, cv, 0.0_dp]
            Fv = this%F_value(u, v)
            Fu = (this%F_value(u + hd, v) - this%F_value(u - hd, v)) / (2.0_dp * hd)
            Fw = (this%F_value(u, v + hd) - this%F_value(u, v - hd)) / (2.0_dp * hd)
            rho = this%C / sqrt(Fv)
            rho_u = -0.5_dp * rho * Fu / Fv
            rho_v = -0.5_dp * rho * Fw / Fv
            dinv = sqrt(this%f)
            r = rho * e / dinv
            N = (rho * rho * su * e - rho * rho_u * su * et - rho * rho_v * ep) * dinv / product(dinv)
        end subroutine epsilon_eval

        pure function epsilon_inside(this, x) result(yes)
            ! rho_t^2 F(theta_t, phi_t) <= C^2 in stretched coordinates
            implicit none
            class(epsilon_shape), intent(in) :: this
            real(dp), intent(in) :: x(3)
            logical :: yes
            real(dp) :: xi(3), rho, th, ph

            xi = sqrt(this%f) * x
            rho = norm2(xi)
            if (rho == 0.0_dp) then
                yes = .true.
                return
            end if
            th = acos(max(-1.0_dp, min(1.0_dp, xi(3) / rho)))
            ph = atan2(xi(2), xi(1))
            yes = rho * rho * this%F_value(th, ph) <= this%C**2
        end function epsilon_inside

        !==========================================================================
        ! shape-independent quantities
        !==========================================================================

        function shape_extent(this) result(ext)
            ! largest |x|, |y|, |z| on the surface (sampled)
            implicit none
            class(surface_shape), intent(in) :: this
            real(dp) :: ext(3), r(3), N(3)
            integer, parameter :: nu = 64, nv = 128
            integer :: i, j
            ext = 0.0_dp
            do j = 0, nv
                do i = 0, nu
                    call this%eval(pi * i / nu, 2.0_dp * pi * j / nv, r, N)
                    ext = max(ext, abs(r))
                end do
            end do
        end function shape_extent

        function shape_volume(this) result(vol)
            ! V = 1/3 oint r . N du dv
            implicit none
            class(surface_shape), intent(in) :: this
            real(dp) :: vol
            vol = surface_integral(this, 1)
        end function shape_volume

        function shape_relative_area(this) result(B_s)
            ! surface area / (4 pi R^2), 4/3 pi R^3 = volume                     (Eq. 70)
            implicit none
            class(surface_shape), intent(in) :: this
            real(dp) :: B_s, vol
            vol = surface_integral(this, 1)
            B_s = surface_integral(this, 2) / (4.0_dp * pi * (3.0_dp * vol / (4.0_dp * pi))**(2.0_dp / 3.0_dp))
        end function shape_relative_area

        function surface_integral(this, what) result(s)
            ! what = 1: volume, 2: area.  Gauss-Legendre over one octant times 8 (D2h)
            implicit none
            class(surface_shape), intent(in) :: this
            integer, intent(in) :: what
            real(dp) :: s, r(3), N(3), u, v
            integer, parameter :: nq = 64
            real(dp) :: tq(nq), wq(nq)
            integer :: i, j
            call gauss_legendre(nq, tq, wq)
            s = 0.0_dp
            do j = 1, nq
                v = 0.25_dp * pi * (tq(j) + 1.0_dp)
                do i = 1, nq
                    u = 0.25_dp * pi * (tq(i) + 1.0_dp)
                    call this%eval(u, v, r, N)
                    if (what == 1) then
                        s = s + wq(i) * wq(j) * dot_product(r, N) / 3.0_dp
                    else
                        s = s + wq(i) * wq(j) * norm2(N)
                    end if
                end do
            end do
            s = 8.0_dp * s * (0.25_dp * pi)**2
        end function surface_integral

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

end module shape_mod
