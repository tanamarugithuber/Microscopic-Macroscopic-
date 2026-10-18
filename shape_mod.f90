module shape_mod
    !---------------------------
    ! Sharp generating shapes of the nucleus
    !
    !   ellipsoid_shape : triaxial ellipsoid with semi-axes R_x, R_y, R_z
    !   epsilon_shape   : Nilsson perturbed-spheroid (epsilon) parameterization with
    !                     eps2, gamma, eps3, eps4, eps6 (Moller et al., ADNDT 109-110 (2016), Eqs. (17)-(25));
    !                     eps1 is fixed by the centre-of-mass condition (Notes added in proof, p. 35)
    !
    ! epsilon_shape is defined as an equipotential of the Nilsson potential (Eq. 21 without
    ! l.s and l^2) in the stretched coordinates
    !   xi_k = sqrt(f_k) x_k,  f = [1 - 2/3 eps2 cos(gamma + 2pi/3), 1 - 2/3 eps2 cos(gamma - 2pi/3),
    !                               1 - 2/3 eps2 cos(gamma)]
    ! (f_k = omega_k / omega0 of ho_basis_mod::frequencies_from_eps_gamma before volume rescaling):
    !   rho_t^2 F(theta_t, phi_t) = C^2,
    !   F = f_x n_x^2 + f_y n_y^2 + f_z n_z^2 + 2 eps1 P1(u) + 2 eps3 P3(u) + 2 eps4 V4(u, v) + 2 eps6 P6(u),
    ! where n = (sin theta_t cos phi_t, sin theta_t sin phi_t, cos theta_t), u = cos theta_t,
    ! v = cos 2 phi_t. The first term equals 1 - 2/3 eps2 cos(gamma) P2(u)
    ! + eps2 sin(gamma) sqrt(8 pi/45) (Y22 + Y2-2), as in Eq. (21). C is fixed by the volume
    ! 4/3 pi R^3 (Eq. 25). With eps1 = eps3 = eps4 = eps6 = 0 the shape is the ellipsoid with semi-axes
    ! proportional to 1/f_k, i.e. the same as ellipsoid_from_eps_gamma.
    !   V4(u, v) = a40 P4(u) + sqrt(4 pi/9) [a42 SY42 + a44 SY44]                  (Eqs. 22, 23)
    !   a40 = (5 cos^2 gamma + 1)/6, a42 = -sqrt(30) sin(2 gamma)/12, a44 = sqrt(70) sin^2(gamma)/12
    !   SY42 = sqrt(45/(32 pi)) (1 - u^2)(7 u^2 - 1) v,  SY44 = sqrt(315/(128 pi)) (1 - u^2)^2 (2 v^2 - 1)
    ! All these shapes are reflection symmetric in x and y. eps1 and eps3 break the z reflection, so
    ! the volume and surface integrals run over the quarter x, y >= 0 (theta in [0, pi],
    ! phi in [0, pi/2]) times 4.
    ! eps1 makes int z dV = 0. With rho_max = C F^(-1/2) in the stretched coordinates,
    !   int z dV = det(D) f_z^(-1/2) C^4/4 int cos(theta_t) F^(-2) dOmega,
    ! so eps1 solves M(eps1) = int cos(theta_t) F^(-2) dOmega = 0, independent of C and R.
    ! dM/d eps1 = -4 int cos^2(theta_t) F^(-3) dOmega < 0, so the root is unique (bisection).
    !
    ! Every shape provides
    !   eval(u, v, r, N) : surface point r(u, v), u in [0, pi], v in [0, 2 pi], and
    !                      N = dr/du x dr/dv (outward normal times area element)
    !   inside(x)        : .true. if the point x is inside the surface
    !---------------------------
    ! 【日本語での説明】
    ! 原子核の鋭い表面をもつ生成形状 (sharp generating shape) を表すモジュール.
    ! FRLDM の B_1, B_3, B_s (frldm_mod) と一粒子ポテンシャル (potential_gh_mod)
    ! は, すべてここで定義した表面点 r と法線 N を使って計算する.
    ! 【物理量】
    !   R          : 同じ体積をもつ球の半径 [fm] (体積 = 4/3 pi R^3)
    !   eps2       : 四重極変形 epsilon_2 (> 0 で長球, < 0 で扁平)
    !   gamma      : 三軸性の角度 [rad] (0 で z 軸まわりに軸対称)
    !   eps3       : 八重極変形 epsilon_3 (z 方向の反転対称性を破る)
    !   eps1       : 重心を原点に置くための epsilon_1 (eps3 から自動で決まる, FRDM2012 p. 35)
    !   eps4, eps6 : 十六重極, 六十四重極変形 epsilon_4, epsilon_6
    !   f(k)       : 伸長係数 f_k = 1 - 2/3 eps2 cos(gamma + 2 pi k/3)
    !                (k = 1, 2, 3 が x, y, z; 調和振動子の振動数比 omega_k / omega0)
    !   C          : 伸長座標での等ポテンシャル面の大きさ [fm] (体積保存で決まる)
    ! 【形状の式】伸長座標 xi_k = sqrt(f_k) x_k の極座標 (rho_t, theta_t, phi_t) で
    !   rho_t^2 F(theta_t, phi_t) = C^2,
    !   F = f_x n_x^2 + f_y n_y^2 + f_z n_z^2 + 2 eps1 P1 + 2 eps3 P3 + 2 eps4 V4 + 2 eps6 P6   (FRDM2012 式 21-23)
    ! 【対称性】どの形状も x, y の各面について鏡映対称.
    ! eps1, eps3 /= 0 では z の鏡映対称性が破れるので,
    ! 体積・表面積分は x, y >= 0 の 1/4 領域 (theta in [0, pi], phi in [0, pi/2]) で行い 4 倍する.
    ! 【eps1 の決め方】重心 int z dV = det(D) f_z^(-1/2) C^4/4 int cos(theta_t) F^(-2) dOmega が 0 になるように,
    ! M(eps1) = int cos(theta_t) F^(-2) dOmega = 0 を二分法で解く (C, R に依らない).
    ! dM/d eps1 = -4 int cos^2(theta_t) F^(-3) dOmega < 0 なので解は 1 つ.
    ! 【構造】抽象型 surface_shape が eval (表面点と法線) と inside (内外判定) を要求し,
    ! ellipsoid_shape と epsilon_shape がそれを実装する. volume, relative_area, extent
    ! は eval だけを使うので, どの形状にも共通に使える.
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    implicit none
    private

    integer, parameter :: dp = real64

    ! すべての形状の親 (抽象型). 子の型は eval と inside を必ず実装する.
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
        real(dp) :: eps1 = 0.0_dp
        real(dp) :: eps2 = 0.0_dp
        real(dp) :: gamma = 0.0_dp      ! radian
        real(dp) :: eps3 = 0.0_dp
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
            ! 【目的】(eps, gamma) から体積 4/3 pi R^3 の三軸楕円体を作る.
            ! 【式】半軸 R_k = R (f_x f_y f_z)^(1/3) / f_k  ->  R_x R_y R_z = R^3
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
            ! 【目的】楕円体表面のパラメータ表示. u は極角, v は方位角.
            ! 【出力】r : 表面点 [fm],  N = dr/du x dr/dv : 外向き法線 x 面積要素 [fm^2]
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
            ! 【目的】点 x [fm] が楕円体の内側かを判定する.  式: sum_k (x_k / R_k)^2 <= 1
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
            ! 【目的】伸長係数 f_k を返す (FRDM2012 式 21 の第 1 項の係数).
            !   f_x = 1 - 2/3 eps cos(gamma + 2 pi/3)
            !   f_y = 1 - 2/3 eps cos(gamma - 2 pi/3)
            !   f_z = 1 - 2/3 eps cos(gamma)
            ! f_k <= 0 では形状が定義できないので停止する.
            implicit none
            real(dp), intent(in) :: eps, gamma
            real(dp) :: f(3)
            f(1) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma + 2.0_dp * pi / 3.0_dp)
            f(2) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma - 2.0_dp * pi / 3.0_dp)
            f(3) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma)
            if (any(f <= 0.0_dp)) error stop "shape_mod: eps2 too large (non-positive stretching factor)"
        end function stretch_factors

        function epsilon_shape_from(R, eps2, gamma, eps3, eps4, eps6) result(shape)
            !---------------------------
            ! epsilon shape of volume 4/3 pi R^3; gamma in radian
            !   V = det(D)/3 int F^(-3/2) dOmega C^3,  D = diag(f^(-1/2))      (Eq. 25)
            ! eps1 from int cos(theta_t) F^(-2) dOmega = 0 (centre of mass at the origin) when eps3 /= 0
            !---------------------------
            ! 【目的】体積が 4/3 pi R^3 で重心が原点にある epsilon 形状を作る (FRDM2012 式 25, p. 35).
            ! 【式】伸長座標での体積 V_xi = C^3/3 int F^(-3/2) dOmega,
            !   実空間の体積 V = det(D) V_xi,  D = diag(f^(-1/2)),  det(D) = (f_x f_y f_z)^(-1/2)
            !   -> C = [4 pi R^3 sqrt(f_x f_y f_z) / int F^(-3/2) dOmega]^(1/3)
            ! 【eps1】eps3 /= 0 のとき, M(eps1) = int cos(theta_t) F^(-2) dOmega = 0 を二分法で解く.
            !   M は eps1 について単調減少なので, [-0.02, 0.02] から幅を倍々に広げて
            !   符号の変わる区間を探す.
            ! 積分は 1/4 領域 (theta in [0, pi], phi in [0, pi/2]) を 128 x 64 点の Gauss-Legendre で行い 4 倍する
            ! (x, y の鏡映対称性). ヤコビアンは (pi/2) (pi/4).
            implicit none
            real(dp), intent(in) :: R, eps2, gamma, eps3, eps4, eps6
            type(epsilon_shape) :: shape
            integer, parameter :: nq = 64
            real(dp) :: tq(nq), wq(nq), tt(2 * nq), wt(2 * nq), th(2 * nq), ph(nq)
            real(dp) :: s, Fv, lo, hi, mid, m_lo, m_hi, m_mid
            integer :: i, j, iter

            shape%eps1 = 0.0_dp
            shape%eps2 = eps2
            shape%gamma = gamma
            shape%eps3 = eps3
            shape%eps4 = eps4
            shape%eps6 = eps6
            shape%f = stretch_factors(eps2, gamma)
            shape%C = 1.0_dp
            ! quarter x, y >= 0: theta in [0, pi], phi in [0, pi/2]
            call gauss_legendre(nq, tq, wq)
            call gauss_legendre(2 * nq, tt, wt)
            th = 0.5_dp * pi * (tt + 1.0_dp)
            ph = 0.25_dp * pi * (tq + 1.0_dp)

            ! eps1: centre of mass at the origin (bisection on the monotonically decreasing M)
            if (eps3 /= 0.0_dp) then
                lo = -0.02_dp
                hi = 0.02_dp
                m_lo = moment(lo)
                m_hi = moment(hi)
                do iter = 1, 20
                    if (m_lo >= 0.0_dp .and. m_hi <= 0.0_dp) exit
                    if (m_lo < 0.0_dp) then
                        lo = 2.0_dp * lo
                        m_lo = moment(lo)
                    else
                        hi = 2.0_dp * hi
                        m_hi = moment(hi)
                    end if
                end do
                if (m_lo < 0.0_dp .or. m_hi > 0.0_dp) error stop "epsilon_shape_from: no root for eps1"
                do iter = 1, 100
                    mid = 0.5_dp * (lo + hi)
                    m_mid = moment(mid)
                    if (m_mid > 0.0_dp) then
                        lo = mid
                    else
                        hi = mid
                    end if
                    if (hi - lo < 1.0e-14_dp) exit
                end do
                shape%eps1 = 0.5_dp * (lo + hi)
            end if

            s = 0.0_dp
            do j = 1, nq
                do i = 1, 2 * nq
                    Fv = shape%F_value(th(i), ph(j))
                    if (Fv <= 0.0_dp) error stop "epsilon_shape_from: F <= 0, deformation too large"
                    s = s + wt(i) * wq(j) * sin(th(i)) * Fv**(-1.5_dp)
                end do
            end do
            s = 4.0_dp * s * (0.5_dp * pi) * (0.25_dp * pi)
            shape%C = (4.0_dp * pi * R**3 * sqrt(product(shape%f)) / s)**(1.0_dp / 3.0_dp)

            contains

                function moment(e1) result(m)
                    ! M(e1) = int cos(theta_t) F^(-2) dOmega / 4 over the quarter (common factors dropped)
                    ! 【目的】eps1 = e1 としたときの重心の積分 M (定数倍は省略).
                    real(dp), intent(in) :: e1
                    real(dp) :: m, Fm
                    type(epsilon_shape) :: trial
                    integer :: ii, jj
                    trial = shape
                    trial%eps1 = e1
                    m = 0.0_dp
                    do jj = 1, nq
                        do ii = 1, 2 * nq
                            Fm = trial%F_value(th(ii), ph(jj))
                            if (Fm <= 0.0_dp) error stop "epsilon_shape_from: F <= 0 while solving for eps1"
                            m = m + wt(ii) * wq(jj) * cos(th(ii)) * sin(th(ii)) / (Fm * Fm)
                        end do
                    end do
                end function moment
        end function epsilon_shape_from

        pure function epsilon_F(this, th, ph) result(Fv)
            ! F(theta_t, phi_t) of the module header
            ! 【目的】形状関数 F(theta_t, phi_t) の値 (FRDM2012 式 21-23).
            ! 【物理量】n : 伸長座標での単位方向ベクトル,  u = cos(theta_t),  v = cos(2 phi_t)
            !   P4, P6           : ルジャンドル多項式
            !   a40, a42, a44    : gamma で決まる epsilon_4 項の係数 (式 23)
            !   SY42, SY44       : 非軸対称な球面調和関数の組 (gamma /= 0 のときに効く)
            !   V4 = a40 P4 + sqrt(4 pi/9) (a42 SY42 + a44 SY44)   (式 22)
            !   P1, P3           : eps1, eps3 の項 (z の向きを区別する奇関数)
            ! 【式】F = sum_k f_k n_k^2 + 2 eps1 P1 + 2 eps3 P3 + 2 eps4 V4 + 2 eps6 P6
            !   (eps1 = eps3 = eps4 = eps6 = 0 なら第 1 項だけ)
            implicit none
            class(epsilon_shape), intent(in) :: this
            real(dp), intent(in) :: th, ph
            real(dp) :: Fv, n(3), u, v, u2, w, P1, P3, P4, P6, a40, a42, a44, SY42, SY44, V4

            n = [sin(th) * cos(ph), sin(th) * sin(ph), cos(th)]
            Fv = sum(this%f * n * n)
            if (this%eps1 /= 0.0_dp .or. this%eps3 /= 0.0_dp .or. this%eps4 /= 0.0_dp .or. this%eps6 /= 0.0_dp) then
                u = n(3)
                u2 = u * u
                v = cos(2.0_dp * ph)
                w = 1.0_dp - u2
                P1 = u
                P3 = 0.5_dp * (5.0_dp * u2 - 3.0_dp) * u
                P4 = (35.0_dp * u2 * u2 - 30.0_dp * u2 + 3.0_dp) / 8.0_dp
                P6 = (231.0_dp * u2**3 - 315.0_dp * u2 * u2 + 105.0_dp * u2 - 5.0_dp) / 16.0_dp
                a40 = (5.0_dp * cos(this%gamma)**2 + 1.0_dp) / 6.0_dp
                a42 = -sqrt(30.0_dp) * sin(2.0_dp * this%gamma) / 12.0_dp
                a44 = sqrt(70.0_dp) * sin(this%gamma)**2 / 12.0_dp
                SY42 = sqrt(45.0_dp / (32.0_dp * pi)) * w * (7.0_dp * u2 - 1.0_dp) * v
                SY44 = sqrt(315.0_dp / (128.0_dp * pi)) * w * w * (2.0_dp * v * v - 1.0_dp)
                V4 = a40 * P4 + sqrt(4.0_dp * pi / 9.0_dp) * (a42 * SY42 + a44 * SY44)
                Fv = Fv + 2.0_dp * this%eps1 * P1 + 2.0_dp * this%eps3 * P3 &
                           + 2.0_dp * this%eps4 * V4 + 2.0_dp * this%eps6 * P6
            end if
        end function epsilon_F

        pure subroutine epsilon_eval(this, u, v, r, N)
            !---------------------------
            ! u = theta_t, v = phi_t (stretched angles).  r = D rho e,  rho = C F^(-1/2), D = diag(f^(-1/2))
            ! N = dr/du x dr/dv = det(D) D^(-1) [rho^2 sin u e - rho rho_u sin u e_theta - rho rho_v e_phi]
            ! (derivatives of F by central differences)
            !---------------------------
            ! 【目的】epsilon 形状表面のパラメータ表示.
            ! u, v は伸長座標の角度 theta_t, phi_t であり, 実空間の極角・方位角ではない.
            ! 【式】rho = C F^(-1/2) (伸長座標での動径),  r = D rho e  (e: (u, v) 方向の単位ベクトル)
            !   d(rho e)/du x d(rho e)/dv = rho^2 sin u e - rho rho_u sin u e_theta - rho rho_v e_phi
            !   実空間では (D a) x (D b) = det(D) D^(-1) (a x b) を使う. dinv = sqrt(f) = D^(-1).
            !   rho_u = -rho F_u / (2F),  rho_v = -rho F_v / (2F)  (F の微分は刻み hd の中心差分)
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
            ! 【目的】点 x [fm] が epsilon 形状の内側かを判定する.
            ! 【式】xi = sqrt(f) x を極座標 (rho_t, theta_t, phi_t) にし, rho_t^2 F(theta_t, phi_t) <= C^2 なら内側.
            ! acos の引数は丸め誤差で |.| > 1 にならないよう [-1, 1] に制限している.
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
            ! 【目的】表面上の |x|, |y|, |z| の最大値 [fm] を 65 x 129 点のサンプリングで求める.
            ! grid_mod が計算箱の大きさを決めるのに使う.
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
            ! 【目的】形状の体積 [fm^3].  式: V = 1/3 oint r . N du dv  (発散定理, div r = 3)
            implicit none
            class(surface_shape), intent(in) :: this
            real(dp) :: vol
            vol = surface_integral(this, 1)
        end function shape_volume

        function shape_relative_area(this) result(B_s)
            ! surface area / (4 pi R^2), 4/3 pi R^3 = volume                     (Eq. 70)
            ! 【目的】相対表面積 B_s = S / (4 pi R^2),  R = (3V / 4 pi)^(1/3)   (FRDM2012 式 70)
            ! S は表面積, V は体積. 平均対ギャップ (式 63-65) と Strutinsky の幅 gamma_S (式 111) に使う.
            implicit none
            class(surface_shape), intent(in) :: this
            real(dp) :: B_s, vol
            vol = surface_integral(this, 1)
            B_s = surface_integral(this, 2) / (4.0_dp * pi * (3.0_dp * vol / (4.0_dp * pi))**(2.0_dp / 3.0_dp))
        end function shape_relative_area

        function surface_integral(this, what) result(s)
            ! what = 1: volume, 2: area.  Gauss-Legendre over the quarter x, y >= 0 times 4
            ! 【目的】what = 1 で体積 V = 1/3 oint r . N du dv,  what = 2 で表面積 S = oint |N| du dv.
            ! 1/4 領域 (u in [0, pi], v in [0, pi/2]) を 128 x 64 点の Gauss-Legendre で積分して 4 倍する.
            ! x, y の各面についての鏡映対称性だけを仮定している (eps3 /= 0 でも使える).
            implicit none
            class(surface_shape), intent(in) :: this
            integer, intent(in) :: what
            real(dp) :: s, r(3), N(3), u, v
            integer, parameter :: nq = 64
            real(dp) :: tq(nq), wq(nq), tt(2 * nq), wt(2 * nq)
            integer :: i, j
            call gauss_legendre(nq, tq, wq)
            call gauss_legendre(2 * nq, tt, wt)
            s = 0.0_dp
            do j = 1, nq
                v = 0.25_dp * pi * (tq(j) + 1.0_dp)
                do i = 1, 2 * nq
                    u = 0.5_dp * pi * (tt(i) + 1.0_dp)
                    call this%eval(u, v, r, N)
                    if (what == 1) then
                        s = s + wt(i) * wq(j) * dot_product(r, N) / 3.0_dp
                    else
                        s = s + wt(i) * wq(j) * norm2(N)
                    end if
                end do
            end do
            s = 4.0_dp * s * (0.5_dp * pi) * (0.25_dp * pi)
        end function surface_integral

        subroutine gauss_legendre(n, x, w)
            ! Gauss-Legendre nodes and weights on [-1, 1] (Newton iteration on P_n)
            ! 【目的】[-1, 1] 上の n 点 Gauss-Legendre 則の節点 x と重み w.
            ! 【式】P_n(z) = 0 をニュートン法で解く (初期値 cos(pi (i - 1/4) / (n + 1/2))).
            !   漸化式 k P_k = (2k - 1) z P_{k-1} - (k - 1) P_{k-2}
            !   微分   P_n' = n (z P_n - P_{n-1}) / (z^2 - 1),  重み w = 2 / ((1 - z^2) P_n'^2)
            ! 節点は 0 について対称なので半分だけ計算し, 残りは鏡映で埋める.
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
