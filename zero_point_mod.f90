module zero_point_mod
    !---------------------------
    ! Zero-point energy of the ground state (FRDM2012 Sec. 2.11, Eqs. 112-118)
    !   E_0,lambda = 1/2 hbar omega_lambda,  omega_lambda = K (C_lambda / B_lambda^irr)^(1/2),  K = frK
    ! for lambda = eps2 and gamma.
    !   B_eps2^irr  = 2/15 (1 + 2/9 eps2^2) / (1 - 2/3 eps2)^2 (1 - 1/3 eps2^2 - 2/27 eps2^3)^(-4/3) M0 R0^2  (116)
    !   B_gamma^irr = 2/15 [(1 - 2/3 eps2)/(1 + 1/3 eps2)]^(2/3) [ln((1 + 1/3 eps2)/(1 - 2/3 eps2))]^2 M0 R0^2 (117)
    ! M0 = nucleus%mass (A m_nucleon), R0 = nucleus%R_a = r0 A^(1/3).
    ! C_lambda = d^2E/dlambda^2 at the ground state from a least-squares parabola through
    ! the seven equally spaced points lambda_gs + k h, k = -3..3 (Eq. 118 and the text below it):
    !   C = [5 (E_3 + E_-3) + 0 (E_2 + E_-2) - 3 (E_1 + E_-1) - 4 E_0] / (42 h^2)
    ! The fit averages over the shell structure; it is not the local second derivative.
    ! With B in MeV fm^2 (masses in MeV) hbar omega = hbar_c sqrt(C / B), C in MeV.
    !---------------------------
    ! 【日本語での説明】
    ! 基底状態の零点エネルギー (FRDM2012 2.11 節). eps2 方向と gamma 方向の振動について
    !   E_0 = 1/2 hbar omega,  omega = K (C / B^irr)^(1/2)   (K = frK = 0.2475)
    ! を計算する. B^irr は非圧縮・渦なし流れの慣性 (式 116, 117), C はポテンシャルのかたさ.
    ! 【C の求め方】基底状態 lambda_gs から等間隔 h で取った 7 点 (k = -3..3) の E に
    ! 最小二乗で放物線をフィットし, C = 2 x (2 次の係数). 等間隔なので重み付き和で書ける.
    !   C = [5 (E_3 + E_-3) - 3 (E_1 + E_-1) - 4 E_0] / (42 h^2)   (E_2, E_-2 の重みは 0)
    ! 殻構造による細かい凹凸をならした曲率なので, 7 点中心差分 (局所的な 2 階微分) とは別物.
    ! 【引数】E_gs_lambda(1:7) は lambda_gs - 3h, ..., lambda_gs + 3h での E_pot [MeV].
    !   gamma の h はラジアンで渡す.
    use iso_fortran_env, only: real64
    use constant_mod, only: hbar_c
    use nucleus_mod, only: nucleus_property
    use micro_constant_mod, only: frK
    implicit none
    private
    integer, parameter :: dp = real64

    ! 零点エネルギーの計算の親 (抽象型). eps2 方向と gamma 方向で慣性の式だけが違う.
    type, abstract, public :: zero_point
        contains
            procedure(zero_point_E), deferred :: cal_energy
    end type zero_point

    abstract interface
        subroutine zero_point_E(this, nucleus, eps2, h, E_gs_lambda, E_0_lambda)
            import :: zero_point, dp, nucleus_property
            class(zero_point), intent(in) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp), intent(in) :: eps2
            real(dp), intent(in) :: h
            real(dp), intent(in) :: E_gs_lambda(7)
            real(dp), intent(out) :: E_0_lambda
        end subroutine zero_point_E
    end interface

    type, extends(zero_point), public :: zero_eps2
        contains
            procedure :: cal_energy => calculate_E_0_eps2
    end type zero_eps2

    type, extends(zero_point), public :: zero_gamma
        contains
            procedure :: cal_energy => calculate_E_0_gamma
    end type zero_gamma

    public :: B_eps2_irr, B_gamma_irr, stiffness_fit

    contains

        subroutine calculate_E_0_eps2(this, nucleus, eps2, h, E_gs_lambda, E_0_lambda)
            ! E_0,eps2 = 1/2 hbar_c K sqrt(C_eps2 / B_eps2^irr)
            ! 【目的】eps2 方向の零点エネルギー E_0,eps2 [MeV]. eps2 は基底状態の値.
            class(zero_eps2), intent(in) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp), intent(in) :: eps2
            real(dp), intent(in) :: h
            real(dp), intent(in) :: E_gs_lambda(7)
            real(dp), intent(out) :: E_0_lambda
            E_0_lambda = zero_point_energy(stiffness_fit(E_gs_lambda, h), B_eps2_irr(nucleus, eps2))
        end subroutine calculate_E_0_eps2

        subroutine calculate_E_0_gamma(this, nucleus, eps2, h, E_gs_lambda, E_0_lambda)
            ! E_0,gamma = 1/2 hbar_c K sqrt(C_gamma / B_gamma^irr)   (h in radian)
            ! 【目的】gamma 方向の零点エネルギー E_0,gamma [MeV].
            ! eps2 は基底状態の値 (B_gamma が依存する).
            class(zero_gamma), intent(in) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp), intent(in) :: eps2
            real(dp), intent(in) :: h
            real(dp), intent(in) :: E_gs_lambda(7)
            real(dp), intent(out) :: E_0_lambda
            E_0_lambda = zero_point_energy(stiffness_fit(E_gs_lambda, h), &
                                           B_gamma_irr(nucleus, eps2))
        end subroutine calculate_E_0_gamma

        function zero_point_energy(C, B) result(E0)
            ! E_0 = 1/2 hbar omega,  hbar omega = K hbar_c sqrt(C / B^irr)   (Eqs. 112-115)
            ! 【目的】かたさ C [MeV] と渦なし慣性 B^irr [MeV fm^2] から E_0 [MeV].
            ! C <= 0 (極小になっていない) なら振動数が定義できないので停止する.
            real(dp), intent(in) :: C, B
            real(dp) :: E0
            if (C <= 0.0_dp) error stop "zero_point_mod: stiffness C <= 0 (not a minimum)"
            E0 = 0.5_dp * frK * hbar_c * sqrt(C / B)
        end function zero_point_energy

        pure function B_eps2_irr(nucleus, eps2) result(B)
            ! Eq. (116), MeV fm^2
            ! 【目的】eps2 方向の渦なし流れの慣性 B_eps2^irr [MeV fm^2] (式 116).
            type(nucleus_property), intent(in) :: nucleus
            real(dp), intent(in) :: eps2
            real(dp) :: B
            B = 2.0_dp / 15.0_dp * (1.0_dp + 2.0_dp / 9.0_dp * eps2**2) / (1.0_dp - 2.0_dp / 3.0_dp * eps2)**2 &
                * (1.0_dp - 1.0_dp / 3.0_dp * eps2**2 - 2.0_dp / 27.0_dp * eps2**3)**(-4.0_dp / 3.0_dp) &
                * nucleus%mass * nucleus%R_a**2
        end function B_eps2_irr

        pure function B_gamma_irr(nucleus, eps2) result(B)
            ! Eq. (117), MeV fm^2 (per radian^2); vanishes like eps2^2 for the sphere
            ! 【目的】gamma 方向の渦なし流れの慣性 B_gamma^irr [MeV fm^2] (式 117). 球で 0 になる.
            type(nucleus_property), intent(in) :: nucleus
            real(dp), intent(in) :: eps2
            real(dp) :: B
            B = 2.0_dp / 15.0_dp &
                * ((1.0_dp - 2.0_dp / 3.0_dp * eps2) / (1.0_dp + 1.0_dp / 3.0_dp * eps2))**(2.0_dp / 3.0_dp) &
                * log((1.0_dp + 1.0_dp / 3.0_dp * eps2) / (1.0_dp - 2.0_dp / 3.0_dp * eps2))**2 &
                * nucleus%mass * nucleus%R_a**2
        end function B_gamma_irr

        pure function stiffness_fit(E, h) result(C)
            ! C = 2 a of the least-squares parabola E = a x^2 + b x + c through x_k = k h, k = -3..3:
            !   with S_n = sum k^n (S_0 = 7, S_2 = 28, S_4 = 196, odd S_n = 0),
            !   a h^2 = (S_0 sum k^2 E_k - S_2 sum E_k) / (S_0 S_4 - S_2^2) = (sum k^2 E_k - 4 sum E_k) / 84
            ! 【目的】等間隔 7 点 E(1:7) (x = -3h..3h) の最小二乗放物線から C = d^2E/dx^2 = 2a.
            !   重み (k = 0, +-1, +-2, +-3) = (-4, -3, 0, 5) / (42 h^2). E = x^2 なら C = 2 が厳密に出る.
            real(dp), intent(in) :: E(7), h
            real(dp) :: C
            C = (5.0_dp * (E(1) + E(7)) - 3.0_dp * (E(3) + E(5)) - 4.0_dp * E(4)) / (42.0_dp * h**2)
        end function stiffness_fit

end module zero_point_mod
