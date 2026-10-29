program main
    !---------------------------
    ! FRLDM + Shell + pairing (FRDM2012 macroscopic-microscopic model, finite-range
    ! liquid-drop macroscopic part) for one nucleus.
    !
    !   M_th = E_mac(Z, N, shape) + E_shell + E_pair + E_zp
    !     E_mac           : FRLDM, Eq. (62), as a mass excess (frldm_mod)
    !     E_shell, E_pair : Strutinsky shell correction and Lipkin-Nogami (or BCS) pairing
    !                       correction (shell_bcs_mod)
    !     E_zp            : zero-point energy of the eps2 and gamma vibrations, Sec. 2.11
    !                       (zero_point_mod; only with zpe = .true.)
    !   binding energy B = Z M_H + N M_n - M_th
    !
    ! The shape is the FRDM epsilon parametrisation (eps2, gamma, eps3, eps4, eps6), Eqs. (17)-(25);
    ! eps1 follows from eps3 by the centre-of-mass condition (shape_mod).
    ! With search = .true. the potential energy E_mac + E_shell + E_pair is minimised and the
    ! final result is given at the minimum. search_method = 'pattern' (default) minimises in
    ! (eps2, eps3, eps4, eps6) on a lattice of step d_eps by a pattern search from several
    ! starting points, as FRDM2012 (Sec. 3, step 3); 'grid' scans the (eps2, eps4) grid.
    ! With zpe = .true. (default) the zero-point energy is added at the final shape: E_pot is
    ! computed at eps2 +- 0.05, 0.10, 0.15 (and gamma +- 5, 10, 15 deg when eps2 >= 0.17) and the
    ! stiffness C is fitted by a least-squares parabola. For eps2 < 0.17, E_0,gamma = E_0,eps2.
    !
    ! Input: namelist &frldm_input in the file given as the first argument (default
    ! input.nml); without the file the defaults below (208Pb, spherical) are used.
    !   Z, N                 : proton and neutron numbers
    !   eps2, gamma_deg, eps3, eps4, eps6 : deformation (gamma in degree)
    !   pairing              : 'LN' (default) or 'BCS'
    !   n_gh                 : Gauss-Hermite points per axis of the single-particle problem
    !   use_grid, h, n_times : B_1, B_3 from the grid + Hockney method instead of the surface
    !                          integrals (cross-check; grid spacing h fm, box n_times x nucleus)
    !   search               : ground-state search
    !   search_method        : 'pattern' (default) or 'grid'
    !   eps2_min, eps2_max, eps3_max, eps4_min, eps4_max, eps6_min, eps6_max : search range
    !                          (eps3 from 0 to eps3_max, since E(-eps3) = E(eps3); min = max fixes
    !                          the variable)
    !   d_eps, d_eps_start   : final and first step of the pattern search
    !   d_eps2, d_eps4       : steps of the 'grid' search
    !   n_gh_search          : Gauss-Hermite points during the search (0: same as n_gh)
    !   zpe                  : add the zero-point energy (6 or 12 more shapes)
    !---------------------------
    ! 【日本語での説明】
    ! 1 つの原子核について FRLDM + 殻補正 + 対補正の質量を計算する.
    !   M_th = E_mac(Z, N, 形状) + E_shell + E_pair + E_zp   (質量超過 [MeV])
    !   B    = Z M_H + N M_n - M_th                    (結合エネルギー [MeV])
    !   E_mac : frldm_mod (FRDM2012 式 62),  E_shell, E_pair : shell_bcs_mod
    ! 入力は namelist &frldm_input (第 1 引数のファイル, 既定 input.nml).
    ! search = .true. なら E_mac + E_shell + E_pair が最小の形状を探し, その形状で詳しい結果を出力する.
    !   search_method = 'pattern' (既定): (eps2, eps3, eps4, eps6) の刻み d_eps の格子上で,
    !     いくつかの初期点からパターン探索 (FRDM2012 3 章の手順 3 と同じ離散最小化) をする.
    !   search_method = 'grid': これまでどおり (eps2, eps4) の格子を全部計算する.
    ! zpe = .true. (既定) なら, 最後の形状で零点エネルギー E_zp (2.11 節, zero_point_mod) を足す.
    !   eps2 +- 0.05, 0.10, 0.15 (eps2 >= 0.17 なら gamma +- 5, 10, 15 度も) で E_pot を計算し,
    !   7 点の最小二乗放物線からかたさ C を求める. eps2 < 0.17 では E_0,gamma = E_0,eps2.
    !$ use omp_lib
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use nucleus_mod, only: nucleus_property
    use grid_mod, only: grid_type
    use shape_mod, only: epsilon_shape, epsilon_shape_from
    use frldm_mod, only: frldm_variables, M_H, M_n
    use micro_constant_mod, only: microscopic_variables
    use shell_bcs_mod, only: shell_bcs_result, shell_bcs_calculation, nucleon_correction
    use zero_point_mod, only: zero_eps2, zero_gamma, stiffness_fit

    implicit none
    integer, parameter :: dp = real64

    ! 入力変数 (namelist &frldm_input) と既定値 (208Pb, 球形):
    !   Z, N        : 陽子数, 中性子数
    !   eps2, gamma_deg, eps3, eps4, eps6 : 変形 (gamma は度). eps1 は eps3 から重心条件で決まる
    !   pairing     : 'LN' (Lipkin-Nogami, 既定) または 'BCS'
    !   n_gh        : 一粒子問題の Gauss-Hermite 点数 (1 軸あたり)
    !   use_grid, h, n_times : B_1, B_3 を格子 + Hockney 法で求める (検算用;
    !                          格子間隔 h [fm], 箱は核の n_times 倍)
    !   search      : 基底状態探索をするか
    !   search_method : 'pattern' (既定) または 'grid'
    !   eps2_min, eps2_max, eps3_max, eps4_min, eps4_max, eps6_min, eps6_max : 探索範囲
    !                 (E(-eps3) = E(eps3) なので eps3 は 0 から eps3_max. min = max ならその変数は固定)
    !   d_eps, d_eps_start : パターン探索の最後と最初の刻み
    !   d_eps2, d_eps4     : 'grid' の刻み
    !   n_gh_search : 探索中の Gauss-Hermite 点数 (0 なら n_gh と同じ)
    !   zpe         : 零点エネルギーを足すか (形状を 6 点または 12 点余分に計算する)
    integer :: Z = 82, N = 126, n_gh = 44, n_times = 3, n_gh_search = 0
    real(dp) :: eps2 = 0.0_dp, gamma_deg = 0.0_dp, eps4 = 0.0_dp, eps6 = 0.0_dp, h = 0.5_dp
    real(dp) :: eps1 = 0.0_dp, eps3 = 0.0_dp
    character(8) :: pairing = 'LN'
    logical :: use_grid = .false., search = .false., zpe = .true.
    logical :: make_table = .false.
    ! data_table の作成 
    integer :: t1, t2, t_rate, t_max, diff
    character(8) :: search_method = 'pattern'
    real(dp) :: eps2_min = -0.30_dp, eps2_max = 0.45_dp, d_eps2 = 0.05_dp
    real(dp) :: eps4_min = -0.12_dp, eps4_max = 0.12_dp, d_eps4 = 0.04_dp
    real(dp) :: eps3_max = 0.20_dp, eps6_min = -0.08_dp, eps6_max = 0.08_dp
    real(dp) :: d_eps = 0.01_dp, d_eps_start = 0.04_dp
    namelist /frldm_input/ Z, N, eps2, gamma_deg, eps3, eps4, eps6, pairing, n_gh, use_grid, h, n_times, &
                           search, search_method, eps2_min, eps2_max, d_eps2, eps4_min, eps4_max, d_eps4, &
                           eps3_max, eps6_min, eps6_max, d_eps, d_eps_start, n_gh_search, zpe, make_table

    character(256) :: file_name
    logical :: exists
    integer :: unit, ios, n_gh_final
    ! cache of pattern_search: lattice points c_x(:, i) with E_pot c_e(i), n_eval shapes computed
    integer, allocatable :: c_x(:,:)
    real(dp), allocatable :: c_e(:)
    integer :: n_eval = 0
    real(dp) :: e_pot
    type(frldm_variables) :: frldm
    type(microscopic_variables) :: mv
    type(shell_bcs_result) :: res

    character(256) :: Z_num, N_num
    real(dp) :: E_bind, B_A, E_zp, M_th, e_mac_out, e_sp_out

    call system_clock(t1)

    file_name = 'input.nml'
    if (command_argument_count() >= 1) call get_command_argument(1, file_name)
    inquire(file=trim(file_name), exist=exists)
    if (exists) then
        open(newunit=unit, file=trim(file_name), status='old', action='read')
        read(unit, nml=frldm_input, iostat=ios)
        close(unit)
        if (ios /= 0) stop "main: error reading namelist &frldm_input"
        print '(a, a)', " input: ", trim(file_name)
    else
        print '(a, a, a)', " input: ", trim(file_name), " not found, defaults used"
    end if

    if (search) then
        n_gh_final = n_gh
        if (n_gh_search > 0) n_gh = n_gh_search
        if (search_method == 'grid') then
            call grid_search()
        else if (search_method == 'pattern') then
            call pattern_search()
        else
            stop "main: search_method must be 'pattern' or 'grid'"
        end if
        n_gh = n_gh_final
    end if

    call calculate(eps2, gamma_deg, eps3, eps4, eps6, .true., e_pot)
    ! E_mac, E_s+p of the final shape: add_zero_point overwrites frldm and res
    e_mac_out = frldm%E_frldm
    e_sp_out = res%e_micro
    if (zpe) then
        call add_zero_point(e_pot, E_zp, M_th, E_bind, B_A)
    else
        E_zp = 0.0_dp
        M_th = e_pot
        E_bind = M_H * Z + M_n * N - M_th
        B_A = E_bind / (Z + N)
    end if

    ! Output the data to a file
    write(*,*) "Writing results to text file..."
    write(Z_num, '(I3.3)') Z
    write(N_num, '(I3.3)') N
    ! the files go to data/ under the current directory; if it does not exist, nothing is written
    open(newunit=unit, file="data/results"//trim(Z_num)//"_"//trim(N_num)//".txt", status="replace", &
         action="write", iostat=ios)
    if (ios /= 0) then
        write(*,*) "warning: cannot open data/results"//trim(Z_num)//"_"//trim(N_num)//".txt (no data directory?)"
    else
        write(unit,*) "Z = ", Z
        write(unit,*) " N = ", N
        write(unit,*) " eps2 = ", eps2
        write(unit,*) " eps3 = ", eps3
        write(unit,*) " eps4 = ", eps4
        write(unit,*) " eps6 = ", eps6
        write(unit,*) " E_mac = ", e_mac_out
        write(unit,*) " E_shell + E_pair = ", e_sp_out
        write(unit,*) " E_zp = ", E_zp
        write(unit,*) " M_th = ", M_th
        write(unit,*) " B = ", E_bind
        write(unit,*) " B / A = ", B_A
        close(unit)
        write(*,*) "Results written to data/results"//trim(Z_num)//"_"//trim(N_num)//".txt"
    end if
    if (make_table) then
        open(newunit=unit, file="data/data_table.txt", status="unknown", position="append", action="write", &
             iostat=ios)
        if (ios /= 0) then
            write(*,*) "warning: cannot open data/data_table.txt (no data directory?)"
        else
            write(unit,*) Z, N, eps2, eps3, eps4, eps6, e_mac_out, e_sp_out, E_zp, M_th, E_bind, B_A
            close(unit)
            write(*,*) "Results appended to data/data_table.txt"
        end if
    end if

    write(*,*) "Finished to calculate the mass of the nucleus."

    call system_clock(t2, t_rate, t_max)
    if (t1 > t2) then
        diff = t_max - t1 + t2 + 1
        print '(a, f8.2, a)', " total time: ", real(diff, dp) / real(t_rate, dp), " s"
    else
        diff = t2 - t1
        print '(a, f8.2, a)', " total time: ", real(diff, dp) / real(t_rate, dp), " s"
    end if

contains

    subroutine grid_search()
        ! E_pot on the full (eps2, eps4) grid, eps3 and eps6 fixed at the input values
        ! 【日本語】(eps2, eps4) の格子を全部計算し, 最小の点を eps2, eps4 に入れる.
        integer :: i2, i4, n2, n4
        real(dp) :: e_best, best(2)

        n2 = nint((eps2_max - eps2_min) / d_eps2)
        n4 = nint((eps4_max - eps4_min) / d_eps4)
        print '(a)', " ground-state search: E_pot = E_mac + E_shell + E_pair [MeV]"
        print '(a)', "     eps2     eps4        E_mac     E_s+p      E_pot"
        e_best = huge(1.0_dp)
        do i2 = 0, n2
            do i4 = 0, n4
                call calculate(eps2_min + i2 * d_eps2, gamma_deg, eps3, eps4_min + i4 * d_eps4, eps6, .false., e_pot)
                print '(2f9.3, 3f11.4)', eps2_min + i2 * d_eps2, eps4_min + i4 * d_eps4, frldm%E_frldm, &
                      res%e_micro, e_pot
                if (e_pot < e_best) then
                    e_best = e_pot
                    best = [eps2_min + i2 * d_eps2, eps4_min + i4 * d_eps4]
                end if
            end do
        end do
        print '(a, 2f8.3)', " minimum on the grid at eps2, eps4 = ", best
        eps2 = best(1)
        eps4 = best(2)
    end subroutine grid_search

    subroutine pattern_search()
        ! Discrete minimisation of E_pot in (eps2, eps3, eps4, eps6) on the lattice eps = d_eps * i,
        ! as FRDM2012 (Sec. 3, step 3). From each starting point: compute the 8 neighbours
        ! x +- s e_k, move to the lowest one if it is below x, otherwise halve s (from
        ! d_eps_start / d_eps down to 1). Computed points are cached, so overlapping paths are
        ! free. Starting points: the input shape, (eps2, eps4) = (+-0.25, +-0.04) and spherical,
        ! each with eps3 = 0 and 0.10 (eps6 = 0). eps3 >= 0 only, since E(-eps3) = E(eps3).
        ! gamma is fixed at gamma_deg. The lowest local minimum is put in eps2, eps3, eps4, eps6.
        ! 【日本語】
        ! 【目的】(eps2, eps3, eps4, eps6) の刻み d_eps の格子上で E_pot の極小を探し,
        ! いちばん低いものを基底状態の形として eps2, eps3, eps4, eps6 に入れる.
        ! 【手順】1. 初期点: 入力の形, (eps2, eps4) = (+-0.25, +-0.04) と球形. それぞれ eps3 = 0 と 0.10.
        !         2. 今の点 x から各軸 +-s の 8 近傍を計算し, x より低い最小の近傍へ移る.
        !            低い近傍がなければ s を半分にし (d_eps_start -> ... -> d_eps), s = d_eps で
        !            下がらなければそこが極小.
        !         3. 計算した点は覚えておき (キャッシュ), 同じ点は計算し直さない.
        ! eps3 は E(-eps3) = E(eps3) なので 0 以上だけを調べる. gamma は gamma_deg に固定.
        integer, parameter :: max_start = 11
        integer :: lo(4), hi(4), x(4), y(4), xb(4), x_glob(4), x_start(4, max_start)
        integer :: s, s0, k, sgn, is, js, n_start, i_s2, i_s4, i_s3
        real(dp) :: e_x, e_y, e_b, e_glob, starts(4, max_start)
        logical :: dup

        lo = [nint(eps2_min / d_eps), 0, nint(eps4_min / d_eps), nint(eps6_min / d_eps)]
        hi = [nint(eps2_max / d_eps), nint(eps3_max / d_eps), nint(eps4_max / d_eps), nint(eps6_max / d_eps)]
        if (any(hi < lo)) stop "main: pattern_search needs min <= max for every eps"
        s0 = max(1, nint(d_eps_start / d_eps))

        ! starting points (rounded to the lattice and clipped to the range; duplicates removed)
        n_start = 1
        starts(:, 1) = [eps2, abs(eps3), eps4, eps6]
        do i_s3 = 0, 1
            do i_s2 = -1, 1, 2
                do i_s4 = -1, 1, 2
                    n_start = n_start + 1
                    starts(:, n_start) = [0.25_dp * i_s2, 0.10_dp * i_s3, 0.04_dp * i_s4, 0.0_dp]
                end do
            end do
            n_start = n_start + 1
            starts(:, n_start) = [0.0_dp, 0.10_dp * i_s3, 0.0_dp, 0.0_dp]
        end do
        js = 0
        do is = 1, n_start
            x = min(max(nint(starts(:, is) / d_eps), lo), hi)
            dup = .false.
            do k = 1, js
                if (all(x_start(:, k) == x)) dup = .true.
            end do
            if (dup) cycle
            js = js + 1
            x_start(:, js) = x
        end do
        n_start = js

        allocate(c_x(4, 0), c_e(0))
        n_eval = 0
        print '(a)', " ground-state search (pattern search): E_pot = E_mac + E_shell + E_pair [MeV]"
        print '(a, f6.3, a, f6.3, a, i3)', " step ", s0 * d_eps, " -> ", d_eps, ",  starting points: ", n_start
        print '(a)', "     eps2     eps3     eps4     eps6        E_mac     E_s+p      E_pot"
        e_glob = huge(1.0_dp)
        x_glob = x_start(:, 1)
        do is = 1, n_start
            x = x_start(:, is)
            print '(a, i3, a, 4f9.3)', " -- start ", is, " at ", x * d_eps
            call energy_at(x, e_x)
            s = s0
            do
                xb = x
                e_b = e_x
                do k = 1, 4
                    if (lo(k) == hi(k)) cycle
                    do sgn = -1, 1, 2
                        y = x
                        y(k) = x(k) + sgn * s
                        if (y(k) < lo(k) .or. y(k) > hi(k)) cycle
                        call energy_at(y, e_y)
                        if (e_y < e_b) then
                            e_b = e_y
                            xb = y
                        end if
                    end do
                end do
                if (any(xb /= x)) then
                    x = xb
                    e_x = e_b
                else if (s > 1) then
                    s = max(1, s / 2)
                else
                    exit
                end if
            end do
            print '(a, 4f9.3, a, f11.4)', " local minimum at", x * d_eps, "   E_pot ", e_x
            if (e_x < e_glob) then
                e_glob = e_x
                x_glob = x
            end if
        end do
        print '(a, i6)', " shapes computed: ", n_eval
        print '(a, 4f8.3, a, f11.4)', " lowest minimum at eps2, eps3, eps4, eps6 = ", x_glob * d_eps, &
              "   E_pot ", e_glob
        eps2 = x_glob(1) * d_eps
        eps3 = x_glob(2) * d_eps
        eps4 = x_glob(3) * d_eps
        eps6 = x_glob(4) * d_eps

    end subroutine pattern_search

    ! 【目的】格子点 xi の E_pot を返す. 計算済みならキャッシュの値を使い, なければ計算して覚える.
    subroutine energy_at(xi, e)
        ! E_pot at the lattice point xi (eps = d_eps * xi), from the cache of pattern_search if computed
        integer, intent(in) :: xi(4)
        real(dp), intent(out) :: e
        integer :: i
        do i = 1, size(c_e)
            if (all(c_x(:, i) == xi)) then
                e = c_e(i)
                return
            end if
        end do
        call calculate(xi(1) * d_eps, gamma_deg, xi(2) * d_eps, xi(3) * d_eps, xi(4) * d_eps, .false., e)
        n_eval = n_eval + 1
        print '(4f9.3, 3f11.4)', xi * d_eps, frldm%E_frldm, res%e_micro, e
        c_x = reshape([c_x, xi], [4, size(c_e) + 1])
        c_e = [c_e, e]
    end subroutine energy_at

    subroutine add_zero_point(e_gs, E_zp, M_th, E_bind, B_A)
        ! E_zp = E_0,eps2 + E_0,gamma at the shape (eps2, gamma_deg, eps3, eps4, eps6), FRDM2012 Sec. 2.11.
        ! E_pot at lambda_gs + k h (k = -3..3, h = 0.05 for eps2, 5 deg for gamma), other
        ! deformations fixed; C from the least-squares parabola (zero_point_mod).
        ! eps2 < 0.17: E_0,gamma = E_0,eps2 (the gamma curvature is not computed).
        ! Triaxial ground states (interpolation of the inertia, Sec. 2.11) are not treated.
        ! The shape is not a minimum when the fitted stiffness C <= 0 (light, soft nuclei):
        ! then E_0,eps2 = 0, or E_0,gamma = E_0,eps2 as for eps2 < 0.17, with a warning.
        ! E_mac and E_s+p of the ground state are taken from the last call of calculate
        ! (the one that gave e_gs) before the scan overwrites them.
        ! 【日本語】
        ! 【目的】最後の形状で零点エネルギー E_zp を求め, M_th + E_zp を表示する.
        ! 【手順】1. eps2 を +-0.05, +-0.10, +-0.15 ずらした 6 点で E_pot を計算 (他の変形は固定).
        !         2. 基底状態の e_gs と合わせた 7 点から zero_eps2 で E_0,eps2.
        !         3. eps2 >= 0.17 なら gamma を +-5, +-10, +-15 度ずらして同様に E_0,gamma.
        !            eps2 < 0.17 なら E_0,gamma = E_0,eps2 (FRDM2012 2.11 節).
        ! かたさ C <= 0 (その方向で極小になっていない. 軽くやわらかい核で起こる) のときは止めずに
        ! 警告を出し, eps2 方向は E_0,eps2 = 0, gamma 方向は E_0,gamma = E_0,eps2 とする.
        ! 基底状態の E_mac, E_s+p は, ずらした形の計算で上書きされる前に控えておく.
        real(dp), intent(in) :: e_gs
        real(dp), parameter :: h_eps = 0.05_dp, h_gam_deg = 5.0_dp, eps2_gamma_min = 0.17_dp
        type(nucleus_property) :: nuc
        type(zero_eps2) :: zp_eps2
        type(zero_gamma) :: zp_gamma
        real(dp) :: E7(7), e0_eps2, e0_gamma, e_mac_gs, e_sp_gs
        real(dp), intent(out) :: E_bind, B_A, M_th, E_zp
        integer :: k

        e_mac_gs = frldm%E_frldm
        e_sp_gs = res%e_micro

        nuc%Z = Z
        nuc%N = N
        nuc%semi1 = nuc%R0 * (Z + N)**(1.0_dp / 3.0_dp)
        call nuc%calculate_fundamental_properties()

        E7(4) = e_gs
        do k = -3, 3
            if (k /= 0) call calculate(eps2 + k * h_eps, gamma_deg, eps3, eps4, eps6, .false., E7(k + 4))
        end do
        print '(a)', " --- zero-point energy (Sec. 2.11) [MeV]"
        print '(a, 7f10.4)', "   E_pot(eps2 - 0.15 .. eps2 + 0.15)     : ", E7
        if (stiffness_fit(E7, h_eps) > 0.0_dp) then
            call zp_eps2%cal_energy(nuc, eps2, h_eps, E7, e0_eps2)
        else
            print '(a)', "   warning: C_eps2 <= 0, the shape is not a minimum in eps2; E_0,eps2 = 0"
            e0_eps2 = 0.0_dp
        end if

        if (eps2 < eps2_gamma_min) then
            e0_gamma = e0_eps2
        else
            do k = -3, 3
                if (k /= 0) call calculate(eps2, gamma_deg + k * h_gam_deg, eps3, eps4, eps6, .false., E7(k + 4))
            end do
            print '(a, 7f10.4)', "   E_pot(gamma - 15 .. gamma + 15 deg)   : ", E7
            if (stiffness_fit(E7, h_gam_deg * pi / 180.0_dp) > 0.0_dp) then
                call zp_gamma%cal_energy(nuc, eps2, h_gam_deg * pi / 180.0_dp, E7, e0_gamma)
            else
                print '(a)', "   warning: C_gamma <= 0, the shape is not a minimum in gamma; E_0,gamma = E_0,eps2"
                e0_gamma = e0_eps2
            end if
        end if

        E_zp = e0_eps2 + e0_gamma
        M_th = e_gs + E_zp
        E_bind = M_H * Z + M_n * N - M_th
        B_A = E_bind / (Z + N)
        print '(a, f12.4)', "   E_0,eps2              : ", e0_eps2
        if (eps2 < eps2_gamma_min) then
            print '(a, f12.4, a)', "   E_0,gamma             : ", e0_gamma, "   (= E_0,eps2, eps2 < 0.17)"
        else
            print '(a, f12.4)', "   E_0,gamma             : ", e0_gamma
        end if
        print '(a, f12.4)', "   E_zp                  : ", E_zp
        print '(a)', " --- result with the zero-point energy [MeV]"
        print '(a, f12.4)', "   E_mac (ground state)  : ", e_mac_gs
        print '(a, f12.4)', "   E_s+p (ground state)  : ", e_sp_gs
        print '(a, f12.4)', "   E_zp                  : ", E_zp
        print '(a, f12.4)', "   mass excess M_th      : ", M_th
        print '(a, f12.4)', "   binding energy B      : ", E_bind
        print '(a, f12.4)', "   B / A                 : ", B_A
    end subroutine add_zero_point

    subroutine calculate(e2_, gam_, e3_, e4_, e6_, verbose, e_pot)
        ! FRLDM + shell + pairing at one deformation; e_pot = M_th (mass excess, MeV)
        ! 【日本語】
        ! 【目的】1 つの変形 (e2_, gam_ [度], e3_, e4_, e6_) で M_th を計算し, e_pot に返す.
        ! verbose = .true. なら各項を表示する.
        ! 【手順】1. nucleus に Z, N と変形を入れ, 半径 R_a の epsilon 形状を作る.
        !         2. frldm_mod で E_mac (use_grid なら格子を用意して B_1, B_3 を格子で計算).
        !         3. shell_bcs_mod で E_shell + E_pair.  4. M_th = E_mac + E_s+p, B = Z M_H + N M_n - M_th.
        real(dp), intent(in) :: e2_, gam_, e3_, e4_, e6_
        logical, intent(in) :: verbose
        real(dp), intent(out) :: e_pot
        type(nucleus_property) :: nucleus
        type(grid_type) :: g_mod
        type(epsilon_shape) :: shape
        real(dp) :: M_th, B
        integer :: t0, t1, rate

        call system_clock(t0, rate)
        nucleus%Z = Z
        nucleus%N = N
        nucleus%semi1 = nucleus%R0 * (Z + N)**(1.0_dp / 3.0_dp)
        nucleus%epsilon2 = e2_
        nucleus%gamma_deg = gam_
        nucleus%epsilon3 = e3_
        nucleus%epsilon4 = e4_
        nucleus%epsilon6 = e6_
        call nucleus%calculate_fundamental_properties()
        shape = epsilon_shape_from(nucleus%R_a, e2_, gam_ * pi / 180.0_dp, e3_, e4_, e6_)
        eps1 = shape%eps1

        ! macroscopic part
        frldm%verbose = .false.
        if (use_grid) then
            g_mod%h_x = h
            g_mod%h_y = h
            g_mod%h_z = h
            g_mod%n_times = n_times
            call g_mod%initialize_grid(nucleus)
            allocate(g_mod%density_index3D(g_mod%n_x_points, g_mod%n_y_points, g_mod%n_z_points))
            call g_mod%inside_outside_nucleus3D(nucleus, g_mod%density_index3D)
            call frldm%calculate_frldm_energy(nucleus, shape, frldm%E_frldm, g_mod)
        else
            call frldm%calculate_frldm_energy(nucleus, shape, frldm%E_frldm)
        end if

        ! microscopic part
        call shell_bcs_calculation(nucleus, mv, res, n_gh, trim(pairing))

        M_th = frldm%E_frldm + res%e_micro
        B = M_H * Z + M_n * N - M_th
        e_pot = M_th
        call system_clock(t1)
        if (.not. verbose) return

        print *
        print '(a, i4, a, i4, a, i4)', " ===== Z = ", Z, "  N = ", N, "  A = ", Z + N
        print '(a, 5f9.4)', " eps2, gamma [deg], eps3, eps4, eps6  : ", e2_, gam_, e3_, e4_, e6_
        print '(a, f12.8)', " eps1 (centre of mass at the origin) : ", eps1
        print '(a, 3f10.6)', " B_s, B_1, B_3                        : ", frldm%B_s, frldm%B_1, frldm%B_3
        print '(a, 2f10.4)', " R_pot [fm], Strutinsky range [MeV]   : ", mv%R_pot, res%gamma_S
        print '(a)', " --- macroscopic energy (FRLDM, Eq. 62) [MeV]"
        print '(a, f12.4)', "   M_H Z + M_n N         : ", frldm%E_mass
        print '(a, f12.4)', "   volume                : ", frldm%E_v
        print '(a, f12.4)', "   surface               : ", frldm%E_s
        print '(a, f12.4)', "   a0                    : ", frldm%E_a0
        print '(a, f12.4)', "   Coulomb               : ", frldm%E_c
        print '(a, f12.4)', "   Coulomb exchange      : ", frldm%E_cex
        print '(a, f12.4)', "   proton form factor    : ", frldm%E_pff
        print '(a, f12.4)', "   charge asymmetry      : ", frldm%E_ca
        print '(a, f12.4)', "   Wigner                : ", frldm%E_W
        print '(a, f12.4)', "   average pairing       : ", frldm%E_pair
        print '(a, f12.4)', "   electron binding      : ", frldm%E_el
        print '(a, f12.4)', "   E_mac                 : ", frldm%E_frldm
        print '(a, a, a)', " --- microscopic corrections (Strutinsky + ", trim(pairing), ") [MeV]"
        call report("neutrons", res%neutron)
        call report("protons ", res%proton)
        print '(a, f12.4)', "   E_shell               : ", res%e_shell
        print '(a, f12.4)', "   E_pair                : ", res%e_pair
        print '(a, f12.4)', "   E_s+p                 : ", res%e_micro
        print '(a)', " --- result [MeV]"
        print '(a, f12.4)', "   mass excess M_th      : ", M_th
        print '(a, f12.4)', "   binding energy B      : ", B
        print '(a, f12.4)', "   B / A                 : ", B / (Z + N)
        print '(a, f8.2, a)', "   time: ", real(t1 - t0, dp) / rate, " s"
    end subroutine calculate

    ! 【目的】中性子または陽子の殻補正・対補正の内訳を 1 行で表示する.
    subroutine report(name, nc)
        character(*), intent(in) :: name
        type(nucleon_correction), intent(in) :: nc
        print '(3x, a, a, f9.3, a, f7.3, a, f7.3, a, f7.3, a, f8.3)', name, ": E_shell ", nc%e_shell, &
              "  Delta_G ", nc%delta_G, "  Delta ", nc%delta, "  lambda2 ", nc%lambda2, "  E_pair ", nc%e_pair
    end subroutine report

end program main
