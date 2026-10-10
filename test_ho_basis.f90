program test_ho_basis
    !---------------------------
    ! Tests for ho_basis_mod
    !   1. orthonormality of the basis functions on the Gauss-Hermite points
    !   2. triaxial T + V_HO must reproduce the analytic oscillator energies
    !   3. spherical Woods-Saxon + spin-orbit must show the (2j+1) degeneracies
    !      1s1/2(2), 1p3/2(4), 1p1/2(2), 1d5/2(6)
    !
    ! build:
    !   gfortran -O2 -fopenmp constant_mod.f90 ho_basis_mod.f90 test_ho_basis.f90 -llapack -o test_ho_basis.exe
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use ho_basis_mod
    implicit none
    integer, parameter :: dp = real64

    type(ho_basis_type) :: basis
    real(dp), allocatable :: T(:,:), Vm(:,:), Axy(:,:), Ayz(:,:), Azx(:,:)
    real(dp), allocatable :: V(:,:,:), dVx(:,:,:), dVy(:,:,:), dVz(:,:,:)
    real(dp), allocatable :: e(:), e_exact(:)
    complex(dp), allocatable :: H(:,:), C(:,:)
    real(dp) :: hw(3), err, x, y, z, r, f, dfdr, kappa
    real(dp), parameter :: V0 = 50.0_dp, R_ws = 1.2_dp * 40.0_dp**(1.0_dp / 3.0_dp), a_ws = 0.65_dp
    integer :: i, j, k, a, ns, n, ng, n_group, group_size(8)
    integer, parameter :: expected(4) = [2, 4, 2, 6]
    logical :: ok_all

    ok_all = .true.

    !---------------------------
    ! 1. orthonormality
    !---------------------------
    hw = frequencies_from_eps_gamma(41.0_dp / 40.0_dp**(1.0_dp / 3.0_dp), 0.25_dp, 20.0_dp * pi / 180.0_dp)
    call basis%initialize(hw, 8, 30)
    err = 0.0_dp
    do j = 0, basis%n_ax_max
        do i = 0, basis%n_ax_max
            f = sum(basis%u(i, :) * basis%u(j, :))
            if (i == j) f = f - 1.0_dp
            err = max(err, abs(f))
        end do
    end do
    call report("orthonormality of u_n(k)", err, 1.0e-12_dp)

    !---------------------------
    ! 2. triaxial harmonic oscillator: V_HO = sum_i 1/2 hw_i xi_i^2
    !---------------------------
    ns = basis%n_states
    ng = basis%n_gh
    allocate(T(ns, ns), Vm(ns, ns), V(ng, ng, ng))
    do k = 1, ng
        do j = 1, ng
            do i = 1, ng
                V(i, j, k) = 0.5_dp * (hw(1) * basis%xi(i)**2 + hw(2) * basis%xi(j)**2 + hw(3) * basis%xi(k)**2)
            end do
        end do
    end do
    call basis%kinetic_matrix(T)
    call basis%potential_matrix(V, Vm)
    allocate(Axy(ns, ns), Ayz(ns, ns), Azx(ns, ns), H(2 * ns, 2 * ns), C(2 * ns, 2 * ns), e(2 * ns))
    Axy = 0.0_dp; Ayz = 0.0_dp; Azx = 0.0_dp
    call basis%build_hamiltonian(T, Vm, Axy, Ayz, Azx, 0.0_dp, H)
    call diagonalize_hamiltonian(H, e, C)

    allocate(e_exact(2 * ns))
    do a = 1, ns
        e_exact(2 * a - 1) = sum(hw * (basis%nq(:, a) + 0.5_dp))
        e_exact(2 * a) = e_exact(2 * a - 1)
    end do
    call sort(e_exact)
    call report("triaxial HO spectrum", maxval(abs(e - e_exact)), 1.0e-9_dp)
    deallocate(T, Vm, V, Axy, Ayz, Azx, H, C, e, e_exact)

    !---------------------------
    ! 3. spherical Woods-Saxon + spin-orbit (neutrons, A = 40)
    !---------------------------
    hw = 41.0_dp / 40.0_dp**(1.0_dp / 3.0_dp)
    call basis%initialize(hw, 10, 40)
    ns = basis%n_states
    ng = basis%n_gh
    allocate(T(ns, ns), Vm(ns, ns), Axy(ns, ns), Ayz(ns, ns), Azx(ns, ns))
    allocate(V(ng, ng, ng), dVx(ng, ng, ng), dVy(ng, ng, ng), dVz(ng, ng, ng))
    do k = 1, ng
        z = basis%point(3, k)
        do j = 1, ng
            y = basis%point(2, j)
            do i = 1, ng
                x = basis%point(1, i)
                r = sqrt(x * x + y * y + z * z)
                f = 1.0_dp / (1.0_dp + exp((r - R_ws) / a_ws))
                V(i, j, k) = - V0 * f
                dfdr = - f * (1.0_dp - f) / a_ws
                ! grad V = -V0 f'(r) r_hat ; r > 0 on all GH points for even n_gh
                dVx(i, j, k) = - V0 * dfdr * x / r
                dVy(i, j, k) = - V0 * dfdr * y / r
                dVz(i, j, k) = - V0 * dfdr * z / r
            end do
        end do
    end do
    call basis%kinetic_matrix(T)
    call basis%potential_matrix(V, Vm)
    call basis%spin_orbit_matrices(dVx, dVy, dVz, Axy, Ayz, Azx)
    kappa = so_strength(0.01875_dp * 40.0_dp + 31.5_dp)
    allocate(H(2 * ns, 2 * ns), C(2 * ns, 2 * ns), e(2 * ns))
    call basis%build_hamiltonian(T, Vm, Axy, Ayz, Azx, kappa, H)
    call diagonalize_hamiltonian(H, e, C)

    print *, "lowest levels [MeV]:"
    n_group = 1
    group_size = 0
    group_size(1) = 1
    do n = 2, 2 * ns
        if (abs(e(n) - e(n - 1)) < 1.0e-2_dp) then
            group_size(n_group) = group_size(n_group) + 1
        else
            if (n_group == size(group_size)) exit
            n_group = n_group + 1
            group_size(n_group) = 1
        end if
    end do
    a = 1
    do n = 1, size(group_size)
        if (group_size(n) == 0) exit
        print '(a, i2, a, f10.4, a, f10.6)', "  degeneracy ", group_size(n), "  e = ", e(a), &
            "  spread = ", e(a + group_size(n) - 1) - e(a)
        a = a + group_size(n)
    end do
    if (all(group_size(1:4) == expected)) then
        print *, "PASS: (2j+1) degeneracies 2,4,2,6 with p3/2 below p1/2"
    else
        print *, "FAIL: degeneracies ", group_size(1:4), " expected ", expected
        ok_all = .false.
    end if

    if (ok_all) then
        print *, "ALL TESTS PASSED"
    else
        print *, "SOME TESTS FAILED"
        error stop 1
    end if

contains
    subroutine report(name, err, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: err, tol
        if (err < tol) then
            print '(a, a, a, es10.2)', "PASS: ", name, "  max error = ", err
        else
            print '(a, a, a, es10.2)', "FAIL: ", name, "  max error = ", err
            ok_all = .false.
        end if
    end subroutine report

    subroutine sort(v)
        real(dp), intent(inout) :: v(:)
        integer :: p, q
        real(dp) :: tmp
        do p = 2, size(v)
            tmp = v(p)
            q = p - 1
            do while (q >= 1)
                if (v(q) <= tmp) exit
                v(q + 1) = v(q)
                q = q - 1
            end do
            v(q + 1) = tmp
        end do
    end subroutine sort
end program test_ho_basis
