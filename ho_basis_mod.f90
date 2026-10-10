module ho_basis_mod
    !---------------------------
    ! Triaxially deformed (Cartesian) harmonic-oscillator basis
    !
    !   phi_{nx ny nz sigma}(r) = psi_nx(x; b_x) psi_ny(y; b_y) psi_nz(z; b_z) chi_sigma
    !   b_i = hbar / sqrt(m * hbar*omega_i)
    !
    ! Single-particle Hamiltonian
    !   h = T + V(r) + V_so,  V_so = i*kappa * sigma . (grad V1 x grad)
    !   kappa = lambda * (hbar/(2 m c))^2
    !
    ! Matrix elements of local operators are computed with Gauss-Hermite quadrature.
    ! On the quadrature points the 1D basis functions are stored as
    !   u_n(k) = sqrt(W_k) * hhat_n(xi_k),   W_k = w_k * exp(xi_k^2)
    !   hhat_n(xi) = H_n(xi) exp(-xi^2/2) / sqrt(2^n n! sqrt(pi))
    ! so that <n'|f(x)|n> = sum_k u_n'(k) u_n(k) f(b xi_k).
    ! u_n(k) does not depend on the axis (only xi is used); the axis enters through b_i.
    !
    ! Potentials (V, grad V1) must be given on the quadrature points
    !   (x_i, y_j, z_k) = (b_x xi_i, b_y xi_j, b_z xi_k)   -> use this%point(axis, i)
    !
    ! Units:
    !   Energy : MeV
    !   Length : fm
    !
    ! Requires LAPACK (dstev, zheev).
    !---------------------------
    !$ use omp_lib
    use iso_fortran_env, only: real64
    use constant_mod, only: pi, hbar_c, m_nucleon
    implicit none
    private

    integer, parameter :: dp = real64

    type, public :: ho_basis_type
        !---------------------------
        ! Oscillator parameters
        !---------------------------
        integer :: N_bas                ! number of (deformed) oscillator shells
        real(dp) :: hbar_omega0         ! MeV, (hw_x hw_y hw_z)^(1/3)
        real(dp) :: hbar_omega(3)       ! MeV, hw_x, hw_y, hw_z
        real(dp) :: b(3)                ! fm, oscillator lengths b_x, b_y, b_z

        !---------------------------
        ! Basis states (spatial part only, spin is added in build_hamiltonian)
        !---------------------------
        integer :: n_states             ! number of spatial basis states
        integer :: n_ax_max             ! max quantum number along any axis
        integer, allocatable :: nq(:,:) ! nq(1:3, a) = (nx, ny, nz) of state a
        integer, allocatable :: parity(:) ! (-1)**(nx+ny+nz)

        !---------------------------
        ! Gauss-Hermite quadrature
        !---------------------------
        integer :: n_gh                 ! number of quadrature points per axis
        real(dp), allocatable :: xi(:)  ! nodes (dimensionless)
        real(dp), allocatable :: u(:,:)  ! u(0:n_ax_max+1, n_gh), basis function values
        real(dp), allocatable :: du(:,:) ! du(0:n_ax_max, n_gh), d/dxi of basis functions

        contains
            procedure :: initialize => initialize_ho_basis
            procedure :: point
            procedure :: kinetic_matrix
            procedure :: potential_matrix
            procedure :: spin_orbit_matrices
            procedure :: build_hamiltonian
    end type ho_basis_type

    public :: frequencies_from_eps_gamma
    public :: gauss_hermite
    public :: so_strength
    public :: diagonalize_hamiltonian

    contains
        function frequencies_from_eps_gamma(hbar_omega0, eps, gamma) result(hw)
            !---------------------------
            ! Nilsson (eps, gamma) parametrisation
            !   omega_x = omega0 [1 - 2/3 eps cos(gamma + 2pi/3)]
            !   omega_y = omega0 [1 - 2/3 eps cos(gamma - 2pi/3)]
            !   omega_z = omega0 [1 - 2/3 eps cos(gamma)]
            ! rescaled so that hw_x hw_y hw_z = hbar_omega0**3 (volume conservation).
            ! gamma in radian.
            !---------------------------
            implicit none
            real(dp), intent(in) :: hbar_omega0, eps, gamma
            real(dp) :: hw(3)

            hw(1) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma + 2.0_dp * pi / 3.0_dp)
            hw(2) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma - 2.0_dp * pi / 3.0_dp)
            hw(3) = 1.0_dp - 2.0_dp / 3.0_dp * eps * cos(gamma)
            if (any(hw <= 0.0_dp)) stop "frequencies_from_eps_gamma: non-positive frequency"
            hw = hw * hbar_omega0 / product(hw)**(1.0_dp / 3.0_dp)
        end function frequencies_from_eps_gamma

        subroutine initialize_ho_basis(this, hbar_omega, N_bas, n_gh)
            implicit none
            class(ho_basis_type), intent(inout) :: this
            real(dp), intent(in) :: hbar_omega(3)
            integer, intent(in) :: N_bas
            integer, intent(in) :: n_gh

            integer :: nx, ny, nz, nmax(3), a, k
            real(dp) :: e_cut
            real(dp), allocatable :: W(:), hh(:)
            print *, "Initializing the deformed harmonic-oscillator basis..."

            this%N_bas = N_bas
            this%hbar_omega = hbar_omega
            this%hbar_omega0 = product(hbar_omega)**(1.0_dp / 3.0_dp)
            this%b = hbar_c / sqrt(m_nucleon * hbar_omega)

            !---------------------------
            ! basis truncation
            !   hw_x(nx+1/2) + hw_y(ny+1/2) + hw_z(nz+1/2) <= (N_bas + 3/2) hw0
            !---------------------------
            e_cut = (real(N_bas, dp) + 1.5_dp) * this%hbar_omega0 * (1.0_dp + 1.0e-12_dp)
            nmax = int((e_cut - 0.5_dp * sum(hbar_omega)) / hbar_omega)

            if (allocated(this%nq)) deallocate(this%nq)
            if (allocated(this%parity)) deallocate(this%parity)

            ! first pass: count, second pass: fill
            this%n_states = 0
            do nz = 0, nmax(3)
                do ny = 0, nmax(2)
                    do nx = 0, nmax(1)
                        if (osc_energy(nx, ny, nz) <= e_cut) this%n_states = this%n_states + 1
                    end do
                end do
            end do
            allocate(this%nq(3, this%n_states), this%parity(this%n_states))
            a = 0
            do nz = 0, nmax(3)
                do ny = 0, nmax(2)
                    do nx = 0, nmax(1)
                        if (osc_energy(nx, ny, nz) <= e_cut) then
                            a = a + 1
                            this%nq(:, a) = [nx, ny, nz]
                            this%parity(a) = 1 - 2 * mod(nx + ny + nz, 2)
                        end if
                    end do
                end do
            end do
            this%n_ax_max = maxval(this%nq)

            !---------------------------
            ! Gauss-Hermite quadrature and basis functions on the nodes
            !---------------------------
            this%n_gh = n_gh
            if (n_gh <= this%n_ax_max + 1) stop "initialize_ho_basis: n_gh must exceed n_ax_max + 1"
            if (allocated(this%xi)) deallocate(this%xi)
            if (allocated(this%u)) deallocate(this%u)
            if (allocated(this%du)) deallocate(this%du)
            allocate(this%xi(n_gh), W(n_gh), hh(0:this%n_ax_max + 1))
            allocate(this%u(0:this%n_ax_max + 1, n_gh), this%du(0:this%n_ax_max, n_gh))

            call gauss_hermite(n_gh, this%xi, W)

            do k = 1, n_gh
                call ho_functions(this%n_ax_max + 1, this%xi(k), hh)
                this%u(:, k) = sqrt(W(k)) * hh(:)
            end do

            ! d/dxi psi_n = sqrt(n/2) psi_{n-1} - sqrt((n+1)/2) psi_{n+1}
            do k = 1, n_gh
                this%du(0, k) = - sqrt(0.5_dp) * this%u(1, k)
                do nx = 1, this%n_ax_max
                    this%du(nx, k) = sqrt(0.5_dp * nx) * this%u(nx - 1, k) &
                                   - sqrt(0.5_dp * (nx + 1)) * this%u(nx + 1, k)
                end do
            end do
            deallocate(W, hh)

            print *, "hbar*omega (x,y,z) [MeV]: ", this%hbar_omega
            print *, "b (x,y,z) [fm]: ", this%b
            print *, "Number of spatial basis states: ", this%n_states
            print *, "Max quantum number per axis: ", this%n_ax_max
            print *, "Number of Gauss-Hermite points per axis: ", this%n_gh
            print *, "Deformed harmonic-oscillator basis initialized."

        contains
            pure function osc_energy(mx, my, mz) result(e)
                integer, intent(in) :: mx, my, mz
                real(dp) :: e
                e = hbar_omega(1) * (mx + 0.5_dp) + hbar_omega(2) * (my + 0.5_dp) &
                  + hbar_omega(3) * (mz + 0.5_dp)
            end function osc_energy
        end subroutine initialize_ho_basis

        pure function point(this, axis, i) result(x)
            ! coordinate (fm) of the i-th quadrature point along axis (1=x, 2=y, 3=z)
            implicit none
            class(ho_basis_type), intent(in) :: this
            integer, intent(in) :: axis, i
            real(dp) :: x
            x = this%b(axis) * this%xi(i)
        end function point

        pure subroutine ho_functions(nmax, x, hh)
            !---------------------------
            ! hh(n) = H_n(x) exp(-x^2/2) / sqrt(2^n n! sqrt(pi)),  n = 0..nmax
            ! by the stable three-term recursion
            !---------------------------
            implicit none
            integer, intent(in) :: nmax
            real(dp), intent(in) :: x
            real(dp), intent(out) :: hh(0:nmax)
            integer :: n

            hh(0) = pi**(-0.25_dp) * exp(-0.5_dp * x * x)
            if (nmax >= 1) hh(1) = sqrt(2.0_dp) * x * hh(0)
            do n = 1, nmax - 1
                hh(n + 1) = sqrt(2.0_dp / (n + 1)) * x * hh(n) - sqrt(real(n, dp) / (n + 1)) * hh(n - 1)
            end do
        end subroutine ho_functions

        subroutine gauss_hermite(n, x, W)
            !---------------------------
            ! Gauss-Hermite nodes x(k) and scaled weights W(k) = w(k) exp(x(k)^2)
            ! for the weight function exp(-x^2).
            !   nodes  : eigenvalues of the Jacobi matrix (Golub-Welsch), polished by Newton
            !   weights: Christoffel formula  W(k) = 1 / sum_{m<n} hhat_m(x_k)^2
            !---------------------------
            implicit none
            integer, intent(in) :: n
            real(dp), intent(out) :: x(n), W(n)
            real(dp), allocatable :: e(:), hh(:)
            real(dp) :: z(1,1), work(1)
            integer :: k, iter, info

            allocate(e(max(n - 1, 1)), hh(0:n))
            x = 0.0_dp
            do k = 1, n - 1
                e(k) = sqrt(0.5_dp * k)
            end do
            call dstev('N', n, x, e, z, 1, work, info)
            if (info /= 0) stop "gauss_hermite: dstev failed"

            do k = 1, n
                do iter = 1, 2
                    call ho_functions(n, x(k), hh)
                    x(k) = x(k) - hh(n) / (sqrt(2.0_dp * n) * hh(n - 1))
                end do
                call ho_functions(n, x(k), hh)
                W(k) = 1.0_dp / sum(hh(0:n - 1)**2)
            end do
            deallocate(e, hh)
        end subroutine gauss_hermite

        subroutine kinetic_matrix(this, T)
            !---------------------------
            ! <n'|p^2/2m|n> = hw/4 [ (2n+1) d_{n',n} - sqrt(n(n-1)) d_{n',n-2}
            !                        - sqrt((n+1)(n+2)) d_{n',n+2} ]   per axis
            !---------------------------
            implicit none
            class(ho_basis_type), intent(in) :: this
            real(dp), intent(out) :: T(this%n_states, this%n_states)
            integer :: a, b, ax, n, dn(3)

            T = 0.0_dp
            !$omp parallel do default(none) private(a, b, ax, n, dn) shared(this, T)
            do b = 1, this%n_states
                do a = 1, this%n_states
                    dn = this%nq(:, a) - this%nq(:, b)
                    do ax = 1, 3
                        if (count(dn /= 0) > 1) exit
                        if (dn(mod(ax, 3) + 1) /= 0 .or. dn(mod(ax + 1, 3) + 1) /= 0) cycle
                        n = this%nq(ax, b)
                        select case (dn(ax))
                        case (0)
                            T(a, b) = T(a, b) + 0.25_dp * this%hbar_omega(ax) * (2 * n + 1)
                        case (-2)
                            T(a, b) = T(a, b) - 0.25_dp * this%hbar_omega(ax) * sqrt(real(n * (n - 1), dp))
                        case (2)
                            T(a, b) = T(a, b) - 0.25_dp * this%hbar_omega(ax) * sqrt(real((n + 1) * (n + 2), dp))
                        end select
                    end do
                end do
            end do
            !$omp end parallel do
        end subroutine kinetic_matrix

        subroutine quad_matrix(this, Lx, Rx, Ly, Ry, Lz, Rz, F, M)
            !---------------------------
            ! M(a,b) = sum_{i,j,k} Lx(n'x,i) Rx(nx,i) Ly(n'y,j) Ry(ny,j) Lz(n'z,k) Rz(nz,k) F(i,j,k)
            !   (n'x,n'y,n'z) = nq(:,a), (nx,ny,nz) = nq(:,b)
            ! The sum is factorised axis by axis: z -> y -> x.
            !---------------------------
            implicit none
            class(ho_basis_type), intent(in) :: this
            real(dp), intent(in) :: Lx(0:,:), Rx(0:,:), Ly(0:,:), Ry(0:,:), Lz(0:,:), Rz(0:,:)
            real(dp), intent(in) :: F(:,:,:)
            real(dp), intent(out) :: M(this%n_states, this%n_states)
            real(dp), allocatable :: A1(:,:,:,:), A2(:,:,:,:,:)
            integer :: nm, ng, i, j, k, p, q, py, qy, a, b
            real(dp) :: s

            nm = this%n_ax_max
            ng = this%n_gh
            if (size(F, 1) /= ng .or. size(F, 2) /= ng .or. size(F, 3) /= ng) &
                stop "quad_matrix: F must be given on the n_gh^3 quadrature points"

            ! A1(i,j,p,q) = sum_k Lz(p,k) Rz(q,k) F(i,j,k)
            allocate(A1(ng, ng, 0:nm, 0:nm))
            A1 = 0.0_dp
            !$omp parallel do collapse(2) default(none) private(p, q, k) shared(A1, Lz, Rz, F, nm, ng)
            do q = 0, nm
                do p = 0, nm
                    do k = 1, ng
                        A1(:, :, p, q) = A1(:, :, p, q) + Lz(p, k) * Rz(q, k) * F(:, :, k)
                    end do
                end do
            end do
            !$omp end parallel do

            ! A2(i,py,qy,p,q) = sum_j Ly(py,j) Ry(qy,j) A1(i,j,p,q)
            allocate(A2(ng, 0:nm, 0:nm, 0:nm, 0:nm))
            A2 = 0.0_dp
            !$omp parallel do collapse(2) default(none) private(p, q, py, qy, j) &
            !$omp shared(A1, A2, Ly, Ry, nm, ng)
            do q = 0, nm
                do p = 0, nm
                    do qy = 0, nm
                        do py = 0, nm
                            do j = 1, ng
                                A2(:, py, qy, p, q) = A2(:, py, qy, p, q) + Ly(py, j) * Ry(qy, j) * A1(:, j, p, q)
                            end do
                        end do
                    end do
                end do
            end do
            !$omp end parallel do
            deallocate(A1)

            ! M(a,b) = sum_i Lx(n'x,i) Rx(nx,i) A2(i,n'y,ny,n'z,nz)
            ! TODO: skip pairs forbidden by the D2h selection rules
            !$omp parallel do default(none) private(a, b, i, s) shared(this, M, A2, Lx, Rx, ng)
            do b = 1, this%n_states
                do a = 1, this%n_states
                    s = 0.0_dp
                    do i = 1, ng
                        s = s + Lx(this%nq(1, a), i) * Rx(this%nq(1, b), i) &
                              * A2(i, this%nq(2, a), this%nq(2, b), this%nq(3, a), this%nq(3, b))
                    end do
                    M(a, b) = s
                end do
            end do
            !$omp end parallel do
            deallocate(A2)
        end subroutine quad_matrix

        subroutine potential_matrix(this, V, M)
            ! M(a,b) = <a|V|b>, V(i,j,k) given on the quadrature points
            implicit none
            class(ho_basis_type), intent(in) :: this
            real(dp), intent(in) :: V(:,:,:)
            real(dp), intent(out) :: M(this%n_states, this%n_states)
            integer :: nm

            nm = this%n_ax_max
            call quad_matrix(this, this%u(0:nm, :), this%u(0:nm, :), this%u(0:nm, :), &
                             this%u(0:nm, :), this%u(0:nm, :), this%u(0:nm, :), V, M)
            M = 0.5_dp * (M + transpose(M))
        end subroutine potential_matrix

        subroutine spin_orbit_matrices(this, dVx, dVy, dVz, Axy, Ayz, Azx)
            !---------------------------
            ! D_jk(a,b) = <a| (d_j V1) d_k |b>
            ! A_jk = D_jk - D_kj  (real antisymmetric)
            ! dVx, dVy, dVz: gradient of V1 (MeV/fm) on the quadrature points
            !---------------------------
            implicit none
            class(ho_basis_type), intent(in) :: this
            real(dp), intent(in) :: dVx(:,:,:), dVy(:,:,:), dVz(:,:,:)
            real(dp), intent(out) :: Axy(this%n_states, this%n_states)
            real(dp), intent(out) :: Ayz(this%n_states, this%n_states)
            real(dp), intent(out) :: Azx(this%n_states, this%n_states)
            real(dp), allocatable :: U(:,:), Dx(:,:), Dy(:,:), Dz(:,:), M(:,:)
            integer :: nm, ns

            nm = this%n_ax_max
            ns = this%n_states
            allocate(U(0:nm, this%n_gh), Dx(0:nm, this%n_gh), Dy(0:nm, this%n_gh), Dz(0:nm, this%n_gh))
            allocate(M(ns, ns))
            U = this%u(0:nm, :)
            Dx = this%du(0:nm, :) / this%b(1)   ! d/dx on the ket
            Dy = this%du(0:nm, :) / this%b(2)   ! d/dy on the ket
            Dz = this%du(0:nm, :) / this%b(3)   ! d/dz on the ket

            ! A_xy = <(dV/dx) d/dy> - <(dV/dy) d/dx>
            call quad_matrix(this, U, U, U, Dy, U, U, dVx, Axy)
            call quad_matrix(this, U, Dx, U, U, U, U, dVy, M)
            Axy = Axy - M

            ! A_yz = <(dV/dy) d/dz> - <(dV/dz) d/dy>
            call quad_matrix(this, U, U, U, U, U, Dz, dVy, Ayz)
            call quad_matrix(this, U, U, U, Dy, U, U, dVz, M)
            Ayz = Ayz - M

            ! A_zx = <(dV/dz) d/dx> - <(dV/dx) d/dz>
            call quad_matrix(this, U, Dx, U, U, U, U, dVz, Azx)
            call quad_matrix(this, U, U, U, U, U, Dz, dVx, M)
            Azx = Azx - M

            ! antisymmetric up to quadrature error
            Axy = 0.5_dp * (Axy - transpose(Axy))
            Ayz = 0.5_dp * (Ayz - transpose(Ayz))
            Azx = 0.5_dp * (Azx - transpose(Azx))
            deallocate(U, Dx, Dy, Dz, M)
        end subroutine spin_orbit_matrices

        pure function so_strength(lambda) result(kappa)
            ! kappa = lambda * (hbar/(2 m c))^2  [fm^2]
            implicit none
            real(dp), intent(in) :: lambda
            real(dp) :: kappa
            kappa = lambda * (hbar_c / (2.0_dp * m_nucleon))**2
        end function so_strength

        subroutine build_hamiltonian(this, T, V, Axy, Ayz, Azx, kappa, H)
            !---------------------------
            ! H in the basis |a, sigma>, index = a (sigma = up), a + n_states (sigma = down)
            !   <up|h|up> = T + V + i kappa A_xy
            !   <dn|h|dn> = T + V - i kappa A_xy
            !   <up|h|dn> = i kappa A_yz + kappa A_zx
            !   <dn|h|up> = i kappa A_yz - kappa A_zx
            ! TODO: block-diagonalise with parity and signature (D2h)
            !---------------------------
            implicit none
            class(ho_basis_type), intent(in) :: this
            real(dp), intent(in) :: T(:,:), V(:,:), Axy(:,:), Ayz(:,:), Azx(:,:)
            real(dp), intent(in) :: kappa
            complex(dp), intent(out) :: H(2 * this%n_states, 2 * this%n_states)
            complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
            integer :: ns

            ns = this%n_states
            H(1:ns, 1:ns) = cmplx(T + V, kappa * Axy, kind=dp)
            H(ns+1:2*ns, ns+1:2*ns) = cmplx(T + V, -kappa * Axy, kind=dp)
            H(1:ns, ns+1:2*ns) = iu * kappa * Ayz + kappa * Azx
            H(ns+1:2*ns, 1:ns) = iu * kappa * Ayz - kappa * Azx
            H = 0.5_dp * (H + conjg(transpose(H)))
        end subroutine build_hamiltonian

        subroutine diagonalize_hamiltonian(H, e, C)
            ! all eigenvalues e (ascending) and eigenvectors C(:,i) of a Hermitian matrix
            implicit none
            complex(dp), intent(in) :: H(:,:)
            real(dp), intent(out) :: e(:)
            complex(dp), intent(out) :: C(:,:)
            complex(dp), allocatable :: work(:)
            real(dp), allocatable :: rwork(:)
            complex(dp) :: wq(1)
            integer :: n, lwork, info

            n = size(H, 1)
            C = H
            allocate(rwork(max(1, 3 * n - 2)))
            call zheev('V', 'U', n, C, n, e, wq, -1, rwork, info)
            lwork = int(real(wq(1)))
            allocate(work(lwork))
            call zheev('V', 'U', n, C, n, e, work, lwork, rwork, info)
            if (info /= 0) stop "diagonalize_hamiltonian: zheev failed"
            deallocate(work, rwork)
        end subroutine diagonalize_hamiltonian

end module ho_basis_mod
