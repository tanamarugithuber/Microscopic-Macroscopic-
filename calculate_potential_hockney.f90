module hockney_mod
    !---------------------------
    ! Free-space convolution potentials on a 3D grid by Hockney's method
    !   R. W. Hockney and J. W. Eastwood, Computer Simulation Using Particles (1988), Sec. 6-5-4
    !
    !   Coulomb : V(r) = int rho(r') / |r - r'| d^3r'                    (Laplacian V = -4 pi rho)
    !   Yukawa  : V(r) = int rho(r') exp(-|r - r'|/a) / |r - r'| d^3r'    ((-Laplacian + 1/a^2) V = 4 pi rho)
    !   exp     : W(r) = int rho(r') exp(-|r - r'|/a) d^3r'               (the non-singular part of the
    !             FRDM kernels: B_1 ~ (2a/r - 1) exp(-r/a), B_3 ~ 1/r - (1/r + 1/(2a)) exp(-r/a))
    ! with V -> 0 at infinity. These are the potentials that CG_method_poisson3D and
    ! CG_method_helmholtz_EQ3D (with 4 pi rho as source) compute with a Dirichlet condition on the box.
    !
    ! Method:
    !   The grid (nx, ny, nz) is embedded in a box of (2nx, 2ny, 2nz) filled with zeros. The kernel
    !   G(r) (exp(-kappa r)/r or exp(-kappa r)) is sampled on the doubled box at the minimum-image offsets. A cyclic
    !   convolution on the doubled box then equals the free-space (aperiodic) sum
    !       V_i = sum_j G(r_i - r_j) rho_j h^3
    !   exactly, because no two points of the original grid are further apart than half the doubled box.
    !   So there is no boundary condition to approximate; only the quadrature of the singular
    !   kernel remains.
    !
    ! Kernel quadrature ("integrated Green's function"): G at each offset is replaced by its
    ! average over the grid cell,  G_m = (1/h^3) int_cell(m) G(r) d^3r.
    !   - self cell:   2.3800773.../h for 1/r, plus Gauss-Legendre for the smooth (exp(-kappa r) - 1)/r
    !   - near cells (|m|_max <= n_near): Gauss-Legendre
    !   - far cells:   G + (h^2/24) Laplacian G  (Laplacian = kappa^2 G for exp(-kappa r)/r,
    !                  (kappa^2 - 2 kappa/r) G for exp(-kappa r))
    ! The result is then exact for a density that is constant in each cell, e.g. the
    ! inside/outside index of grid_mod used for B_1 and B_3.
    ! For a smooth density sampled at the grid points, pass smooth = .true.: the density is
    ! first corrected to cell averages, rho - (h^2/24) Laplacian rho, which makes the error O(h^4).
    !
    ! Grid: rho(i,j,k) at x = (i - (nx+1)/2) h etc. (same spacing h on all axes; nx, ny, nz may differ).
    ! Memory: 8 times the grid in real + complex arrays. Requires FFTW3.
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
    private

    integer, parameter :: dp = real64
    ! int over the unit cube [-1/2, 1/2]^3 of 1/r
    real(dp), parameter :: cube_inv_r = 2.3800773639795_dp
    integer, parameter :: n_near = 4   ! cells with |m|_max <= n_near get a Gauss-Legendre average
    integer, parameter :: kernel_yukawa = 1   ! exp(-kappa r)/r  (kappa = 0: Coulomb)
    integer, parameter :: kernel_exp = 2      ! exp(-kappa r)

    public :: hockney_coulomb, hockney_yukawa, hockney_exp

    contains

        subroutine hockney_coulomb(rho, h, V, smooth)
            ! V = int rho(r') / |r - r'| d^3r'
            implicit none
            real(dp), intent(in) :: rho(:,:,:), h
            real(dp), intent(out) :: V(:,:,:)
            logical, intent(in), optional :: smooth
            call convolve(rho, h, 0.0_dp, kernel_yukawa, V, smooth)
        end subroutine hockney_coulomb

        subroutine hockney_yukawa(rho, h, a, V, smooth)
            ! V = int rho(r') exp(-|r - r'|/a) / |r - r'| d^3r'
            implicit none
            real(dp), intent(in) :: rho(:,:,:), h, a
            real(dp), intent(out) :: V(:,:,:)
            logical, intent(in), optional :: smooth
            call convolve(rho, h, 1.0_dp / a, kernel_yukawa, V, smooth)
        end subroutine hockney_yukawa

        subroutine hockney_exp(rho, h, a, W, smooth)
            ! W = int rho(r') exp(-|r - r'|/a) d^3r'
            implicit none
            real(dp), intent(in) :: rho(:,:,:), h, a
            real(dp), intent(out) :: W(:,:,:)
            logical, intent(in), optional :: smooth
            call convolve(rho, h, 1.0_dp / a, kernel_exp, W, smooth)
        end subroutine hockney_exp

        subroutine convolve(rho, h, kappa, kernel, V, smooth)
            !---------------------------
            ! V = int rho(r') G(|r - r'|) d^3r'
            !   kernel = kernel_yukawa: G = exp(-kappa r)/r (kappa = 0: Coulomb);  kernel_exp: G = exp(-kappa r)
            ! smooth (default .false.): rho is a smooth function sampled at the grid points
            !   (.false.: rho is constant in each cell)
            ! FFTW is row-major: a Fortran array (n1, n2, n3) is passed with dims (n3, n2, n1),
            ! and the r2c output has shape (n1/2+1, n2, n3).
            !---------------------------
            !$ use omp_lib
            implicit none
            real(dp), intent(in) :: rho(:,:,:), h, kappa
            integer, intent(in) :: kernel
            real(dp), intent(out) :: V(:,:,:)
            logical, intent(in), optional :: smooth
            real(C_DOUBLE), allocatable :: work(:,:,:)
            complex(C_DOUBLE_COMPLEX), allocatable :: rho_k(:,:,:), G_k(:,:,:)
            type(C_PTR) :: plan_f, plan_b
            integer :: nx, ny, nz, n1, n2, n3, i, j, k, mx, my, mz
            real(dp) :: r, dx, dy, dz, norm
            real(dp) :: xg(16), wg(16)
            logical :: smooth_rho

            nx = size(rho, 1)
            ny = size(rho, 2)
            nz = size(rho, 3)
            if (any(shape(V) /= shape(rho))) stop "hockney_mod: V and rho must have the same shape"
            n1 = 2 * nx
            n2 = 2 * ny
            n3 = 2 * nz
            smooth_rho = .false.
            if (present(smooth)) smooth_rho = smooth

            allocate(work(n1, n2, n3), rho_k(n1 / 2 + 1, n2, n3), G_k(n1 / 2 + 1, n2, n3))
            plan_f = fftw_plan_dft_r2c_3d(n3, n2, n1, work, rho_k, FFTW_ESTIMATE)
            plan_b = fftw_plan_dft_c2r_3d(n3, n2, n1, rho_k, work, FFTW_ESTIMATE)

            !---------------------------
            ! kernel at the minimum-image offsets of the doubled box
            !---------------------------
            call gauss_legendre(16, xg, wg)
            !$omp parallel do collapse(2) default(none) private(i, j, k, mx, my, mz, dx, dy, dz, r) &
            !$omp shared(work, n1, n2, n3, h, kappa, kernel, xg, wg)
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1
                        mx = offset(i, n1)
                        my = offset(j, n2)
                        mz = offset(k, n3)
                        if (max(abs(mx), abs(my), abs(mz)) <= n_near) then
                            work(i, j, k) = cell_average(mx, my, mz, h, kappa, kernel, xg, wg)
                        else
                            dx = h * mx
                            dy = h * my
                            dz = h * mz
                            r = sqrt(dx * dx + dy * dy + dz * dz)
                            if (kernel == kernel_yukawa) then
                                work(i, j, k) = (1.0_dp + (kappa * h)**2 / 24.0_dp) * exp(-kappa * r) / r
                            else
                                work(i, j, k) = (1.0_dp + h * h / 24.0_dp * (kappa * kappa - 2.0_dp * kappa / r)) &
                                              * exp(-kappa * r)
                            end if
                        end if
                    end do
                end do
            end do
            !$omp end parallel do
            call fftw_execute_dft_r2c(plan_f, work, G_k)

            !---------------------------
            ! zero-padded density
            !---------------------------
            work = 0.0_dp
            work(1:nx, 1:ny, 1:nz) = rho
            if (smooth_rho) call to_cell_average(rho, work)
            call fftw_execute_dft_r2c(plan_f, work, rho_k)

            ! convolution theorem; h^3 for the quadrature, 1/(n1 n2 n3) for the unnormalised c2r
            norm = h**3 / (real(n1, dp) * real(n2, dp) * real(n3, dp))
            rho_k = rho_k * G_k * norm
            call fftw_execute_dft_c2r(plan_b, rho_k, work)
            V = work(1:nx, 1:ny, 1:nz)

            call fftw_destroy_plan(plan_f)
            call fftw_destroy_plan(plan_b)
            deallocate(work, rho_k, G_k)
        end subroutine convolve

        pure function offset(i, n) result(m)
            ! cyclic index i = 1..n  ->  offset in -n/2 .. n/2 - 1 (i = 1 is offset 0)
            implicit none
            integer, intent(in) :: i, n
            integer :: m
            m = i - 1
            if (m >= n / 2) m = m - n
        end function offset

        pure function cell_average(mx, my, mz, h, kappa, kernel, xg, wg) result(G)
            !---------------------------
            ! (1/h^3) int over the cell centred at h (mx, my, mz) of the kernel
            ! self cell of exp(-kappa r)/r: cube_inv_r / h + Gauss-Legendre of the smooth (exp(-kappa r) - 1)/r
            ! other cells: Gauss-Legendre, 16 points per axis next to the origin, 8 further out
            !---------------------------
            implicit none
            integer, intent(in) :: mx, my, mz, kernel
            real(dp), intent(in) :: h, kappa, xg(16), wg(16)
            real(dp) :: G
            real(dp) :: x8(8), w8(8), s, r, px, py, pz
            integer :: i, j, k, ng

            if (kernel == kernel_yukawa .and. mx == 0 .and. my == 0 .and. mz == 0) then
                G = cube_inv_r / h
                if (kappa == 0.0_dp) return
                s = 0.0_dp
                do k = 1, 16
                    do j = 1, 16
                        do i = 1, 16
                            r = 0.5_dp * h * sqrt(xg(i)**2 + xg(j)**2 + xg(k)**2)
                            s = s + wg(i) * wg(j) * wg(k) * expm1_over_r(kappa, r)
                        end do
                    end do
                end do
                G = G + s / 8.0_dp
                return
            end if

            if (max(abs(mx), abs(my), abs(mz)) <= 1) then  ! includes the self cell of exp(-kappa r)
                ng = 16
            else
                ng = 8
            end if
            call gauss_legendre8(x8, w8)
            s = 0.0_dp
            do k = 1, ng
                do j = 1, ng
                    do i = 1, ng
                        if (ng == 16) then
                            px = h * (mx + 0.5_dp * xg(i))
                            py = h * (my + 0.5_dp * xg(j))
                            pz = h * (mz + 0.5_dp * xg(k))
                            r = sqrt(px * px + py * py + pz * pz)
                            s = s + wg(i) * wg(j) * wg(k) * kernel_value(kernel, kappa, r)
                        else
                            px = h * (mx + 0.5_dp * x8(i))
                            py = h * (my + 0.5_dp * x8(j))
                            pz = h * (mz + 0.5_dp * x8(k))
                            r = sqrt(px * px + py * py + pz * pz)
                            s = s + w8(i) * w8(j) * w8(k) * kernel_value(kernel, kappa, r)
                        end if
                    end do
                end do
            end do
            G = s / 8.0_dp
        end function cell_average

        pure function kernel_value(kernel, kappa, r) result(G)
            implicit none
            integer, intent(in) :: kernel
            real(dp), intent(in) :: kappa, r
            real(dp) :: G
            if (kernel == kernel_yukawa) then
                G = exp(-kappa * r) / r
            else
                G = exp(-kappa * r)
            end if
        end function kernel_value

        subroutine to_cell_average(rho, work)
            !---------------------------
            ! cell average of a smooth density from its point values: rho - (h^2/24) Laplacian rho
            ! (7-point Laplacian, rho = 0 outside the grid), written into work(1:nx, 1:ny, 1:nz)
            !---------------------------
            implicit none
            real(dp), intent(in) :: rho(:,:,:)
            real(C_DOUBLE), intent(inout) :: work(:,:,:)
            integer :: nx, ny, nz, i, j, k
            real(dp) :: lap

            nx = size(rho, 1)
            ny = size(rho, 2)
            nz = size(rho, 3)
            !$omp parallel do collapse(2) default(none) private(i, j, k, lap) shared(rho, work, nx, ny, nz)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        lap = - 6.0_dp * rho(i, j, k)
                        if (i > 1)  lap = lap + rho(i - 1, j, k)
                        if (i < nx) lap = lap + rho(i + 1, j, k)
                        if (j > 1)  lap = lap + rho(i, j - 1, k)
                        if (j < ny) lap = lap + rho(i, j + 1, k)
                        if (k > 1)  lap = lap + rho(i, j, k - 1)
                        if (k < nz) lap = lap + rho(i, j, k + 1)
                        ! (h^2/24) * lap / h^2
                        work(i, j, k) = rho(i, j, k) - lap / 24.0_dp
                    end do
                end do
            end do
            !$omp end parallel do
        end subroutine to_cell_average

        pure function expm1_over_r(kappa, r) result(f)
            ! (exp(-kappa r) - 1) / r without cancellation for small kappa r
            implicit none
            real(dp), intent(in) :: kappa, r
            real(dp) :: f, t
            t = kappa * r
            if (t < 1.0e-3_dp) then
                f = -kappa * (1.0_dp - 0.5_dp * t + t * t / 6.0_dp - t**3 / 24.0_dp)
            else
                f = (exp(-t) - 1.0_dp) / r
            end if
        end function expm1_over_r

        pure subroutine gauss_legendre8(x, w)
            ! 8-point Gauss-Legendre nodes and weights on [-1, 1]
            implicit none
            real(dp), intent(out) :: x(8), w(8)
            x(5:8) = [0.1834346424956498_dp, 0.5255324099163290_dp, 0.7966664774136267_dp, 0.9602898564975363_dp]
            w(5:8) = [0.3626837833783620_dp, 0.3137066458778873_dp, 0.2223810344533745_dp, 0.1012285362903763_dp]
            x(1:4) = - x(8:5:-1)
            w(1:4) = w(8:5:-1)
        end subroutine gauss_legendre8

        subroutine gauss_legendre(n, x, w)
            ! nodes and weights on [-1, 1]
            implicit none
            integer, intent(in) :: n
            real(dp), intent(out) :: x(n), w(n)
            real(dp) :: z, z1, p1, p2, p3, pp
            integer :: i, j

            do i = 1, (n + 1) / 2
                z = cos(pi * (i - 0.25_dp) / (n + 0.5_dp))
                do
                    p1 = 1.0_dp
                    p2 = 0.0_dp
                    do j = 1, n
                        p3 = p2
                        p2 = p1
                        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
                    end do
                    pp = n * (z * p1 - p2) / (z * z - 1.0_dp)
                    z1 = z
                    z = z1 - p1 / pp
                    if (abs(z - z1) < 1.0e-15_dp) exit
                end do
                x(i) = -z
                x(n + 1 - i) = z
                w(i) = 2.0_dp / ((1.0_dp - z * z) * pp * pp)
                w(n + 1 - i) = w(i)
            end do
        end subroutine gauss_legendre

end module hockney_mod
