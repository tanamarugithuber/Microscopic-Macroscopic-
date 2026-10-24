module parameter
    implicit none
    double precision :: FRLDM_MICRO_energy
    double precision :: E_p_c_shell

    ! Common parameters used across multiple modules and programs

    double precision, parameter :: pi = 3.14159265358979323846
    double precision, parameter :: hbar_c = 197.32891d0 ! hbar*c in MeV*fm
    double precision, parameter :: hbar = 6.582119569d-22 ! hbar in MeV*s
    double precision, parameter :: dx = 1.0d-1
    double precision, parameter :: dy = 1.0d-1
    double precision, parameter :: dz = 1.0d-1
    double precision, parameter :: dV = dx * dy * dz
    double precision, parameter :: dx2 = dx*dx
    double precision, parameter :: inv_dx = 1.0d0/dx
    double precision, parameter :: inv_dy = 1.0d0/dy
    double precision, parameter :: inv_dz = 1.0d0/dz
    double precision, parameter :: inv_dx2 = 1.0d0/(dx*dx)
    double precision, parameter :: inv_dy2 = 1.0d0/(dy*dy)
    double precision, parameter :: inv_dz2 = 1.0d0/(dz*dz)
    double precision, parameter :: r_0 = 1.16d0 ! nuclear radius constant in fm

    ! Add other common parameters as needed

    !!!!!! Parameters for FRLDM calculations
    ! First category parameters
    double precision, parameter :: M_H = 7.289034d0  ! hydrogen-atom mass excess MeV/c^2
    double precision, parameter :: M_n = 8.071431d0  ! neutron mass excess MeV/c^2
    double precision, parameter :: e_2 = 1.4399764d0 ! e^2 in MeV*fm


    ! Second category parameters
    double precision, parameter :: a_Y = 0.68d0     ! range of yukawa plus exponential potential in fm
    double precision, parameter :: r_p = 0.80d0     ! proton radius constant in fm
    double precision, parameter :: a_el = 1.433d-5  ! electronic binding constant in MeV
    double precision, parameter :: a_den = 0.70d0   ! range of yukawa function used to generate  nuclear charge distribution in fm

    ! Third category parameters
    double precision, parameter :: r_mac = 4.80d0   ! average paring gap constant in MeV
    double precision, parameter :: h_np = 6.6d0     ! neutron-proton interaction constant in MeV
    double precision, parameter :: W_const = 30.0d0 ! wigner energy constant in MeV
    double precision, parameter :: a_d = 0.9d0      ! wigner damping constant in MeV

    ! Fourth category parameters
    double precision, parameter :: a_v = 16.022835d0 ! volume energy constant in MeV
    double precision, parameter :: k_v = 1.927910d0  ! volume asymmetry constant in MeV
    double precision, parameter :: a_s = 21.26946d0  ! surface energy constant in MeV
    double precision, parameter :: k_s = 2.388587d0  ! surface asymmetry constant in MeV
    double precision, parameter :: a_0 = 2.649971d0  ! a0 energy constant in MeV
    double precision, parameter :: c_a = 0.055673d0  ! charge asymmetry constant in MeV
    double precision :: B_s


    ! Coefficients calculated from the above parameters
    double precision :: I_n
    double precision :: c_1 = 0.6d0*e_2/r_0
    double precision :: c_4 = 0.8d0*(3.0d0/(2.0d0*pi))**(2.0d0/3.0d0)*0.6d0*e_2/r_0
    double precision :: f_kfr_p 
    double precision :: k_f
    double precision :: x_0 ! nucleus size parameter
    double precision :: y_0 ! nucleus deformation parameter
    double precision :: gap_p ! proton pairing gap
    double precision :: gap_n ! neutron pairing gap
    double precision :: gap_np ! neutron-proton pairing gap


    ! FRLDM energy components
    double precision :: M_ZN              ! mass excess of Z hydrogens and N neutrons
    double precision :: VE                ! volume energy
    double precision :: SE                ! surface energy
    double precision :: A0E               ! a0 energy
    double precision :: CE                ! coulomb energy
    double precision :: CEE               ! coulomb exchange energy
    double precision :: PF2CE             ! proton form factor correction to the coulomb energy
    double precision :: CAE               ! charge asymmetry energy
    double precision :: WE                ! Wigner energy
    double precision :: APE               ! average pairing energy
    double precision :: EoBE              ! energy of bound electrons
    double precision :: FRLDM_energy      ! total FRLDM energy
    double precision :: FRLDM_energy_ball ! total FRLDM energy for spherical shape

    !!!!!! Parameters for Microscopic Method

    ! First category parameter
    double precision, parameter :: m_nuc = 938.90595d0 ! nucleon mass in MeV

    ! Second category parameter
    double precision, parameter :: C_cur = 41d0      ! basis curvature constant in MeV
    double precision, parameter :: V_s = 52.5d0      ! symmetric potential-depth constant in MeV
    double precision, parameter :: V_a = 48.7d0      ! asymmetric potential-depth constant in MeV
    double precision, parameter :: bigA_den = 0.82d0 ! potential radius correction constant in fm
    double precision, parameter :: bigB_den = 0.56d0 ! potential radius curvature-correction constant in fm**2
    double precision, parameter :: a_pot = 0.8d0     ! potential diffuseness constant in fm
    double precision, parameter :: k_p = 0.025d0     ! proton spin–orbit A coefficient 
    double precision, parameter :: l_p = 28.0d0      ! proton spin–orbit constant
    double precision, parameter :: k_n = 0.01875d0   ! neutron spin–orbit A coefficient
    double precision, parameter :: l_n = 31.5d0      !  neutron spin–orbit constant

    ! Third category parameter 
    double precision, parameter :: N_bas = 12.0d0 ! number of basis functions
    double precision, parameter :: p = 8.0d0      ! order of Strutinsky shell correction
    double precision, parameter :: C_s = 1.0d0     ! Strutinsky range coefficient

    ! Fourth category parameter
    double precision, parameter :: r_mic = 3.2d0  ! LN effective-interaction pairing-gap constant
    double precision, parameter :: frK = 0.2475d0 ! Zero-point energy constant
    double precision, parameter :: a_2 = 22.00d0  ! surface-energy constant in MeV
    double precision, parameter :: bigJ = 35.0d0  ! symmetry-energy constant in MeV
    double precision, parameter :: bigL = 99.0d0  ! density-symmetry constant in MeV
    double precision, parameter :: bigQ = 25.0d0  ! effective surface-stiffness constant in MeV
    double precision, parameter :: bigK = 300.0d0 ! compressibility constant in MeV

    ! single particle energy 
    character(len=100) :: potential_type = "Yukawa" ! = "Yukawa" or "Woods-Saxon"
    double precision, allocatable :: V_particle(:,:,:)
    double precision, allocatable :: V_1(:,:,:)
    double precision, allocatable :: nabla_V_1_x(:,:,:)
    double precision, allocatable :: nabla_V_1_y(:,:,:)
    double precision, allocatable :: nabla_V_1_z(:,:,:)
    double precision, allocatable :: V_Coulomb(:,:,:)
    complex(kind(0d0)), allocatable :: V_so_up(:,:,:,:)
    complex(kind(0d0)), allocatable :: V_so_dn(:,:,:,:)
    double precision, allocatable :: E_single_particle(:)
    double precision, allocatable :: E_single_particle_proton(:)
    double precision, allocatable :: E_single_particle_neutron(:)
    double precision, allocatable :: delta_Etot(:)
    double precision, allocatable :: Resi_norm(:)
    double precision, allocatable :: norm(:)
    integer, allocatable :: active(:)
    integer :: n_active
    double precision :: V_0 
    double precision :: R_pot
    double precision :: R_den
    double precision :: V_p
    double precision :: V_n
    double precision :: epsilon_bar
    double precision :: gap_bar
    double precision :: lambda_p
    double precision :: lambda_n
    double precision :: rho_c
    double precision :: omega_0

    ! Imaginary Time Evolution
    ! wavefunctions of protons
    character(len=100) :: identifier ! = "proton" ! "proton" or "neutron"
    complex(kind(0d0)), allocatable :: psi_up(:,:,:,:),h_psi_up(:,:,:,:), psi_dn(:,:,:,:),h_psi_dn(:,:,:,:)
    complex(kind(0d0)), allocatable :: Laplacian_psi_up(:,:,:,:), Laplacian_psi_dn(:,:,:,:)
    complex(kind(0d0)), allocatable :: dpsi_dx_up(:,:,:,:), dpsi_dy_up(:,:,:,:), dpsi_dz_up(:,:,:,:)
    complex(kind(0d0)), allocatable :: dpsi_dx_dn(:,:,:,:), dpsi_dy_dn(:,:,:,:), dpsi_dz_dn(:,:,:,:)

    double precision, parameter :: delta_tau = 0.0001d0 ! imaginary time step
    double precision, allocatable :: norm_factor(:)
    integer, parameter :: max_iter = 100000
    integer :: iter
    integer :: i_x_min 
    integer :: i_x_max
    integer :: i_y_min
    integer :: i_y_max
    integer :: i_z_min
    integer :: i_z_max
    complex(kind(0d0)), parameter :: I_cmplx = cmplx(0.0d0, 1.0d0, kind(0d0))

    ! Lipkin-Nogami
    integer :: N_tot
    integer, parameter :: N_extra =20 ! number of extra levels for pairing calculation
    integer :: n_levels
    double precision, allocatable :: eps(:)
    integer, allocatable :: g(:)
    double precision, allocatable :: v2(:)
    double precision, allocatable :: u2(:)
    double precision, allocatable :: Ek(:)
    double precision :: hbar_omega_0 ! cutoff energy in MeV
    integer :: k_lo, k_hi
    integer :: N_1
    integer :: N_2
    double precision :: bigG
    double precision :: v_k2
    double precision :: epsilon_k
    double precision :: e_n_k  
    double precision :: lambda_2
    double precision :: u_k2
    double precision :: bigE_k
    double precision :: E_p_c
    double precision :: n_k
    double precision :: theta_odd
    double precision :: delta
    double precision :: delta_E_p_c_proton
    double precision :: delta_E_p_c_neutron

    ! effective-interaction paring gap
    double precision :: gap_Gn
    double precision :: gap_Gp
    double precision :: rho_tilde
    double precision :: y_1
    double precision :: y_2
    double precision :: bigG_eipg
    double precision :: E_p_c_tilde
    double precision :: lambda_wave
    double precision :: bigA
    double precision :: bigB
    double precision :: bigC
    double precision :: lambda_tilde_2

    !!!!!! Strutinsky shell correction parameters
    double precision :: g_s
    double precision :: E_Fermi      ! Fermi energy
    double precision :: E_Fermi_tilde ! smoothed Fermi energy
    double precision, parameter :: delta_E = 0.01d0  ! smoothing parameter for Strutinsky shell correction in MeV
    ! double precision, allocatable :: g_shell(:)        ! shell density of single particle states
    double precision :: gamma_shell                     ! Parameter for shell correction lamda*A_num**(-2.0d0/3.0d0)
    double precision :: lambda_shell                  ! fermi level for shell correction
    ! double precision, allocatable :: g_tilde_shell(:)      ! smoothed level density lamda*A_num**(-1/3)
    ! double precision, allocatable :: delta_g_shell(:)  ! difference between shell level density and smoothed level density
    ! if delta_g_shell is negative, density of states is low
    double precision :: E_shell         ! sum of single particle energy levels under Fermi level
    double precision :: E_tilde_shell   ! sum of smoothed single particle energy levels under Smoothed Fermi level
    double precision :: delta_E_shell_proton
    double precision :: delta_E_shell_neutron
       ! difference between single particle energy levels and smoothed single particle energy levels

    ! zero point energy
    double precision :: E_0_lambda
    double precision :: omega_lambda
    double precision :: B_lambda
    double precision :: B_lambda_irr
    double precision :: omega_lambda_irr
    double precision :: B_ep2_irr
    double precision :: B_gamma_irr
    double precision :: C_gamma
    double precision :: C_ep2

    !!!!!! Nuclear numbers
    integer, parameter :: Z_num = 8
    integer, parameter :: N_num = 8

    integer :: A_num 

    !!!!!! calculate shape parameters
    double precision :: semi_a            !spheroidal semi axis
    double precision :: semi_a_pot        !spheroidal semi axis for potential
    double precision :: semi_b            !spheroidal semi axis
    double precision :: semi_b_pot        !spheroidal semi axis for potential
    double precision :: ecc               !spheroidal eccentricity
    double precision :: obrateness        !spheroidal oblateness
    double precision :: surface_area      ! nuclear surface area
    double precision :: volume_nucleus    ! nuclear volume
    double precision :: volume_potential  ! volume used for potential calcul
    integer, allocatable :: i_shape(:,:,:)
    integer, allocatable :: i_shape_pot(:,:,:)
    integer :: j_x_max, j_y_max, j_z_max  ! max indices for nuclear shape
    integer :: j_x_min, j_y_min, j_z_min  ! min indices for nuclear shape


    !!!!!! Exponential function arrays
    double precision, allocatable :: exp_a_Y(:, :, :)
    double precision, allocatable :: exp_a_den(:, :, :)
    double precision, allocatable :: exp_a_pot(:, :, :)
    double precision, allocatable :: r_relative(:, :, :)         ! distance between two points
    double precision, allocatable :: r_relative_inverse(:, :, :) ! inverse of distance between two points



    !!!!!! Nuclear radius
    double precision :: R_A

end module 

program potential
    !$ use omp_lib
    use parameter
    implicit none
    ! integer :: ios,unit,ios2,i
    ! character(len=200) :: line
    ! integer :: state
    ! double precision :: energy
    ! character(len=20) :: dummy
    call compute_constants
    call allocate_arrays_for_FRLDM
    call compute_shape_Frequently_used_functions
    call compute_FRLDM_only_for_spheroidal
    ! call compute_V_1_and_V_cloulomb
    ! call compute_nabla_on_V_1
    ! identifier = "proton"
    ! call allocate_arrays_for_SPE
    ! call prepare_initial_wavefunction
    ! call Imaginary_Time_Evolution
    ! allocate( E_single_particle_proton( Z_num + N_extra ) )
    ! allocate( eps( Z_num + N_extra ) )
    ! unit = 10
    ! open(unit,file="Energy_Levels_protonof_8-8y.dat",status="old")
    ! do state = 1, Z_num + N_extra
    !     read(unit,*) energy
    !     E_single_particle_proton(state) = energy
    !     print *, "Proton State ", state, ": Energy =", energy
    ! end do

    ! close(unit)
    ! do i = 1, 10
    !     gamma_shell = i 
    !     print *, "Gamma_shell =", gamma_shell
    !     call compute_smoothed_fermi_energy
    ! end do
    ! gamma_shell = 6.0d0
    ! call compute_smoothed_fermi_energy
    ! call effective_BCS_pairing
    ! call shell_and_BCS
    ! call deallocate_arrays_for_SPE
    ! identifier = "neutron"
    ! call allocate_arrays_for_SPE
    ! call prepare_initial_wavefunction
    ! call Imaginary_Time_Evolution
    ! call compute_smoothed_fermi_energy
    ! call effective_BCS_pairing
    ! call shell_and_BCS

    ! call sum_FRLDM_MICRO

    

    ! Deallocate arrays
    ! deallocate( exp_a_Y )
    ! deallocate( exp_a_den )
    ! deallocate( r_relative )
    ! deallocate( r_relative_inverse )
    ! deallocate( exp_a_pot )
    ! call deallocate_arrays
end program potential

!!!!!!!!!!!!!!!!!!!!!!!!preparation subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! to allocate arrays for single particle potential calculation
subroutine allocate_arrays_for_FRLDM ! 2nd subroutine
    use parameter
    implicit none
    integer :: j_num

    allocate( V_1( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max ) )
    allocate( V_coulomb( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max ) )
    allocate( nabla_V_1_x( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max ) )
    allocate( nabla_V_1_y( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max ) )
    allocate( nabla_V_1_z( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max ) )
    allocate( exp_a_pot( 2*i_x_min:2*i_x_max, 2*i_y_min:2*i_y_max, 2*i_z_min:2*i_z_max ) )
    allocate( r_relative( 2*i_x_min:2*i_x_max, 2*i_y_min:2*i_y_max, 2*i_z_min:2*i_z_max ) )
    allocate( r_relative_inverse( 2*i_x_min:2*i_x_max, 2*i_y_min:2*i_y_max, 2*i_z_min:2*i_z_max ) )
    allocate( i_shape( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max ) )
    allocate( i_shape_pot( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max ) )
    allocate( exp_a_Y( 2*i_x_min:2*i_x_max, 2*i_y_min:2*i_y_max, 2*i_z_min:2*i_z_max ) )
    allocate( exp_a_den( 2*i_x_min:2*i_x_max, 2*i_y_min:2*i_y_max, 2*i_z_min:2*i_z_max ) )

    print *, "Arrays for single particle potential calculation allocated."
end subroutine allocate_arrays_for_FRLDM

subroutine allocate_arrays_for_SPE
    use parameter
    implicit none
    integer :: j_num
    if (identifier == "proton" ) then
        j_num = Z_num + N_extra
    else if (identifier == "neutron") then
        j_num = N_num + N_extra
    else
        print *, "Error: identifier must be 'proton' or 'neutron'."
        stop
    end if
    allocate( norm_factor( j_num ) )
    allocate( V_so_up( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate( V_so_dn( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate( psi_up( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate( h_psi_up( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate( psi_dn( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate( h_psi_dn( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate( Laplacian_psi_up( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max,j_num ) )
    allocate( Laplacian_psi_dn( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max,j_num ) )
    allocate ( dpsi_dx_up( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate ( dpsi_dy_up( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate ( dpsi_dz_up( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate ( dpsi_dx_dn( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate ( dpsi_dy_dn( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate ( dpsi_dz_dn( i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max, j_num ) )
    allocate( E_single_particle( j_num ) )
end subroutine allocate_arrays_for_SPE

subroutine deallocate_arrays_for_SPE ! deallocation subroutine
    use parameter
    implicit none

    deallocate( psi_up )
    deallocate( h_psi_up )
    deallocate( psi_dn )
    deallocate( h_psi_dn )
    deallocate( Laplacian_psi_up )
    deallocate( Laplacian_psi_dn )
    deallocate ( dpsi_dx_up )
    deallocate ( dpsi_dy_up )
    deallocate ( dpsi_dz_up )
    deallocate ( dpsi_dx_dn )
    deallocate ( dpsi_dy_dn )
    deallocate ( dpsi_dz_dn )
    deallocate( norm_factor )
    deallocate( E_single_particle )
    deallocate( V_so_up )
    deallocate( V_so_dn )
    print *, "Arrays for single particle potential calculation deallocated."
end subroutine deallocate_arrays_for_SPE


subroutine deallocate_arrays ! deallocation subroutine
    use parameter
    implicit none

    deallocate( V_1 )
    deallocate( V_coulomb )
    deallocate( V_so_up )
    deallocate( V_so_dn )
    deallocate( nabla_V_1_x )
    deallocate( nabla_V_1_y )
    deallocate( nabla_V_1_z )
    deallocate( exp_a_pot )
    deallocate( r_relative )
    deallocate( r_relative_inverse )
    deallocate( i_shape )

    deallocate( psi_up )
    deallocate( h_psi_up )
    deallocate( psi_dn )
    deallocate( h_psi_dn )
    deallocate( Laplacian_psi_up )
    deallocate( Laplacian_psi_dn )
    deallocate ( dpsi_dx_up )
    deallocate ( dpsi_dy_up )
    deallocate ( dpsi_dz_up )
    deallocate ( dpsi_dx_dn )
    deallocate ( dpsi_dy_dn )
    deallocate ( dpsi_dz_dn )
    deallocate( norm_factor )
    deallocate( E_single_particle )
    print *, "Arrays for single particle potential calculation deallocated."
end subroutine deallocate_arrays

subroutine compute_shape_Frequently_used_functions
    !$ use omp_lib
    use parameter
    implicit none
    integer :: j_x, j_y, j_z
    double precision :: x, y, z,r
    !$omp parallel default(none) &
    !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
    !$omp        i_shape, exp_a_pot, r_relative_inverse, r_relative,&
    !$omp        semi_a,semi_b,exp_a_Y, exp_a_den, semi_a_pot, semi_b_pot, &
    !$omp        i_shape_pot)&
    !$omp private(j_x,j_y,j_z,x,y,z,r)
    !$omp do collapse(3) schedule(static) 
    do j_z = 2*i_z_min, 2*i_z_max
        do j_y = 2*i_y_min, 2*i_y_max
            do j_x = 2*i_x_min, 2*i_x_max
                x = dble(j_x)*dx
                y = dble(j_y)*dy
                z = dble(j_z)*dz
                if (j_x==0 .and. j_y==0 .and. j_z==0) then
                    r_relative(j_x,j_y,j_z) = 1.0d-10 ! to avoid division by zero
                    r_relative_inverse(j_x,j_y,j_z) = 1.0d10
                    exp_a_pot(j_x,j_y,j_z) = exp(- r_relative(j_x,j_y,j_z)/ a_pot)
                else
                r = sqrt( x**2 + y**2 + z**2 )
                r_relative(j_x,j_y,j_z) = r
                r_relative_inverse(j_x,j_y,j_z) = 1.0d0 / r
                exp_a_Y(j_x,j_y,j_z) = exp(- r/ a_Y)
                exp_a_den(j_x,j_y,j_z) = exp(- r / a_den)
                ! exp_a_pot(j_x,j_y,j_z) = exp(- r/ a_pot)
                end if
                if(j_x>i_x_max .or. j_x<i_x_min .or. &
                   j_y>i_y_max .or. j_y<i_y_min .or. &
                   j_z>i_z_max .or. j_z<i_z_min) then
                    cycle
                else if ( (x**2)/(semi_a**2) + (y**2)/(semi_b**2) + (z**2)/(semi_b**2) <= 1.0d0 ) then
                    i_shape(j_x,j_y,j_z) = 1 ! inside nucleus
                else

                    i_shape(j_x,j_y,j_z) = 0 ! outside nucleus
                end if
                if(j_x>i_x_max .or. j_x<i_x_min .or. &
                   j_y>i_y_max .or. j_y<i_y_min .or. &
                   j_z>i_z_max .or. j_z<i_z_min) then
                    cycle
                else if ( (x**2)/(semi_a_pot**2) + (y**2)/(semi_b_pot**2) + (z**2)/(semi_b_pot**2) <= 1.0d0 ) then
                    i_shape_pot(j_x,j_y,j_z) = 1 ! inside potential volume
                else
                    i_shape_pot(j_x,j_y,j_z) = 0 ! outside potential volume
                end if
            end do  
        end do
    end do
    !$omp end do
    !$omp end parallel
    print *, "sum(i_shape)=", sum(i_shape), "  max(i_shape)=", maxval(i_shape)
    print *, "sum(i_shape_pot)=", sum(i_shape_pot), "  max(i_shape_pot)=", maxval(i_shape_pot)

    print *, "Nuclear shape and frequently used functions computed."

end subroutine compute_shape_Frequently_used_functions

subroutine compute_constants ! 1st subroutine
    use parameter
    implicit none
    A_num = Z_num + N_num ! mass number

    R_A = r_0 * (dble(A_num))**(1.0d0/3.0d0) ! average nuclear radius in fm
    semi_a = R_A !spheroidal semi axis
    volume_nucleus = (4.0d0/3.0d0) * pi *R_A**3 ! nuclear volume
    semi_b = sqrt( volume_nucleus / ( (4.0d0/3.0d0) * pi * semi_a ) ) !spheroidal semi axis


    ecc = sqrt( abs(1.0d0 - (min(semi_a,semi_b)**2)/(max(semi_a,semi_b)**2)) ) !spheroidal eccentricity
    obrateness = 1.0d0 - semi_b/semi_a !spheroidal oblateness


    ! calculate surface area
    if (semi_a>semi_b)then  
        surface_area = 2.0d0 * pi *( semi_b**2 + semi_a*semi_b*asin(ecc)/ ecc )
    else if (semi_a==semi_b) then
        surface_area = 4.0d0 * pi * R_A**2
        B_s = A_num**(-2.0d0/3.0d0)*surface_area/(4.0d0*pi*r_0**2)
    else
        surface_area = 2.0d0 * pi * (semi_b**2 + semi_a**2 *atanh(ecc)/ ecc )
        B_s = A_num**(-2.0d0/3.0d0)*surface_area/(4.0d0*pi*r_0**2)
    end if

    j_x_max = int( semi_a/ dx )
    j_y_max = int( semi_b/ dy )
    j_z_max = int( semi_b/ dz )

    j_x_min = -j_x_max
    j_y_min = -j_y_max
    j_z_min = -j_z_max

    i_x_max = int(2*j_x_max)
    i_x_min = -i_x_max
    i_y_max = int(2*j_y_max)
    i_y_min = -i_y_max
    i_z_max = int(2*j_z_max)
    i_z_min = -i_z_max

    A_num = Z_num + N_num ! mass number
    I_n = dble(N_num - Z_num)/dble(A_num)
    gap_bar = (I_n + 3.0d0*c_1*Z_num**2/(8.0d0*bigQ*A_num**(5.0d0/3.0d0)) ) / ( 1.0d0 + &
    9.0d0*bigJ/(4.0d0*bigQ*A_num**(1.0d0/3.0d0)) )

    epsilon_bar = (-2.0d0*a_2/(A_num**(1.0d0/3.0d0)) + bigL*gap_bar**2 + c_1*Z_num**2/(A_num**(4.0d0/3.0d0)) ) / bigK
    R_den = r_0*A_num**(1.0d0/3.0d0) * ( 1.0d0 + epsilon_bar )
    R_pot = R_den + bigA_den - bigB_den / ( R_den )
    rho_c = Z_num * sqrt(e_2) / ( (4.0d0/3.0d0)*pi*R_pot**3 )
    
    V_p = V_s + V_a * gap_bar ! potential depth for protons
    V_n = V_s - V_a * gap_bar ! potential depth for neutrons

    lambda_p = k_p * dble(A_num) + l_p ! spin–orbit coupling constant for protons
    lambda_n = k_n * dble(A_num) + l_n ! spin–orbit coupling constant for neutrons
    gamma_shell = C_s * 41.0d0 / ( A_num**(1.0d0/3.0d0) ) ! smoothing parameter for Strutinsky shell correction
    gap_Gn = r_mic* B_s * N_num**(-1.0d0/3.0d0) ! effective-interaction pairing gap for neutrons
    gap_Gp = r_mic* B_s * Z_num**(-1.0d0/3.0d0) ! effective-interaction pairing gap for protons
    volume_potential = (4.0d0/3.0d0) * pi * R_pot**3 ! volume used for potential calculation
    semi_a_pot = R_pot !spheroidal semi axis for potential
    semi_b_pot = sqrt( volume_potential / ( (4.0d0/3.0d0) * pi * semi_a_pot ) ) !spheroidal semi axis for potential
    print *, "R_A =", R_A, " fm"
    print *, "R_pot =", R_pot, " fm"
    print *, "Constants for single particle potential calculation computed."
end subroutine compute_constants

!!!!!!!!!!!!!!!!!!!! Single Particle Potential Calculation Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_V_1_and_V_cloulomb
    !$ use omp_lib
    use parameter
    implicit none
    character(len=256) :: filename
    integer :: j_x, j_y, j_z,j_x_prime, j_y_prime, j_z_prime, j_x2x, j_y2y, j_z2z
    double precision :: v1_local, vc_local,exp_r
    print *, "Computing single particle potentials V_1 and V_coulomb ..."
    V_1(:,:,:) = 0.0d0
    V_coulomb(:,:,:) = 0.0d0
    if(potential_type == "Yukawa") then
    !$omp parallel do collapse(3) default(none) schedule(static) &
    !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
    !$omp        j_x_min,j_x_max,j_y_min,j_y_max,j_z_min,j_z_max, &
    !$omp        i_shape_pot, exp_a_pot, r_relative_inverse, V_1, V_coulomb,V_0,identifier,potential_type) &
    !$omp private(j_x,j_y,j_z,j_x_prime,j_y_prime,j_z_prime,j_x2x,j_y2y,j_z2z,v1_local,vc_local)
    do j_z = i_z_min, i_z_max
        do j_y = i_y_min, i_y_max
                do j_x = i_x_min, i_x_max
                    v1_local = 0.0d0
                    vc_local = 0.0d0
                        do j_z_prime = j_z_min, j_z_max                           
                            do j_y_prime = j_y_min, j_y_max
                                    do j_x_prime = j_x_min, j_x_max
                                    if (i_shape_pot(j_x_prime,j_y_prime,j_z_prime) == 0) cycle ! outside the nucleus
                                    j_x2x = j_x - j_x_prime
                                    j_y2y = j_y - j_y_prime
                                    j_z2z = j_z - j_z_prime
                                    if (j_x2x == 0 .and. j_y2y == 0 .and. j_z2z == 0) cycle ! to avoid singularity
                                    
                                    v1_local = v1_local + &
                                    exp_a_pot(j_x2x,j_y2y,j_z2z)*r_relative_inverse(j_x2x,j_y2y,j_z2z)
                               
                                    vc_local = vc_local + &
                                    r_relative_inverse(j_x2x,j_y2y,j_z2z)
                                    end do
                            end do
                        end do
                    
                            V_1(j_x,j_y,j_z) = v1_local
                        
                        V_coulomb(j_x,j_y,j_z) = vc_local

                end do
        end do
        ! write (*,*) "Progress in calculation: j_z=", j_z, "/", i_z_max
    end do
    !$omp end parallel do

        V_1(:,:,:) = - 1/(4.0d0 * pi * a_pot**2) * V_1(:,:,:)*dV
    else if (potential_type == "Woods-Saxon") then
        V_0 = (-51.0d0 + 33.0d0 * (N_num-Z_num)/A_num ) ! potential depth in MeV for Woods-Saxon potential
        print *, "Using Woods-Saxon potential."
        do j_z = i_z_min, i_z_max
            do j_y = i_y_min, i_y_max
                do j_x = i_x_min, i_x_max
                    exp_r = exp( ( sqrt( r_relative(j_x,j_y,j_z) ) - R_A ) / 0.6d0 )
                    if ( exp_r > 700.0d0 ) then
                        V_1(j_x,j_y,j_z) = 0.0d0
                        cycle
                    end if
                    V_1(j_x,j_y,j_z) = - 1 / ( 1.0d0 + exp_r )*V_0
                end do
            end do
        end do
    end if
    print *, "min/max V1 :", minval(V_1), maxval(V_1)
    V_coulomb(:,:,:) = sqrt(e_2)*rho_c* V_coulomb(:,:,:)*dV
    print *, "min/max Vc :", minval(V_coulomb), maxval(V_coulomb)
    ! output V_1 and V_coulomb for checking
    open(unit=10, file="temp_V_1.dat", status='unknown')

    do j_x = i_x_min, i_x_max
        write(10,*) j_x*dx, V_1(j_x,0,0)*V_p
    end do
    close(10)
    ! open(unit=20, file="temp_V_coulomb.dat", status='unknown')

    ! do j_x = i_x_min, i_x_max
    !     write(20,*)j_x*dx, V_coulomb(j_x,0,0)
    ! end do

    ! close(20)

    print *, "Single particle potentials V_1 and V_coulomb computed."



end subroutine compute_V_1_and_V_cloulomb



subroutine compute_nabla_on_V_1
    !$ use omp_lib
    use parameter
    implicit none
    integer :: i_x, i_y, i_z
    double precision :: inv_840
    inv_840 = 1.0d0 / 840.0d0
    ! Compute first derivatives using central difference
    !$omp parallel default(none) &
    !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
    !$omp        V_1, nabla_V_1_x, nabla_V_1_y, nabla_V_1_z,inv_840) &
    !$omp private(i_x,i_y,i_z)

    !$omp do schedule(static) 
    do i_x = i_x_min, i_x_max !9点差分
        if (i_x == i_x_min  )then
            nabla_V_1_x(i_x,:,:) = (-3*V_1(i_x+4,:,:) + 32*V_1(i_x+3,:,:) - 168*V_1(i_x+2,:,:) &
            + 672*V_1(i_x+1,:,:) ) * inv_840 * inv_dx
        else if (i_x == i_x_min +  1) then
            nabla_V_1_x(i_x,:,:) = (-3*V_1(i_x+4,:,:) + 32*V_1(i_x+3,:,:) - 168*V_1(i_x+2,:,:) &
            + 672*V_1(i_x +1,:,:) -672*V_1(i_x -1,:,:)) * inv_840 * inv_dx
        else if (i_x == i_x_min + 2) then
            nabla_V_1_x(i_x,:,:) = (-3*V_1(i_x+4,:,:) + 32*V_1(i_x+3,:,:) - 168*V_1(i_x+2,:,:) &
            + 672*V_1(i_x+1,:,:) - 672*V_1(i_x-1,:,:) + 168*V_1(i_x-2,:,:)) * inv_840 * inv_dx
        else if (i_x == i_x_min + 3) then
            nabla_V_1_x(i_x,:,:) = (-3*V_1(i_x+4,:,:) + 32*V_1(i_x+3,:,:) - 168*V_1(i_x+2,:,:) &
            + 672*V_1(i_x+1,:,:) - 672*V_1(i_x-1,:,:) + 168*V_1(i_x-2,:,:) - 32*V_1(i_x-3,:,:)) * inv_840 * inv_dx
        else if (i_x == i_x_max - 3) then
            nabla_V_1_x(i_x,:,:) = (3*V_1(i_x-4,:,:) - 32*V_1(i_x-3,:,:) + 168*V_1(i_x-2,:,:) &
            - 672*V_1(i_x-1,:,:) + 672*V_1(i_x+1,:,:) - 168*V_1(i_x+2,:,:) + 32*V_1(i_x+3,:,:)) * inv_840 * inv_dx
        else if (i_x == i_x_max - 2) then
            nabla_V_1_x(i_x,:,:) = (3*V_1(i_x-4,:,:) - 32*V_1(i_x-3,:,:) + 168*V_1(i_x-2,:,:) &
            - 672*V_1(i_x-1,:,:) + 672*V_1(i_x+1,:,:) - 168*V_1(i_x+2,:,:) ) * inv_840 * inv_dx
        else if(i_x == i_x_max - 1) then
            nabla_V_1_x(i_x,:,:) = (3*V_1(i_x-4,:,:) - 32*V_1(i_x-3,:,:) + 168*V_1(i_x-2,:,:) &
            - 672*V_1(i_x-1,:,:) + 672*V_1(i_x+1,:,:) ) * inv_840 * inv_dx
        else if(i_x == i_x_max) then
            nabla_V_1_x(i_x,:,:) = (3*V_1(i_x-4,:,:) - 32*V_1(i_x-3,:,:) + 168*V_1(i_x-2,:,:) &
            - 672*V_1(i_x-1,:,:) ) * inv_840 * inv_dx
        else
        nabla_V_1_x(i_x,:,:) = (3*V_1(i_x-4,:,:) - 32*V_1(i_x-3,:,:) + 168*V_1(i_x-2,:,:) &
        - 672*V_1(i_x-1,:,:) + 672*V_1(i_x+1,:,:) - 168*V_1(i_x+2,:,:) + 32*V_1(i_x+3,:,:) - 3*V_1(i_x+4,:,:)) * inv_840 * inv_dx
        end if
    end do
    !$omp end do
    !$omp do schedule(static)
    do i_y = i_y_min, i_y_max
        if (i_y == i_y_min  )then
            nabla_V_1_y(:,i_y,:) = (-3*V_1(:,i_y+4,:) + 32*V_1(:,i_y+3,:) - 168*V_1(:,i_y+2,:) &
            + 672*V_1(:,i_y+1,:) ) * inv_840 * inv_dy
        else if (i_y == i_y_min +  1) then
            nabla_V_1_y(:,i_y,:) = (-3*V_1(:,i_y+4,:) + 32*V_1(:,i_y+3,:) - 168*V_1(:,i_y+2,:) &
            + 672*V_1(:,i_y +1,:) -672*V_1(:,i_y -1,:)) * inv_840 * inv_dy
        else if (i_y == i_y_min + 2) then
            nabla_V_1_y(:,i_y,:) = (-3*V_1(:,i_y+4,:) + 32*V_1(:,i_y+3,:) - 168*V_1(:,i_y+2,:) &
            + 672*V_1(:,i_y+1,:) - 672*V_1(:,i_y-1,:) + 168*V_1(:,i_y-2,:)) * inv_840 * inv_dy
        else if (i_y == i_y_min + 3) then
            nabla_V_1_y(:,i_y,:) = (-3*V_1(:,i_y+4,:) + 32*V_1(:,i_y+3,:) - 168*V_1(:,i_y+2,:) &
            + 672*V_1(:,i_y+1,:) - 672*V_1(:,i_y-1,:) + 168*V_1(:,i_y-2,:) - 32*V_1(:,i_y-3,:)) * inv_840 * inv_dy
        else if (i_y == i_y_max - 3) then
            nabla_V_1_y(:,i_y,:) = (3*V_1(:,i_y-4,:) - 32*V_1(:,i_y-3,:) + 168*V_1(:,i_y-2,:) &
            - 672*V_1(:,i_y-1,:) + 672*V_1(:,i_y+1,:) - 168*V_1(:,i_y+2,:) + 32*V_1(:,i_y+3,:)) * inv_840 * inv_dy
        else if (i_y == i_y_max - 2) then
            nabla_V_1_y(:,i_y,:) = (3*V_1(:,i_y-4,:) - 32*V_1(:,i_y-3,:) + 168*V_1(:,i_y-2,:) &
            - 672*V_1(:,i_y-1,:) + 672*V_1(:,i_y+1,:) - 168*V_1(:,i_y+2,:) ) * inv_840 * inv_dy
        else if(i_y == i_y_max - 1) then
            nabla_V_1_y(:,i_y,:) = (3*V_1(:,i_y-4,:) - 32*V_1(:,i_y-3,:) + 168*V_1(:,i_y-2,:) &
            - 672*V_1(:,i_y-1,:) + 672*V_1(:,i_y+1,:) ) * inv_840 * inv_dy
        else if(i_y == i_y_max) then
            nabla_V_1_y(:,i_y,:) = (3*V_1(:,i_y-4,:) - 32*V_1(:,i_y-3,:) + 168*V_1(:,i_y-2,:) &
            - 672*V_1(:,i_y-1,:) ) * inv_840 * inv_dy
        else
        nabla_V_1_y(:,i_y,:) = (3*V_1(:,i_y-4,:) - 32*V_1(:,i_y-3,:) + 168*V_1(:,i_y-2,:) &
        - 672*V_1(:,i_y-1,:) + 672*V_1(:,i_y+1,:) - 168*V_1(:,i_y+2,:) + 32*V_1(:,i_y+3,:) - 3*V_1(:,i_y+4,:)) * inv_840 * inv_dy
        end if
    end do
    !$omp end do
    !$omp do schedule(static)
    do i_z = i_z_min, i_z_max
        if (i_z == i_z_min  )then
            nabla_V_1_z(:,:,i_z) = (-3*V_1(:,:,i_z+4) + 32*V_1(:,:,i_z+3) - 168*V_1(:,:,i_z+2) &
            + 672*V_1(:,:,i_z+1) ) * inv_840 * inv_dz
        else if (i_z == i_z_min +  1) then
            nabla_V_1_z(:,:,i_z) = (-3*V_1(:,:,i_z+4) + 32*V_1(:,:,i_z+3) - 168*V_1(:,:,i_z+2) &
            + 672*V_1(:,:,i_z +1) -672*V_1(:,:,i_z -1)) * inv_840 * inv_dz
        else if (i_z == i_z_min + 2) then
            nabla_V_1_z(:,:,i_z) = (-3*V_1(:,:,i_z+4) + 32*V_1(:,:,i_z+3) - 168*V_1(:,:,i_z+2) &
            + 672*V_1(:,:,i_z+1) - 672*V_1(:,:,i_z-1) + 168*V_1(:,:,i_z-2)) * inv_840 * inv_dz
        else if (i_z == i_z_min + 3) then
            nabla_V_1_z(:,:,i_z) = (-3*V_1(:,:,i_z+4) + 32*V_1(:,:,i_z+3) - 168*V_1(:,:,i_z+2) &
            + 672*V_1(:,:,i_z+1) - 672*V_1(:,:,i_z-1) + 168*V_1(:,:,i_z-2) - 32*V_1(:,:,i_z-3)) * inv_840 * inv_dz
        else if (i_z == i_z_max - 3) then
            nabla_V_1_z(:,:,i_z) = (3*V_1(:,:,i_z-4) - 32*V_1(:,:,i_z-3) + 168*V_1(:,:,i_z-2) &
            - 672*V_1(:,:,i_z-1) + 672*V_1(:,:,i_z+1) - 168*V_1(:,:,i_z+2) + 32*V_1(:,:,i_z+3)) * inv_840 * inv_dz
        else if (i_z == i_z_max - 2) then
            nabla_V_1_z(:,:,i_z) = (3*V_1(:,:,i_z-4) - 32*V_1(:,:,i_z-3) + 168*V_1(:,:,i_z-2) &
            - 672*V_1(:,:,i_z-1) + 672*V_1(:,:,i_z+1) - 168*V_1(:,:,i_z+2) ) * inv_840 * inv_dz
        else if(i_z == i_z_max - 1) then
            nabla_V_1_z(:,:,i_z) = (3*V_1(:,:,i_z-4) - 32*V_1(:,:,i_z-3) + 168*V_1(:,:,i_z-2) &
            - 672*V_1(:,:,i_z-1) + 672*V_1(:,:,i_z+1) ) * inv_840 * inv_dz
        else if(i_z == i_z_max) then
            nabla_V_1_z(:,:,i_z) = (3*V_1(:,:,i_z-4) - 32*V_1(:,:,i_z-3) + 168*V_1(:,:,i_z-2) &
            - 672*V_1(:,:,i_z-1) ) * inv_840 * inv_dz
        else
        nabla_V_1_z(:,:,i_z) = (3*V_1(:,:,i_z-4) - 32*V_1(:,:,i_z-3) + 168*V_1(:,:,i_z-2) &
        - 672*V_1(:,:,i_z-1) + 672*V_1(:,:,i_z+1) - 168*V_1(:,:,i_z+2) + 32*V_1(:,:,i_z+3) - 3*V_1(:,:,i_z+4)) * inv_840 * inv_dz
        end if
    end do
    !$omp end do
    !$omp end parallel
    print *, "Nabla on V_1 computed."


end subroutine compute_nabla_on_V_1

subroutine prepare_initial_wavefunction
    !$ use omp_lib
    use parameter
    implicit none
    integer :: j_x, j_y, j_z, j_num,j_num_max! index for wavefunction number
    double precision :: norm_factor_local
    double precision, parameter :: sigma2 = 25.0d0 ! width of the Gaussian wave packet
    double precision :: x2,y2,z2,amp
    norm_factor = 0.0d0
    if (identifier == "proton") then
        j_num_max = Z_num + N_extra
    else if (identifier == "neutron") then
        j_num_max = N_num + N_extra
    end if

    do j_num = 1, j_num_max
    !$omp parallel do collapse(3) default(none) schedule(static) &
    !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
    !$omp        psi_up,psi_dn,j_num) &
    !$omp private(j_x,j_y,j_z,x2,y2,z2,amp) reduction(+:norm_factor_local)
        do j_z = i_z_min, i_z_max
            do j_y = i_y_min, i_y_max
                do j_x = i_x_min, i_x_max
                    x2 = (dble(j_x)*dx)**2
                    y2 = (dble(j_y)*dy)**2
                    z2 = (dble(j_z)*dz)**2
                    amp = exp( - ( x2 + y2 + z2 ) / sigma2 )
                    psi_up(j_x,j_y,j_z,j_num) = cmplx(amp, 0.0d0, kind(0d0) )
                    psi_dn(j_x,j_y,j_z,j_num) = cmplx(amp, 0.0d0, kind(0d0) )
                    norm_factor_local = norm_factor_local + amp**2 * dV *2.0 ! factor 2 for spin up and down
                end do
            end do
        end do
    !$omp end parallel do
        norm_factor(j_num) = sqrt(norm_factor_local)
    !$omp parallel do collapse(3) default(none) schedule(static) &
    !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
    !$omp        psi_up,psi_dn,norm_factor,j_num) &
    !$omp private(j_x,j_y,j_z)
        do j_z = i_z_min, i_z_max
            do j_y = i_y_min, i_y_max
                do j_x = i_x_min, i_x_max
                    psi_up(j_x,j_y,j_z,j_num) = psi_up(j_x,j_y,j_z,j_num) / norm_factor(j_num)
                    psi_dn(j_x,j_y,j_z,j_num) = psi_dn(j_x,j_y,j_z,j_num) / norm_factor(j_num)
                end do
                
                
            end do
        end do
    !$omp end parallel do
    end do

    print *, "Initial wavefunction prepared and normalized."
end subroutine prepare_initial_wavefunction

subroutine compute_nabla_on_psi
    !$ use omp_lib
    use parameter
    implicit none

    integer :: i_x, i_y, i_z,ia, j_num

    !$omp parallel do default(none) collapse(4) schedule(static) &
    !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
    !$omp        psi_up, psi_dn,&
    !$omp        dpsi_dx_up, dpsi_dy_up, dpsi_dz_up, &
    !$omp        dpsi_dx_dn, dpsi_dy_dn, dpsi_dz_dn,active,n_active) &
    !$omp private(i_x,i_y,i_z,ia,j_num)
    do ia = 1, n_active
        do i_z = i_z_min, i_z_max
            do i_y = i_y_min, i_y_max
                do i_x = i_x_min, i_x_max
                    j_num = active(ia)
                    if (i_x == i_x_min) then
                        dpsi_dx_up(i_x,i_y,i_z,j_num) = ( psi_up(i_x+1,i_y,i_z,j_num) ) *inv_dx*0.5d0
                        dpsi_dx_dn(i_x,i_y,i_z,j_num) = ( psi_dn(i_x+1,i_y,i_z,j_num) ) *inv_dx*0.5d0
                    else if (i_x == i_x_max) then
                        dpsi_dx_up(i_x,i_y,i_z,j_num) = (  - psi_up(i_x-1,i_y,i_z,j_num) ) *inv_dx*0.5d0
                        dpsi_dx_dn(i_x,i_y,i_z,j_num) = (  - psi_dn(i_x-1,i_y,i_z,j_num) ) *inv_dx*0.5d0
                    else if (i_x >= i_x_min+4 .and. i_x <= i_x_max-4) then
                        dpsi_dx_up(i_x,i_y,i_z,j_num) = ( &
                            3.0d0*psi_up(i_x-4,i_y,i_z,j_num) &
                            -32.0d0*psi_up(i_x-3,i_y,i_z,j_num) &
                        +168.0d0*psi_up(i_x-2,i_y,i_z,j_num) &
                        -672.0d0*psi_up(i_x-1,i_y,i_z,j_num) &
                        +672.0d0*psi_up(i_x+1,i_y,i_z,j_num) &
                        -168.0d0*psi_up(i_x+2,i_y,i_z,j_num) &
                            +32.0d0*psi_up(i_x+3,i_y,i_z,j_num) &
                            -3.0d0*psi_up(i_x+4,i_y,i_z,j_num) ) * (inv_dx/840.0d0)

                        dpsi_dx_dn(i_x,i_y,i_z,j_num) = ( &
                            3.0d0*psi_dn(i_x-4,i_y,i_z,j_num) &
                            -32.0d0*psi_dn(i_x-3,i_y,i_z,j_num) &
                        +168.0d0*psi_dn(i_x-2,i_y,i_z,j_num) &
                        -672.0d0*psi_dn(i_x-1,i_y,i_z,j_num) &
                        +672.0d0*psi_dn(i_x+1,i_y,i_z,j_num) &
                        -168.0d0*psi_dn(i_x+2,i_y,i_z,j_num) &
                            +32.0d0*psi_dn(i_x+3,i_y,i_z,j_num) &
                            -3.0d0*psi_dn(i_x+4,i_y,i_z,j_num) ) * (inv_dx/840.0d0)
                    else
                        dpsi_dx_up(i_x,i_y,i_z,j_num) = ( psi_up(i_x+1,i_y,i_z,j_num) - &
                                                            psi_up(i_x-1,i_y,i_z,j_num) ) * inv_dx*0.5d0
                        dpsi_dx_dn(i_x,i_y,i_z,j_num) = ( psi_dn(i_x+1,i_y,i_z,j_num) - &
                                                            psi_dn(i_x-1,i_y,i_z,j_num) ) * inv_dx*0.5d0
                    end if

                    if (i_y == i_y_min) then
                        dpsi_dy_up(i_x,i_y,i_z,j_num) = ( psi_up(i_x,i_y+1,i_z,j_num) ) * inv_dy * 0.5d0
                        dpsi_dy_dn(i_x,i_y,i_z,j_num) = ( psi_dn(i_x,i_y+1,i_z,j_num) ) * inv_dy * 0.5d0
                    else if (i_y == i_y_max) then
                        dpsi_dy_up(i_x,i_y,i_z,j_num) = (  - psi_up(i_x,i_y-1,i_z,j_num) ) * inv_dy * 0.5d0
                        dpsi_dy_dn(i_x,i_y,i_z,j_num) = (  - psi_dn(i_x,i_y-1,i_z,j_num) ) * inv_dy * 0.5d0
                    else if (i_y >= i_y_min+4 .and. i_y <= i_y_max-4) then
                        dpsi_dy_up(i_x,i_y,i_z,j_num) = ( &
                            3.0d0*psi_up(i_x,i_y-4,i_z,j_num) &
                            -32.0d0*psi_up(i_x,i_y-3,i_z,j_num) &
                        +168.0d0*psi_up(i_x,i_y-2,i_z,j_num) &
                        -672.0d0*psi_up(i_x,i_y-1,i_z,j_num) &
                        +672.0d0*psi_up(i_x,i_y+1,i_z,j_num) &
                        -168.0d0*psi_up(i_x,i_y+2,i_z,j_num) &
                            +32.0d0*psi_up(i_x,i_y+3,i_z,j_num) &
                            -3.0d0*psi_up(i_x,i_y+4,i_z,j_num) ) * (inv_dy/840.0d0)

                        dpsi_dy_dn(i_x,i_y,i_z,j_num) = ( &
                            3.0d0*psi_dn(i_x,i_y-4,i_z,j_num) &
                            -32.0d0*psi_dn(i_x,i_y-3,i_z,j_num) &
                        +168.0d0*psi_dn(i_x,i_y-2,i_z,j_num) &
                        -672.0d0*psi_dn(i_x,i_y-1,i_z,j_num) &
                        +672.0d0*psi_dn(i_x,i_y+1,i_z,j_num) &
                        -168.0d0*psi_dn(i_x,i_y+2,i_z,j_num) &
                            +32.0d0*psi_dn(i_x,i_y+3,i_z,j_num) &
                            -3.0d0*psi_dn(i_x,i_y+4,i_z,j_num) ) * (inv_dy/840.0d0)
                    else
                        dpsi_dy_up(i_x,i_y,i_z,j_num) = ( psi_up(i_x,i_y+1,i_z,j_num) - &
                                                            psi_up(i_x,i_y-1,i_z,j_num) ) * inv_dy * 0.5d0
                        dpsi_dy_dn(i_x,i_y,i_z,j_num) = ( psi_dn(i_x,i_y+1,i_z,j_num) - &
                                                            psi_dn(i_x,i_y-1,i_z,j_num) ) * inv_dy * 0.5d0
                    end if

                    if (i_z == i_z_min) then
                        dpsi_dz_up(i_x,i_y,i_z,j_num) = ( psi_up(i_x,i_y,i_z+1,j_num) ) * inv_dz * 0.5d0
                        dpsi_dz_dn(i_x,i_y,i_z,j_num) = ( psi_dn(i_x,i_y,i_z+1,j_num) ) * inv_dz * 0.5d0
                    else if (i_z == i_z_max) then
                        dpsi_dz_up(i_x,i_y,i_z,j_num) = ( - psi_up(i_x,i_y,i_z-1,j_num) ) * inv_dz * 0.5d0
                        dpsi_dz_dn(i_x,i_y,i_z,j_num) = ( - psi_dn(i_x,i_y,i_z-1,j_num) ) * inv_dz * 0.5d0
                    else if (i_z >= i_z_min+4 .and. i_z <= i_z_max-4) then
                        dpsi_dz_up(i_x,i_y,i_z,j_num) = ( &
                            3.0d0*psi_up(i_x,i_y,i_z-4,j_num) &
                            -32.0d0*psi_up(i_x,i_y,i_z-3,j_num) &
                        +168.0d0*psi_up(i_x,i_y,i_z-2,j_num) &
                        -672.0d0*psi_up(i_x,i_y,i_z-1,j_num) &
                        +672.0d0*psi_up(i_x,i_y,i_z+1,j_num) &
                        -168.0d0*psi_up(i_x,i_y,i_z+2,j_num) &
                            +32.0d0*psi_up(i_x,i_y,i_z+3,j_num) &
                            -3.0d0*psi_up(i_x,i_y,i_z+4,j_num) ) * (inv_dz/840.0d0)
                        dpsi_dz_dn(i_x,i_y,i_z,j_num) = ( &
                            3.0d0*psi_dn(i_x,i_y,i_z-4,j_num) &
                            -32.0d0*psi_dn(i_x,i_y,i_z-3,j_num) &
                        +168.0d0*psi_dn(i_x,i_y,i_z-2,j_num) &
                        -672.0d0*psi_dn(i_x,i_y,i_z-1,j_num) &
                        +672.0d0*psi_dn(i_x,i_y,i_z+1,j_num) &
                        -168.0d0*psi_dn(i_x,i_y,i_z+2,j_num) &
                            +32.0d0*psi_dn(i_x,i_y,i_z+3,j_num) &
                            -3.0d0*psi_dn(i_x,i_y,i_z+4,j_num) ) * (inv_dz/840.0d0)
                    else
                        dpsi_dz_up(i_x,i_y,i_z,j_num) = ( psi_up(i_x,i_y,i_z+1,j_num) - &
                                                            psi_up(i_x,i_y,i_z-1,j_num) ) * inv_dz * 0.5d0
                        dpsi_dz_dn(i_x,i_y,i_z,j_num) = ( psi_dn(i_x,i_y,i_z+1,j_num) - &
                                                            psi_dn(i_x,i_y,i_z-1,j_num) ) * inv_dz * 0.5d0
                    end if
                end do
            end do
        end do
    end do
    !$omp end parallel do

end subroutine compute_nabla_on_psi


subroutine compute_V_so
    !$ use omp_lib
    use parameter
    implicit none
    integer :: i_x, i_y, i_z, j_num,j_num_max,ia
    double precision :: pref

    call compute_nabla_on_psi
    if (potential_type == "Yukawa") then
        if (identifier == "proton") then
            pref = - (lambda_p * hbar_c**2 /(4.0d0 * m_nuc**2 ))*V_p
        else if (identifier == "neutron") then
            pref = - (lambda_n * hbar_c**2 /(4.0d0 * m_nuc**2 ))*V_n
        else
            print *, "Error: identifier must be 'proton' or 'neutron'."
            stop
        end if
    else if (potential_type == "Woods-Saxon") then
        pref = - 0.44*V_0
    else
        print *, "Error: potential_type must be 'Yukawa' or 'Woods-Saxon'."
        stop
    end if

    !$omp parallel default(none) &
    !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
    !$omp        n_active, V_so_up,V_so_dn, &
    !$omp        nabla_V_1_x,nabla_V_1_y,nabla_V_1_z, &
    !$omp        dpsi_dx_up,dpsi_dy_up,dpsi_dz_up, &
    !$omp        dpsi_dx_dn,dpsi_dy_dn,dpsi_dz_dn, pref, active) &
    !$omp private(i_x,i_y,i_z,ia,j_num)
    !$omp  do collapse(4) schedule(static)
    do ia = 1, n_active
        do i_z = i_z_min, i_z_max
            do i_y = i_y_min, i_y_max
                do i_x = i_x_min, i_x_max
                    j_num = active(ia)
                    ! ! Compute gradients of V_1 using central difference
                    V_so_up(i_x,i_y,i_z,j_num) = -pref*I_cmplx*(&
                    ( ( nabla_V_1_y(i_x,i_y,i_z) * dpsi_dz_dn(i_x,i_y,i_z,j_num) ) &
                    -(nabla_V_1_z(i_x,i_y,i_z) * dpsi_dy_dn(i_x,i_y,i_z,j_num) ) ) &
                    +( ( nabla_V_1_z(i_x,i_y,i_z) * dpsi_dx_dn(i_x,i_y,i_z,j_num)  &
                    - nabla_V_1_x(i_x,i_y,i_z) * dpsi_dz_dn(i_x,i_y,i_z,j_num) ) )*(-I_cmplx)&
                    +( ( nabla_V_1_x(i_x,i_y,i_z) * dpsi_dy_up(i_x,i_y,i_z,j_num) ) &
                    - ( nabla_V_1_y(i_x,i_y,i_z) * dpsi_dx_up(i_x,i_y,i_z,j_num) ) ))

                    V_so_dn(i_x,i_y,i_z,j_num) =-pref*I_cmplx*( &
                        ( ( nabla_V_1_y(i_x,i_y,i_z) * dpsi_dz_up(i_x,i_y,i_z,j_num) ) &
                    - ( nabla_V_1_z(i_x,i_y,i_z) * dpsi_dy_up(i_x,i_y,i_z,j_num) ) ) &
                    + ( ( nabla_V_1_z(i_x,i_y,i_z) * dpsi_dx_up(i_x,i_y,i_z,j_num) ) &
                    - ( nabla_V_1_x(i_x,i_y,i_z) * dpsi_dz_up(i_x,i_y,i_z,j_num) ) )*(I_cmplx) &
                    - ( ( nabla_V_1_x(i_x,i_y,i_z) * dpsi_dy_dn(i_x,i_y,i_z,j_num) ) &
                    - ( nabla_V_1_y(i_x,i_y,i_z) * dpsi_dx_dn(i_x,i_y,i_z,j_num) ) ) )
                end do
                
            end do
        end do

    end do
    !$omp end do
    !$omp end parallel

    
end subroutine compute_V_so

subroutine compute_Laplacian
    !$ use omp_lib
    use parameter
    implicit none
    integer :: i_x, i_y, i_z, j_num,j_num_max, ia
    double precision :: inv_5040
    complex(kind(0d0)):: d2psi_up_dx2, d2psi_up_dy2, d2psi_up_dz2
    complex(kind(0d0)):: d2psi_dn_dx2, d2psi_dn_dy2, d2psi_dn_dz2

    ! Compute Laplacian using central difference
    inv_5040 = 1.0d0 / 5040.0d0

    !$omp parallel default(none) &
    !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
    !$omp        n_active, psi_up, psi_dn, &  
    !$omp        laplacian_psi_up, laplacian_psi_dn, active,inv_5040) &
    !$omp private(i_x,i_y,i_z,ia,j_num, &
    !$omp        d2psi_up_dx2, d2psi_up_dy2, d2psi_up_dz2, &
    !$omp        d2psi_dn_dx2, d2psi_dn_dy2, d2psi_dn_dz2)
    !$omp do collapse(4) schedule(static)
    do ia = 1, n_active
        do i_z = i_z_min, i_z_max
            do i_y = i_y_min, i_y_max
                do i_x = i_x_min, i_x_max
                    j_num = active(ia)
                        if (i_x == i_x_min)then
                            d2psi_up_dx2 = ( psi_up(i_x+1,i_y,i_z,j_num) &
                            - 2.0d0*psi_up(i_x,i_y,i_z,j_num) ) 
                            d2psi_dn_dx2 = ( psi_dn(i_x+1,i_y,i_z,j_num) &
                            - 2.0d0*psi_dn(i_x,i_y,i_z,j_num) ) 
                        else if (i_x == i_x_max) then
                            d2psi_up_dx2 = ( - 2.0d0*psi_up(i_x,i_y,i_z,j_num) &
                            + psi_up(i_x-1,i_y,i_z,j_num) ) 
                            d2psi_dn_dx2 = ( - 2.0d0*psi_dn(i_x,i_y,i_z,j_num) &
                            + psi_dn(i_x-1,i_y,i_z,j_num) ) 
                        else if (i_x >= i_x_min+4 .and. i_x <= i_x_max-4) then
    ! 9-point (8th-order) second derivative
                            d2psi_up_dx2 = ( &
                                -  9.0d0*psi_up(i_x-4,i_y,i_z,j_num) &
                                +128.0d0*psi_up(i_x-3,i_y,i_z,j_num) &
                                -1008.0d0*psi_up(i_x-2,i_y,i_z,j_num) &
                                +8064.0d0*psi_up(i_x-1,i_y,i_z,j_num) &
                                -14350.0d0*psi_up(i_x  ,i_y,i_z,j_num) &
                                +8064.0d0*psi_up(i_x+1,i_y,i_z,j_num) &
                                -1008.0d0*psi_up(i_x+2,i_y,i_z,j_num) &
                                +128.0d0*psi_up(i_x+3,i_y,i_z,j_num) &
                                -  9.0d0*psi_up(i_x+4,i_y,i_z,j_num) ) *inv_5040

                            d2psi_dn_dx2 = ( &
                                -  9.0d0*psi_dn(i_x-4,i_y,i_z,j_num) &
                                +128.0d0*psi_dn(i_x-3,i_y,i_z,j_num) &
                                -1008.0d0*psi_dn(i_x-2,i_y,i_z,j_num) &
                                +8064.0d0*psi_dn(i_x-1,i_y,i_z,j_num) &
                                -14350.0d0*psi_dn(i_x  ,i_y,i_z,j_num) &
                                +8064.0d0*psi_dn(i_x+1,i_y,i_z,j_num) &
                                -1008.0d0*psi_dn(i_x+2,i_y,i_z,j_num) &
                                +128.0d0*psi_dn(i_x+3,i_y,i_z,j_num) &
                                -  9.0d0*psi_dn(i_x+4,i_y,i_z,j_num) ) *inv_5040
                        else
                            d2psi_up_dx2 = ( psi_up(i_x+1,i_y,i_z,j_num) &
                            - 2.0d0*psi_up(i_x,i_y,i_z,j_num) + &
                            psi_up(i_x-1,i_y,i_z,j_num) ) 
                            d2psi_dn_dx2 = ( psi_dn(i_x+1,i_y,i_z,j_num) &
                            - 2.0d0*psi_dn(i_x,i_y,i_z,j_num) + &
                            psi_dn(i_x-1,i_y,i_z,j_num) ) 
                        end if 
                        d2psi_dn_dx2 = d2psi_dn_dx2 * inv_dx**2
                        d2psi_up_dx2 = d2psi_up_dx2 * inv_dx**2

                    
                        if (i_y == i_y_min)then
                            d2psi_up_dy2 = ( psi_up(i_x,i_y+1,i_z,j_num) &
                            - 2.0d0*psi_up(i_x,i_y,i_z,j_num) ) 
                            d2psi_dn_dy2 = ( psi_dn(i_x,i_y+1,i_z,j_num) &
                            - 2.0d0*psi_dn(i_x,i_y,i_z,j_num) ) 
                        else if (i_y == i_y_max) then
                            d2psi_up_dy2 = ( - 2.0d0*psi_up(i_x,i_y,i_z,j_num) &
                            + psi_up(i_x,i_y-1,i_z,j_num) ) 
                            d2psi_dn_dy2 = ( - 2.0d0*psi_dn(i_x,i_y,i_z,j_num) &
                            + psi_dn(i_x,i_y-1,i_z,j_num) ) 
                        else if (i_y >= i_y_min+4 .and. i_y <= i_y_max-4) then
    ! 9-point (8th-order) second derivative
                            d2psi_up_dy2 = ( &
                                -  9.0d0*psi_up(i_x,i_y-4,i_z,j_num) &
                                +128.0d0*psi_up(i_x,i_y-3,i_z,j_num) &
                                -1008.0d0*psi_up(i_x,i_y-2,i_z,j_num) &
                                +8064.0d0*psi_up(i_x,i_y-1,i_z,j_num) &
                                -14350.0d0*psi_up(i_x,i_y  ,i_z,j_num) &
                                +8064.0d0*psi_up(i_x,i_y+1,i_z,j_num) &
                                -1008.0d0*psi_up(i_x,i_y+2,i_z,j_num) &
                                +128.0d0*psi_up(i_x,i_y+3,i_z,j_num) &
                                -  9.0d0*psi_up(i_x,i_y+4,i_z,j_num) ) *inv_5040

                            d2psi_dn_dy2 = ( &
                                -  9.0d0*psi_dn(i_x,i_y-4,i_z,j_num) &
                                +128.0d0*psi_dn(i_x,i_y-3,i_z,j_num) &
                                -1008.0d0*psi_dn(i_x,i_y-2,i_z,j_num) &
                                +8064.0d0*psi_dn(i_x,i_y-1,i_z,j_num) &
                                -14350.0d0*psi_dn(i_x,i_y  ,i_z,j_num) &
                                +8064.0d0*psi_dn(i_x,i_y+1,i_z,j_num) &
                                -1008.0d0*psi_dn(i_x,i_y+2,i_z,j_num) &
                                +128.0d0*psi_dn(i_x,i_y+3,i_z,j_num) &
                                -  9.0d0*psi_dn(i_x,i_y+4,i_z,j_num) ) *inv_5040
                        else
                            d2psi_up_dy2 = ( psi_up(i_x,i_y+1,i_z,j_num) &
                            - 2.0d0*psi_up(i_x,i_y,i_z,j_num) + &
                            psi_up(i_x,i_y-1,i_z,j_num) ) 
                            d2psi_dn_dy2 = ( psi_dn(i_x,i_y+1,i_z,j_num) &
                            - 2.0d0*psi_dn(i_x,i_y,i_z,j_num) + &
                            psi_dn(i_x,i_y-1,i_z,j_num) ) 
                        end if
                        d2psi_dn_dy2 = d2psi_dn_dy2 * inv_dy**2
                        d2psi_up_dy2 = d2psi_up_dy2 * inv_dy**2

                        if (i_z == i_z_min)then
                            d2psi_up_dz2 = ( psi_up(i_x,i_y,i_z+1,j_num) &
                            - 2.0d0*psi_up(i_x,i_y,i_z,j_num) ) 
                            d2psi_dn_dz2 = ( psi_dn(i_x,i_y,i_z+1,j_num) &
                            - 2.0d0*psi_dn(i_x,i_y,i_z,j_num) ) 
                        else if (i_z == i_z_max) then
                            d2psi_up_dz2 = ( - 2.0d0*psi_up(i_x,i_y,i_z,j_num) &
                            + psi_up(i_x,i_y,i_z-1,j_num) ) 
                            d2psi_dn_dz2 = ( - 2.0d0*psi_dn(i_x,i_y,i_z,j_num) &
                            + psi_dn(i_x,i_y,i_z-1,j_num) ) 
                        else if (i_z >= i_z_min+4 .and. i_z <= i_z_max-4) then
    ! 9-point (8th-order) second derivative
                            d2psi_up_dz2 = ( &
                                -  9.0d0*psi_up(i_x,i_y,i_z-4,j_num) &
                                +128.0d0*psi_up(i_x,i_y,i_z-3,j_num) &
                                -1008.0d0*psi_up(i_x,i_y,i_z-2,j_num) &
                                +8064.0d0*psi_up(i_x,i_y,i_z-1,j_num) &
                                -14350.0d0*psi_up(i_x,i_y,i_z  ,j_num) &
                                +8064.0d0*psi_up(i_x,i_y,i_z+1,j_num) &
                                -1008.0d0*psi_up(i_x,i_y,i_z+2,j_num) &
                                +128.0d0*psi_up(i_x,i_y,i_z+3,j_num) &
                                -  9.0d0*psi_up(i_x,i_y,i_z+4,j_num) ) *inv_5040

                            d2psi_dn_dz2 = ( &
                                -  9.0d0*psi_dn(i_x,i_y,i_z-4,j_num) &
                                +128.0d0*psi_dn(i_x,i_y,i_z-3,j_num) &
                                -1008.0d0*psi_dn(i_x,i_y,i_z-2,j_num) &
                                +8064.0d0*psi_dn(i_x,i_y,i_z-1,j_num) &
                                -14350.0d0*psi_dn(i_x,i_y,i_z  ,j_num) &
                                +8064.0d0*psi_dn(i_x,i_y,i_z+1,j_num) &
                                -1008.0d0*psi_dn(i_x,i_y,i_z+2,j_num) &
                                +128.0d0*psi_dn(i_x,i_y,i_z+3,j_num) &
                                -  9.0d0*psi_dn(i_x,i_y,i_z+4,j_num) ) *inv_5040
                        else
                            d2psi_up_dz2 = ( psi_up(i_x,i_y,i_z+1,j_num) &
                            - 2.0d0*psi_up(i_x,i_y,i_z,j_num) + &
                            psi_up(i_x,i_y,i_z-1,j_num) ) 
                            d2psi_dn_dz2 = ( psi_dn(i_x,i_y,i_z+1,j_num) &
                            - 2.0d0*psi_dn(i_x,i_y,i_z,j_num) + &
                            psi_dn(i_x,i_y,i_z-1,j_num) ) 
                        end if
                        d2psi_dn_dz2 = d2psi_dn_dz2 * inv_dz**2
                        d2psi_up_dz2 = d2psi_up_dz2 * inv_dz**2

                        Laplacian_psi_up(i_x,i_y,i_z,j_num) = d2psi_up_dx2 + d2psi_up_dy2 &
                        + d2psi_up_dz2

                        Laplacian_psi_dn(i_x,i_y,i_z,j_num) = d2psi_dn_dx2 + d2psi_dn_dy2 &
                        + d2psi_dn_dz2

                end do
                
                
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

end subroutine compute_Laplacian

subroutine  Imaginary_Time_Evolution
    !$ use omp_lib
    use parameter
    implicit none
    integer :: j_x, j_y, j_z,j_num,j_num_max, ia, old_n_active
    character(len=256) :: filename
    double precision, allocatable :: Etot_prev(:)
    double precision :: pref

    pref = - hbar_c**2/(2.0d0 * m_nuc)
    
    ! delta_tau = 0.01d0 ! imaginary time step
    ! max_iter = 100000
    print *, "Starting Imaginary Time Evolution for ", trim(identifier), "s."
    if (potential_type == "Yukawa") then
        print *, "Using Yukawa Potential."
            if (identifier == "proton") then
            j_num_max = Z_num + N_extra
            allocate( E_single_particle_proton(j_num_max) )
            if(allocated(V_particle)) then
                deallocate(V_particle)
            end if
            allocate( V_particle(i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max) )
            V_particle(:,:,:) = V_1(:,:,:)*V_p+ V_coulomb(:,:,:)
        else if (identifier == "neutron") then
            j_num_max = N_num + N_extra
            allocate( E_single_particle_neutron(j_num_max) )
            if(allocated(V_particle)) then
                deallocate(V_particle)
            end if
            allocate( V_particle(i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max) )
            V_particle(:,:,:) = V_1(:,:,:)*V_n
        else
            print *, "Error: identifier must be 'proton' or 'neutron'."
            stop
        end if
    else if (potential_type == "Woods-Saxon") then
        print *, "Using Woods-Saxon Potential."
        if(allocated(V_particle)) then
            deallocate(V_particle)
        end if
        allocate( V_particle(i_x_min:i_x_max, i_y_min:i_y_max, i_z_min:i_z_max) )
        V_particle(:,:,:) = V_1(:,:,:)
    else
        print *, "Error: potential_type must be 'Yukawa' or 'Woods-Saxon'."
        stop
    end if
    
    if (allocated(Etot_prev)) then
        deallocate(Etot_prev)
    end if
    if (allocated(delta_Etot)) then
        deallocate(delta_Etot)
    end if
    if (allocated(active)) then
        deallocate(active)
    end if
    allocate(Etot_prev(j_num_max),delta_Etot(j_num_max),&
    active(j_num_max))
    
    Etot_prev(:) = 0.0d0
    n_active = j_num_max
    active(:) = (/(ia, ia = 1, j_num_max)/)

    
    do iter = 1, max_iter
        call compute_Laplacian
        call compute_V_so
        ! Compute Hamiltonian acting on wavefunctions
        !$omp parallel do default(none) collapse(4) schedule(static)&
        !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
        !$omp        n_active, psi_up, psi_dn, &
        !$omp        Laplacian_psi_up, Laplacian_psi_dn, &
        !$omp        V_particle, V_so_up, V_so_dn, &
        !$omp        h_psi_up, h_psi_dn,pref, active) &
        !$omp private(j_x,j_y,j_z,ia,j_num)
        do ia = 1, n_active
            do j_z = i_z_min, i_z_max
                do j_y = i_y_min, i_y_max
                    do j_x = i_x_min, i_x_max
                        j_num = active(ia)
                        h_psi_up(j_x,j_y,j_z,j_num) =&
                         ( pref *Laplacian_psi_up(j_x,j_y,j_z,j_num) + &
                        V_particle(j_x,j_y,j_z)*psi_up(j_x,j_y,j_z,j_num) ) + &
                        V_so_up(j_x,j_y,j_z,j_num)

                        h_psi_dn(j_x,j_y,j_z,j_num) = &
                        ( pref *Laplacian_psi_dn(j_x,j_y,j_z,j_num) + &
                        V_particle(j_x,j_y,j_z)*psi_dn(j_x,j_y,j_z,j_num) ) + &
                        V_so_dn(j_x,j_y,j_z,j_num)
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! Compute total energy for each state
        do ia = 1, n_active
            j_num = active(ia)
            E_single_particle(j_num) = 0.0d0
        end do
        !$omp parallel do default(none) collapse(4) schedule(static) &
        !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
        !$omp        n_active, psi_up, psi_dn, h_psi_up, h_psi_dn, E_single_particle, active)&
        !$omp private(j_x,j_y,j_z,j_num,ia)
        do ia = 1, n_active
            do j_z = i_z_min, i_z_max
                do j_y = i_y_min, i_y_max
                    do j_x = i_x_min, i_x_max
                        j_num = active(ia)
                        !$omp atomic
                        E_single_particle(j_num) = E_single_particle(j_num) + dble(conjg(psi_up(j_x,j_y,j_z,j_num)) *&
                         h_psi_up(j_x,j_y,j_z,j_num)) * dV
                        ! norm(j_num) = norm(j_num) + real(psi_up(j_x,j_y,j_z,j_num))**2 * dV +&
                        !  aimag(psi_up(j_x,j_y,j_z,j_num))**2 * dV
                        ! h_psi_norm(j_num) = h_psi_norm(j_num) + real(h_psi_up(j_x,j_y,j_z,j_num))**2 * dV +&
                        !  aimag(h_psi_up(j_x,j_y,j_z,j_num))**2 * dV
                        !$omp atomic
                        E_single_particle(j_num) = E_single_particle(j_num) + dble(conjg(psi_dn(j_x,j_y,j_z,j_num)) *&
                         h_psi_dn(j_x,j_y,j_z,j_num)) * dV
                        ! norm(j_num) = norm(j_num) + real(psi_dn(j_x,j_y,j_z,j_num))**2 * dV +&
                        !  aimag(psi_dn(j_x,j_y,j_z,j_num))**2 * dV
                        ! h_psi_norm(j_num) = h_psi_norm(j_num) + real(h_psi_dn(j_x,j_y,j_z,j_num))**2 * dV +&
                        !  aimag(h_psi_dn(j_x,j_y,j_z,j_num))**2 * dV
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        if (mod(iter, 1000) == 0) then
            print *, "Iteration: ", iter
            do j_num = 1, j_num_max
                print *, "Energy for state ", j_num, ": ", E_single_particle(j_num), " MeV"
            end do
            print *, "Number of active states: ", n_active
        end if

        ! Check for convergence
        if (iter > 1) then
            ! Resi_norm(:) = 0.0d0
            delta_Etot(:) = abs(E_single_particle(:) - Etot_prev(:))

            ! ! Check if all states have converged
            ! !$omp parallel do default(none) collapse(4) schedule(static) &
            ! !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
            ! !$omp        n_active, psi_up, psi_dn, h_psi_up, h_psi_dn, &
            ! !$omp        E_single_particle, Resi_norm, active) &
            ! !$omp private(j_x,j_y,j_z,ia,j_num)
            ! do ia = 1, n_active
            !     do j_z = i_z_min, i_z_max
            !         do j_y = i_y_min, i_y_max
            !             do j_x = i_x_min, i_x_max
            !                 j_num = active(ia)
            !                 ! Compute residual norm for convergence check
            !                 !$omp atomic
            !                 Resi_norm(j_num) = Resi_norm(j_num) + &
            !                 abs( h_psi_up(j_x,j_y,j_z,j_num) - E_single_particle(j_num)*&
            !                     psi_up(j_x,j_y,j_z,j_num) )**2 *dV

            !                 !$omp atomic
            !                 Resi_norm(j_num) = Resi_norm(j_num) + &
            !                 abs( h_psi_dn(j_x,j_y,j_z,j_num) - E_single_particle(j_num)*&
            !                     psi_dn(j_x,j_y,j_z,j_num) )**2 *dV
                            
            !             end do
            !         end do
            !     end do
            ! end do
            ! !$omp end parallel do
            old_n_active = n_active
            n_active = 0
            ! Resi_norm(:) = sqrt(Resi_norm(:))
            ! norm(:) = sqrt(norm(:))
            ! h_psi_norm(:) = sqrt(h_psi_norm(:))
            do ia = 1, old_n_active
                j_num = active(ia)
                ! Resi_norm(j_num) = Resi_norm(j_num) / ( norm(j_num) *abs(E_single_particle(j_num)) +&
                !  h_psi_norm(j_num) )
                ! if (Resi_norm(j_num) > 1.0d-5 .and. delta_Etot(j_num) > 1.0d-12) then
                if ( delta_Etot(j_num) > 1.0d-10) then
                    n_active = n_active + 1
                    active(n_active) = j_num
                end if
            end do

            if ( n_active == 0 .or. iter == 50000) then
                print *, "Converged after ", iter, " iterations."
                ! Output final energies 
                write(filename, '(A,A,A,I0,A,I0,A)') 'Energy_Levels_', trim(identifier),'of_',Z_num,'-',N_num, '.dat'
                open(unit=30, file=trim(filename), status='unknown')
                do j_num = 1, j_num_max
                    print *, "Final energy for state ", j_num, ": ", E_single_particle(j_num), " MeV"
                    ! write(filename, '(A,I0,A,A,A,I0,A,I0,A)') 'Final_Wavefunctions_', j_num, '_', &
                    ! trim(identifier),'of_',Z_num,'-',N_num, '.dat'
                    ! write energy levels to one file
                    ! open(unit=30, file=trim(filename), status='unknown')
                    ! write energy levels at the top of the file
                    write(30,*) '# Final energy for state ', j_num, ': ', E_single_particle(j_num), ' MeV'
                    ! do j_z = i_z_min, i_z_max
                    !     do j_y = i_y_min, i_y_max
                    !         do j_x = i_x_min, i_x_max
                    !             write(30,*) j_x*dx, j_y*dy, j_z*dz,j_num, real(psi_up(j_x,j_y,j_z,j_num)), &
                    !             aimag(psi_up(j_x,j_y,j_z,j_num)), &
                    !             real(psi_dn(j_x,j_y,j_z,j_num)), aimag(psi_dn(j_x,j_y,j_z,j_num)),E_single_particle(j_num)
                    !         end do
                    !     end do
                    ! end do
                    !close(30)
                end do
                write(30,*) "Converged after ", iter, " iterations."
                close(30)
                if (identifier == "proton") then
                    E_single_particle_proton = E_single_particle
                else if (identifier == "neutron") then
                    E_single_particle_neutron = E_single_particle
                end if
                ! check orthogonality of final wavefunctions
                ! print *, "Checking orthogonality of final wavefunctions..."
                ! print *, "Overlap matrix:"
                ! call Check_Orthogonality
                deallocate(Etot_prev,delta_Etot,active)
                exit
            end if
            
        end if
        Etot_prev = E_single_particle
        
        ! Update wavefunctions
        !$omp parallel do default(none) collapse(4) schedule(static)&
        !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max,&
        !$omp        n_active, psi_up, psi_dn, h_psi_up, h_psi_dn, active)&
        !$omp private(j_x,j_y,j_z,j_num,ia)
        do ia = 1, n_active
            do j_z = i_z_min, i_z_max
                do j_y = i_y_min, i_y_max
                    do j_x = i_x_min, i_x_max
                        j_num = active(ia)
                        psi_up(j_x,j_y,j_z,j_num) = psi_up(j_x,j_y,j_z,j_num) - delta_tau * h_psi_up(j_x,j_y,j_z,j_num)
                        psi_dn(j_x,j_y,j_z,j_num) = psi_dn(j_x,j_y,j_z,j_num) - delta_tau * h_psi_dn(j_x,j_y,j_z,j_num)
                        
                    end do
                    
                end do
            end do
        end do
        !$omp end parallel do
        call Gram_Schmidt_Orthogonalization
        
    end do



end subroutine Imaginary_Time_Evolution

subroutine Check_Orthogonality
    use parameter
    implicit none
    integer :: i_x, i_y, i_z, j_num, k_num, ia, ja, j_num_max
    double precision :: overlap_2
    complex(kind(0d0)) :: overlap
    if (identifier == "proton") then
        j_num_max = Z_num + N_extra
    else if (identifier == "neutron") then
        j_num_max = N_num + N_extra
    else
        print *, "Error: identifier must be 'proton' or 'neutron'."
        stop
    end if
    do ia = 1, j_num_max
        do ja = 1, j_num_max
            if (ja < ia) then
                overlap = (0.0d0, 0.0d0)
                !$omp parallel default(none) &
                !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
                !$omp        psi_up,psi_dn,ia,ja,overlap)&
                !$omp private(i_x, i_y, i_z)
                !$omp do collapse(3) schedule(static) reduction(+:overlap)
                do i_z = i_z_min, i_z_max
                    do i_y = i_y_min, i_y_max
                        do i_x = i_x_min, i_x_max
                            overlap = overlap + ( conjg(psi_up(i_x,i_y,i_z,ja))*psi_up(i_x,i_y,i_z,ia) + &
                                            conjg(psi_dn(i_x,i_y,i_z,ja))*psi_dn(i_x,i_y,i_z,ia) ) * dV
                        end do
                    end do
                end do
                !$omp end do
                !$omp end parallel
                overlap_2 = dble(overlap*conjg(overlap))
                if (overlap_2 > 1.0d-6) then
                    print *, "Warning: Significant overlap detected between state ", ia, " and state ", ja, ": ", overlap_2
                end if
            end if
        end do
    end do

end subroutine Check_Orthogonality

subroutine Gram_Schmidt_Orthogonalization 
    !$ use omp_lib
    use parameter
    implicit none

    integer :: i_x, i_y, i_z, j_num, k_num,ia
    complex(kind(0d0)) :: proj,tmp
    double precision :: norm2, invnorm
    double precision :: proj_real, proj_imag

    do ia = 1, n_active
        j_num = active(ia)
        ! --- まず過去の状態に対して直交化 ---

        do k_num = 1, j_num - 1
            proj_imag = 0.0d0
            proj_real = 0.0d0
            !$omp parallel default(none) &
            !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max, &
            !$omp        psi_up,psi_dn,j_num,k_num,proj_real,proj_imag)&
            !$omp private(i_x, i_y, i_z,tmp,proj)
            !$omp do collapse(3) schedule(static) reduction(+:proj_real, proj_imag)
            do i_z = i_z_min, i_z_max
                do i_y = i_y_min, i_y_max
                    do i_x = i_x_min, i_x_max
                        tmp = ( conjg(psi_up(i_x,i_y,i_z,k_num))*psi_up(i_x,i_y,i_z,j_num) + &
                                        conjg(psi_dn(i_x,i_y,i_z,k_num))*psi_dn(i_x,i_y,i_z,j_num) ) * dV
                        proj_real = proj_real + dble( tmp)
                        proj_imag = proj_imag + dimag( tmp)
                    
                    end do
                end do
            end do
            !$omp end do
            !$omp do collapse(3) schedule(static)
            do i_z = i_z_min, i_z_max
                do i_y = i_y_min, i_y_max
                    do i_x = i_x_min, i_x_max
                        proj = cmplx(proj_real, proj_imag, kind(0d0))
                        psi_up(i_x,i_y,i_z,j_num) = psi_up(i_x,i_y,i_z,j_num) - proj*psi_up(i_x,i_y,i_z,k_num)
                        psi_dn(i_x,i_y,i_z,j_num) = psi_dn(i_x,i_y,i_z,j_num) - proj*psi_dn(i_x,i_y,i_z,k_num)
                    end do
                    
                end do
            end do
            !$omp end do
            !$omp end parallel
        end do


        ! --- 直交化後に必ず正規化 ---
        norm2 = 0.0d0
        !$omp parallel default(none) &
        !$omp shared(i_x_min,i_x_max,i_y_min,i_y_max,i_z_min,i_z_max,psi_up,psi_dn,j_num,invnorm,norm2) &
        !$omp private(i_x,i_y,i_z)
        !$omp do collapse(3) schedule(static) reduction(+:norm2)
        do i_z = i_z_min, i_z_max
            do i_y = i_y_min, i_y_max
                do i_x = i_x_min, i_x_max
                    norm2 = norm2 + ( abs(psi_up(i_x,i_y,i_z,j_num))**2 + abs(psi_dn(i_x,i_y,i_z,j_num))**2 ) * dV
                end do
            end do
        end do
        !$omp end do
        !$omp single
        if (norm2 <= 0.0d0) stop "GS failed: norm2 <= 0"

        invnorm = 1.0d0 / sqrt(norm2)
        !$omp end single

        !$omp do collapse(3) schedule(static)
        do i_z = i_z_min, i_z_max
            do i_y = i_y_min, i_y_max
                do i_x = i_x_min, i_x_max
                    psi_up(i_x,i_y,i_z,j_num) = psi_up(i_x,i_y,i_z,j_num) * invnorm
                    psi_dn(i_x,i_y,i_z,j_num) = psi_dn(i_x,i_y,i_z,j_num) * invnorm
                end do
                
            end do
        end do
        !$omp end do
        !$omp end parallel
    end do

end subroutine Gram_Schmidt_Orthogonalization
!!!!!!!!!!!!!!!FRLDM!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_FRLDM_only_for_spheroidal
    use parameter
    implicit none
    double precision :: b_1 
    double precision :: b_3
    double precision :: b_w
    ! double precision :: b_s
    
    ! ! double precision :: ecc !spheroidal eccentricity

    ! b_s = surface_area

    ! Calculate coefficients
    gap_p = r_mac*b_s/ (dble(Z_num))**(1.0d0/3.0d0)
    gap_n = r_mac*b_s/ (dble(N_num))**(1.0d0/3.0d0)
    gap_np = h_np/ b_s/(dble(A_num))**(2.0d0/3.0d0)

    ! calculate mass excess of Z hydrogens and N neutrons
    M_ZN = dble(Z_num)*M_H + dble(N_num)*M_n

    ! calculate volume energy
    I_n = dble(N_num - Z_num)/dble(A_num)
    VE = - a_v*(1.0d0 - k_v*I_n**2)*dble(A_num)

    ! calculate surface energy
    ! compute b_1 
    call compute_b_1_b_3_integration(b_1, b_3)
    SE = a_s*(1.0d0 - k_s*I_n**2)*b_1 &
    /(8*pi**2*r_0**2*a_Y**4)
    write(*,*) "SE=", SE

    ! x_0 = r_0*A_num**(1.0d0/3.0d0)/a_Y
    ! y_0 = r_0*A_num**(1.0d0/3.0d0)/a_den
    ! b_1 = 1.0d0 - 3.0d0/x_0**2 + (1.0d0 + x_0)*(2.0d0 + 3.0d0/x_0 + 3.0d0/x_0**2)*exp(-2.0d0*x_0)
    ! SE = a_s*(1.0d0 - k_s*I_n**2)*dble(A_num)**(2.0d0/3.0d0)*b_1


    ! calculate a0 energy
    !! compute b_w from deformation parameters
    ! then now, i can't understand parameter c_3, so I temporarily use b_w=1.0d0
    b_w = 1.0d0
    !A0E = a_0*b_w

    ! calculate coulomb energy
    c_1 = 3*e_2/(5.0d0*r_0)
    CE = c_1 * Z_num**2 / dble(A_num)**(1.0d0/3.0d0) * b_3*&
    15/32/pi**2 * A_num**(-5.0d0/3.0d0) /r_0**5
    write(*,*) "CE=", CE

    ! calculate coulomb exchange energy
    c_4 = 5.0d0/4.0d0*(3.0d0/(2.0d0*pi))**(2.0d0/3.0d0)*c_1
    CEE = - c_4 * Z_num**(4.0d0/3.0d0) / dble(A_num)**(1.0d0/3.0d0)

    ! calculate proton form factor correction to the coulomb energy
    k_f = (9*pi*Z_num/(4.0d0*A_num))**(1.0d0/3.0d0)/r_0
    f_kfr_p = - 1d0/8d0 * r_p**2 *e_2 / r_0**3 * &
    (145.0d0/48.0d0 - 327d0/2880d0 *(k_f*r_p)**2 + 1527d0/1209600d0 *(k_f*r_p)**4)

    PF2CE = f_kfr_p * Z_num**2 / A_num

    ! calculate charge asymmetry energy
    CAE = - c_a * ( N_num - Z_num)

    ! calculate Wigner energy
    if (mod(Z_num,2) ==1 .and. mod(N_num,2) == 1 .and. Z_num -N_num ==0) then
        WE = W_const *( abs(I_n)*b_w +1/A_num)
    else
        WE = W_const * abs(I_n)*b_w 
    end if

    ! calculate average pairing energy
    If (mod(Z_num,2) == 1 .and. mod(N_num,2) == 1) then
        APE = gap_n + gap_p - gap_np
    else if (mod(Z_num,2) == 1 .and. mod(N_num,2) == 0) then
        APE = gap_p
    else if (mod(Z_num,2) == 0 .and. mod(N_num,2) == 1) then
        APE = gap_n
    else
        APE = 0.0d0
    end if

    ! calculate electronic binding energy
    EoBE = - a_el * Z_num**2.39d0

    ! calculate FRLDM energy
    !FRLDM_energy = M_ZN + VE + SE + A0E + CE + CEE + PF2CE + CAE + WE + APE + EoBE
    FRLDM_energy =  VE + SE + CE + CEE + PF2CE + CAE + WE + APE + EoBE
    print *, "FRLDM energy=", FRLDM_energy

    deallocate(exp_a_Y)
    deallocate(exp_a_den)

end subroutine compute_FRLDM_only_for_spheroidal


! Compute b_1 coefficient for Coulomb energy calculation
subroutine compute_b_1_b_3_integration(b_1, b_3)
    !$ use omp_lib
    use parameter
    implicit none
    double precision, intent(out) :: b_1
    double precision, intent(out) :: b_3
    double precision :: b_1_sum
    double precision :: b_3_sum
    double precision :: b_1_1 ! first integration result
    double precision :: b_3_1 ! first integration result

    integer :: j_x, j_y, j_z, j_x_dot, j_y_dot, j_z_dot
    integer :: j_x2x, j_y2y, j_z2z
    double precision :: on6,a_Y_inverse, a_den_inverse ! temporary variables to optimize calculation

    integer :: vec_dot(3)
    on6 = 1.0d0 / 6.0d0
    a_Y_inverse = 1.0d0 / a_Y
    a_den_inverse = 1.0d0 / a_den

    b_1_sum = 0.0d0
    b_3_sum = 0.0d0
    !$omp parallel do default(none) collapse(3) schedule(guided) &
    !$omp shared(j_x_min,j_x_max,j_y_min,j_y_max,j_z_min,j_z_max, &
    !$omp        i_shape, r_relative, r_relative_inverse, exp_a_Y, exp_a_den, &
    !$omp        on6,  a_Y_inverse, a_den_inverse) &
    !$omp private(j_x,j_y,j_z,j_x_dot,j_y_dot,j_z_dot,j_x2x,j_y2y,j_z2z, &
    !$omp         b_1_1,b_3_1)&
    !$omp reduction(+:b_1_sum,b_3_sum)
    do j_z = j_z_min, j_z_max
        do j_y = j_y_min, j_y_max
            do j_x = j_x_min, j_x_max
                if (i_shape(j_x,j_y,j_z) == 0)cycle ! skip calculation for spherical shape
                b_1_1 = 0.0d0
                b_3_1 = 0.0d0
                do j_z_dot = j_z_min, j_z_max
                    do j_y_dot = j_y_min, j_y_max
                        do j_x_dot = j_x_min, j_x_max
                            if (i_shape(j_x_dot,j_y_dot,j_z_dot)==0) cycle ! skip calculation for spherical shape
                            j_x2x = j_x - j_x_dot
                            j_y2y = j_y - j_y_dot
                            j_z2z = j_z - j_z_dot
                            if (r_relative(j_x2x,j_y2y,j_z2z) > 0.0d0) then
                                b_1_1 = b_1_1 + (2-r_relative(j_x2x,j_y2y,j_z2z)*a_Y_inverse) * exp_a_Y(j_x2x,j_y2y,j_z2z) &
                                *r_relative_inverse(j_x2x,j_y2y,j_z2z) * a_Y 
                            else if (r_relative(j_x2x,j_y2y,j_z2z) == 0.0d0) then
                                ! b_1_1 = b_1_1 + 4*pi*a_Y*(dx**2 + dy**2 + dz**2) + &
                                ! 2*(-1 + r_relative * a_Y_inverse *0.5 - r_relative**2 * a_Y_inverse**2 * on6)* &
                                ! dV - exp_r_a_Y* dV
                                continue
                            end if                                
                            if (r_relative(j_x2x,j_y2y,j_z2z) > 0.001d0) then
                                b_3_1 = b_3_1 + ((1- exp_a_den(j_x2x,j_y2y,j_z2z))* &
                                r_relative_inverse(j_x2x,j_y2y,j_z2z) + &
                                exp_a_den(j_x2x,j_y2y,j_z2z)*0.5 * a_den_inverse )
                            else
                                b_3_1 = b_3_1 - (exp_a_den(j_x2x,j_y2y,j_z2z)*0.5 * a_den_inverse + &
                                a_den_inverse -&
                                0.5*a_den_inverse**2 * r_relative(j_x2x,j_y2y,j_z2z) + &
                                on6 * a_den_inverse**3 * r_relative(j_x2x,j_y2y,j_z2z)**2 )
                            end if
                        end do 
                    end do
                end do
                b_1_sum = b_1_sum + b_1_1
                b_3_sum = b_3_sum + b_3_1
            end do
        end do
    end do
    !$omp end parallel do
    b_1 = b_1_sum * dV**2
    b_3 = b_3_sum * dV**2

end subroutine compute_b_1_b_3_integration

! !!!!!!!!!!!!!!!!!!!!! Strutinsky calculation subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bubbleSort(vec)
    use parameter
    implicit none
    double precision, intent(inout) :: vec(int(n_levels/2))
    integer :: i,j,n
    double precision :: temp
    n = size(vec)
    do i=1,n-1
        do j=2,n-i+1
            if( vec(j-1) > vec(j) ) then
                temp = vec(j-1)
                vec(j-1) = vec(j)
                vec(j) = temp
            end if
        end do
    end do
end subroutine bubbleSort


subroutine compute_smoothed_fermi_energy
    !$ use omp_lib
    use parameter
    implicit none
    double precision :: e_lower, e_upper, e_mid, n_levels_mid,e_lowest! variables for Fermi energy search
    integer :: j 
    integer :: max_j
    double precision :: e_min, e_max, L, N_upper
    if(allocated(eps)) then
        deallocate(eps)
    end if
    if (identifier == "proton") then
        n_levels = Z_num + N_extra
        allocate(eps(int(n_levels/2)))
        do j = 1, int(n_levels/2)
            eps(j) = E_single_particle_proton(2*j-1)  ! take only one spin state
        end do
        n_tot = Z_num
    else if (identifier == "neutron") then
        n_levels = N_num + N_extra
        allocate(eps(int(n_levels/2)))
        do j = 1, int(n_levels/2)
            eps(j) = E_single_particle_neutron(2*j-1)  ! take only one spin state
        end do
        n_tot = N_num
    end if
    call bubbleSort(eps)
    ! print *, "Sorted single-particle energy levels for ", trim(identifier), ":"
    ! do j = 1, int(n_levels/2)
    !     print *, "Level ", j, ": ", eps(j), " MeV"
    ! end do
    ! Find Fermi energy by solving N = integral of smoothed level density up to Fermi energy
    L = 20.0d0   ! 10だと足りないことがあるのでまず20
    e_min = minval(eps)
    e_max = maxval(eps)

    e_lowest = -100.0d0  !十分に低い値に設定
    e_lower  = e_lowest
    e_upper  = e_max + L*gamma_shell

    ! do while (abs(e_upper - e_lower) > 1.0d-6)
    !     e_mid = 0.5d0 * (e_lower + e_upper)
    !     n_levels_mid = 0.0d0
    !     max_j =  int((e_mid - e_lowest)/delta_E + 1.0d0)
 
    !     do j = 1, max_j
    !         n_levels_mid = n_levels_mid + 2*smoothed_level_density(e_lowest + &
    !         dble(j-1)*delta_E) * delta_E
    !     end do

    !     if (abs(n_levels_mid - dble(n_tot)) < 1.0d-4) exit
    !     if (n_levels_mid < dble(n_tot)) then
    !         e_lower = e_mid
    !     else
    !         e_upper = e_mid
    !     end if
    ! end do
    E_Fermi_tilde = -2.0d0  
    print *, "Smoothed Fermi energy for ", trim(identifier), ": ", E_Fermi_tilde, " MeV"

    ! n_levels_mid = 0.0d0
    ! max_j = int((E_Fermi_tilde - e_lowest)/delta_E + 1.0d0)
    ! do j = 1, max_j
    !     n_levels_mid = n_levels_mid + 2*smoothed_level_density(e_lowest + &
    !     dble(j-1)*delta_E) * delta_E
    ! end do
    ! print *, "Number of levels up to smoothed Fermi energy for ", trim(identifier), ": ", n_levels_mid


    !!!!!!!!!!!!!!!!!!calculate sum of smoothed energy levels !!!!!!!!!!!!!!!!!!
    ! E_tilde_shell = 0.0d0

    ! max_j = int((E_Fermi_tilde - e_lowest)/delta_E + 1.0d0)
    !
    ! do j = 1, max_j
    !     E_tilde_shell = E_tilde_shell + 2*(e_lowest + dble(j-1)*delta_E) * &
    !     smoothed_level_density(e_lowest + dble(j-1)*delta_E) * delta_E
    ! end do
    
    print *, "Smoothed total energy levels for ", trim(identifier), ": ", E_tilde_shell, " MeV"
    open(unit=40, file='file__'//trim(identifier)//'.dat', status='unknown')
    write(40,*) 'Fermi_Energy', E_Fermi_tilde
    ! write(40,*) 'Smoothed_Total_Energy_Levels', E_tilde_shell
    close(40)


    contains 
    double precision function smoothed_level_density(e)
        use parameter
        implicit none
        double precision, intent(in) :: e

        integer :: j
        double precision ::  coeff


        smoothed_level_density = 0.0d0
        coeff = 1.0d0 / (gamma_shell*sqrt(pi))

        do j = 1, int(n_levels/2)
            smoothed_level_density = smoothed_level_density + &
            exp(- (e- eps(j))**2 / (gamma_shell**2) )
        end do
        smoothed_level_density = smoothed_level_density * coeff

    end function smoothed_level_density
end subroutine compute_smoothed_fermi_energy




subroutine effective_BCS_pairing
    !$ use omp_lib
    use parameter
    implicit none
    integer :: j
    double precision :: inv_bigG, gap_G,delta_bar ! delta_bar is pairing gap
    double precision, allocatable :: eps_2(:)
    print *, "Calculating effective BCS pairing energy for ", trim(identifier)
    hbar_omega_0 = 41.0d0 / (A_num)**(1.0d0/3.0d0)
    
    if (identifier == "proton") then
        gap_G = gap_Gp
        delta_bar = gap_p
    else if (identifier == "neutron") then
        gap_G = gap_Gn
        delta_bar = gap_n
    end if
    if (allocated(eps_2)) then
            deallocate(eps_2)
    end if
    ! allocate(eps_2(int(n_levels/2)))
    ! do j = 1, int(n_levels/2)
    !     eps_2(j) = eps(2*j-1)
    ! end do
    N_1 = 0
    N_2 = 0
    ! print *, "Estimated Fermi energy for ", trim(identifier), ": ", E_Fermi_tilde, " MeV"
    do j = 1, int(n_levels/2)
        if ( eps(j) < E_Fermi_tilde - hbar_omega_0 ) then
            N_1 = N_1 + 1
        else if ( E_Fermi_tilde - hbar_omega_0  <= eps(j)  .and. eps(j) <= E_Fermi_tilde + hbar_omega_0 ) then
            N_2 = N_2 + 1
        end if
    end do
    N_2 = N_2 + N_1
    N_1 = N_1 + 1

    ! rho_tilde = smoothed_level_density(E_Fermi_tilde)
    y_1 = (-0.5d0*N_tot + N_1 -1) / rho_tilde
    y_2 = (0.5d0*N_tot + N_2) / rho_tilde
    inv_bigG = A_num / 17.0d0  ! initial guess
    inv_bigG = 0.5d0 * rho_tilde * log( (y_2 + sqrt(y_2**2 + gap_G**2)) - &
     (y_1 + sqrt(y_1**2 + gap_G**2)) )
    bigG = 1.0d0 / inv_bigG
    ! E_p_c_tilde = -0.25d0*(N_2-N_1 +1)**2 / (rho_tilde) *((1.0d0+ (2*rho_tilde*gap_G/(N_2 -N_1 +1))**2 )**0.5 -1)+ &
    ! 0.5d0*gap_G*rho_tilde *bigG*atan( (N_2 -N_1 +1) / (2*rho_tilde*gap_G) )
    print *, " Effective paring gap  ", trim(identifier), ": ", gap_G, " MeV"
    print *, " Pairing strength  ", trim(identifier), ": ", bigG, " MeV"
    bigA = (rho_tilde * bigG *0.25d0)**2 * ((2+inv_bigG/rho_tilde) **2 -&
    (0.5d0*log(y_2**2 + gap_G**2) - 0.5d0*log (y_1**2 + gap_G**2) )**2 )
    bigB = gap_G**2 * rho_tilde**2 /16.0d0 * ( atan(y_2/gap_G) - atan(y_1/gap_G) )**2
    bigC = gap_G * rho_tilde /32.0d0 * (gap_G *(y_2/(y_2**2 + gap_G**2) - y_1/(y_1**2 + gap_G**2)) + &
    atan(y_2/gap_G) - atan(y_1/gap_G) ) 
    lambda_tilde_2 = bigG * 0.25d0* (bigA -bigC) / (bigB - bigC)
    E_p_c_tilde = 0.5d0*rho_tilde*((y_2-bigG)*(y_2 - sqrt(y_2**2 + gap_G**2)) + &
    (y_1 - bigG)*(sqrt(y_1**2 + gap_G**2) + y_1) ) + 0.25d0*(bigG - 4*lambda_tilde_2)* &
    ( atan(y_2/gap_G) - atan(y_1/gap_G) ) 
    print *, "Effective pairing energy for ", trim(identifier), ": ", E_p_c_tilde, " MeV"

    ! if (mod(N_tot,2) == 1) then
    !     E_p_c_tilde = E_p_c_tilde + delta_bar 
    ! end if
    open(unit=40, file='file__'//trim(identifier)//'.dat', position='append')
    write(40,*) 'Effective_Pairing_Gap', gap_G
    write(40,*) 'Effective_Pairing_Energy', E_p_c_tilde
    close(40)


    contains 
    double precision function smoothed_level_density(e)
        use parameter
        implicit none
        double precision, intent(in) :: e
        integer :: j
        double precision ::  coeff

        smoothed_level_density = 0.0d0
        coeff = 1.0d0 / (gamma_shell*sqrt(pi))

        do j = 1, int(n_levels/2)
            smoothed_level_density = smoothed_level_density + &
            exp(- (e-eps(j))**2 / (gamma_shell**2) )
        end do
        smoothed_level_density = smoothed_level_density * coeff

    end function smoothed_level_density
    
end subroutine effective_BCS_pairing

subroutine Shell_and_BCS
    !$ use omp_lib
    use parameter
    implicit none
    double precision :: E_fermi_lower, E_fermi_upper, E_fermi_mid, E_fermi_mid_prev  ! variables for Fermi energy search
    double precision :: delta_lower, delta_upper, delta_mid, delta_mid_prev ! variables for pairing gap search
    double precision :: N_calc ! calculated particle number
    double precision :: gap_eq_right ! pairing gap from eqref(95) of FRLDM2012
    ! double precision,allocatable :: inv_sqrt_energy(:) ! eqref(95) of FRLDM2012

    integer :: j,N
    integer :: max_j
    ! allocate(inv_sqrt_energy(N_tot))
    E_p_c = 0.0d0
    E_fermi_lower = -100.0d0
    E_fermi_upper = 50.0d0
    delta_lower = -20.0d0
    delta_upper = 20.0d0
    E_fermi_mid_prev = 0.0d0
    delta_mid_prev = 0.0d0
    E_fermi_mid = 0.5d0 * (E_fermi_lower + E_fermi_upper)
    delta_mid = 0.5d0 * (delta_lower + delta_upper)


    do while (abs(E_fermi_mid - E_fermi_mid_prev) > 1.0d-6 .or. abs(delta_mid - delta_mid_prev) > 1.0d-6)
        E_fermi_lower = -100.0d0
        E_fermi_upper = 50.0d0
        delta_lower = 0.0d0
        delta_upper = 20.0d0
        E_fermi_mid_prev = E_fermi_mid
        delta_mid_prev = delta_mid
        E_fermi_mid = 0.5d0 * (E_fermi_lower + E_fermi_upper)
        ! calculate N_1 : number of levels below Fermi energy not contributing to pairing 
        ! calculate N_2 : number of levels above Fermi energy
   
        ! compute E_fermi_mid and delta_mid dependent quantities
        do while (E_fermi_upper - E_fermi_lower > 1.0d-6)
            N_calc = 0.0d0
            !$omp parallel do default(none) schedule(guided) &
            !$omp shared(N_tot, eps, E_fermi_mid, delta_mid) &
            !$omp private(j) &
            !$omp reduction(+:N_calc)
            do j = N_1, N_2
                N_calc = N_calc + (1.0d0 - (eps(j) - E_fermi_mid) / &
                sqrt( (eps(j) - E_fermi_mid)**2 + delta_mid**2 ))
            end do
            !$omp end parallel do
            N_calc = N_calc + 2*(N_1-1)
            if (abs(N_calc - dble(N_tot)) < 1.0d-4) exit

            if (N_calc < dble(N_tot)) then
                E_fermi_lower = E_fermi_mid
            else
                E_fermi_upper = E_fermi_mid
            end if
            E_fermi_mid = 0.5d0 * (E_fermi_lower + E_fermi_upper)
        end do
        ! update delta_mid
        delta_mid = 0.5d0 * (delta_lower + delta_upper)
        do while (delta_upper - delta_lower > 1.0d-6)
            gap_eq_right = 0.0d0
            !$omp parallel do default(none) schedule(guided) &
            !$omp shared(N_tot, eps, E_fermi_mid, delta_mid) &
            !$omp private(j) &
            !$omp reduction(+:gap_eq_right)
            do j = N_1, N_2
                gap_eq_right = gap_eq_right +  1.0d0 / &
                sqrt( (eps(j) - E_fermi_mid)**2 + delta_mid**2 )
            end do
            !$omp end parallel do
            gap_eq_right = gap_eq_right * 0.5d0 * bigG
            if(abs(gap_eq_right - 1.0d0) < 1.0d-4) exit
            if (gap_eq_right > 1) then
                delta_lower = delta_mid
            else
                delta_upper = delta_mid
            end if
            delta_mid = 0.5d0 * (delta_lower + delta_upper)
        end do
    end do
    E_Fermi = E_fermi_mid
    delta = delta_mid
    print *, "BCS Fermi energy for ", trim(identifier), ": ", E_Fermi, " MeV"
    print *, "BCS pairing gap for ", trim(identifier), ": ", delta, " MeV"
    E_shell = 0.0d0
    !$omp parallel do default(none) schedule(guided) &
    !$omp shared(N_tot, eps, E_Fermi, delta) &
    !$omp private(j) &
    !$omp reduction(+:E_shell)
    do j = 1, int(n_tot/2)
        if ( eps(j) > E_Fermi ) cycle
        E_shell = E_shell + 2*eps(j)
    end do
    !$omp end parallel do
    print *, "Total energy levels for ", trim(identifier), ": ", E_shell, " MeV"

    E_p_c = 0.0d0
    !$omp parallel do default(none) schedule(guided) &
    !$omp shared(N_tot, eps, E_Fermi, delta,bigG) &
    !$omp private(j) &
    !$omp reduction(+:E_p_c)
    do j = N_1, N_2
        E_p_c = E_p_c + eps(j)*(1.0d0 - (eps(j) - E_Fermi) / &
        sqrt( (eps(j) - E_Fermi)**2 + delta**2 )) -&
         bigG /4.0d0 * &
        ( 1 - ((eps(j) - E_Fermi) / &
        sqrt( (eps(j) - E_Fermi)**2 + delta**2 )) )**2
    end do
    !$omp end parallel do
    
    N = int((N_1 + N_2) / 2)

    do j = N_1, N
        E_p_c = E_p_c - 2*eps(j) + bigG
    end do

    E_p_c = E_p_c - delta * delta / bigG
    print *, "Total pairing energy for ", trim(identifier), ": ", E_p_c, " MeV"
    print *, "N_1=", N_1, " N_2=", N_2
    print *, "N_calc=", N_calc, " N_tot=", N_tot
    print *, "gamma_shell=", gamma_shell


    if (identifier == "proton") then
        delta_E_shell_proton = E_shell - E_tilde_shell
        delta_E_p_c_proton = E_p_c - E_p_c_tilde
        print *, "Proton delta_E_shell=", delta_E_shell_proton
        print *, "Proton delta_E_p_c=", delta_E_p_c_proton
    else if (identifier == "neutron") then
        delta_E_shell_neutron = E_shell - E_tilde_shell
        delta_E_p_c_neutron = E_p_c - E_p_c_tilde
        print *, "Neutron delta_E_shell=", delta_E_shell_neutron
        print *, "Neutron delta_E_p_c=", delta_E_p_c_neutron
    end if
    open(unit=40, file='file__'//trim(identifier)//'.dat', position='append')
    write(40,*) 'BCS_Fermi_Energy', E_Fermi
    write(40,*) 'BCS_Pairing_Gap', delta
    write(40,*) 'BCS_Total_Energy_Levels', E_shell
    write(40,*) 'BCS_Pairing_Energy', E_p_c
    write(40,*) 'Delta_E_shell', E_shell - E_tilde_shell
    write(40,*) 'Delta_E_pairing', E_p_c - E_p_c_tilde
    write(40,*) 'N_1', N_1
    write(40,*) 'N_2', N_2
    write(40,*) 'N_calc', N_calc
    write(40,*) 'N_tot', N_tot
    write(40,*) 'gamma_shell', gamma_shell
    close(40)


end subroutine Shell_and_BCS

subroutine sum_FRLDM_MICRO
    use parameter
    implicit none
    
    print *, "Calculating total FRLDM + Microscopic energy..."
    FRLDM_MICRO_energy = FRLDM_energy + delta_E_shell_proton + delta_E_shell_neutron + &
    delta_E_p_c_proton + delta_E_p_c_neutron
    print *, "Total FRLDM + Microscopic energy=", FRLDM_MICRO_energy
    E_p_c_shell = delta_E_p_c_proton + delta_E_p_c_neutron + delta_E_shell_proton + delta_E_shell_neutron

    

end subroutine sum_FRLDM_MICRO

