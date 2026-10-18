module zero_point_mod
    use iso_fortran_env, only: real64 
    use constant_mod, only: pi, hbar_c
    use nucleus_mod, only: nucleus_property
    use micro_constant_mod
  implicit none
  private
  integer, parameter :: dp = real64

  type, abstract, public :: zero_point
    contains
      procedure(zero_point_E), deferred :: cal_energy
  end type zero_point

  abstract interface
    subroutine zero_point_E(this, nucleus, eps2, h, E_gs_lambda, E_0_lambda)
        import :: zero_point, dp, nucleus_property, microscopic_variables
        class(zero_point), intent(in) :: this
        type(nucleus_property), intent(in) :: nucleus
        real(dp), intent(in) :: h
        real(dp), intent(in) :: eps2
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

  

  contains
    subroutine calculate_E_0_eps2(this, nucleus, eps2, h, E_gs_lambda, E_0_lambda)
        class(zero_eps2), intent(in) :: this
        type(nucleus_property), intent(in) :: nucleus
        real(dp), intent(in) :: eps2
        real(dp), intent(in) :: h
        real(dp), intent(in) :: E_gs_lambda(7)
        real(dp), intent(out) :: E_0_lambda
        real(dp) :: R0, M0, eps2_2, eps2_3, B_eps2, C_eps2, M0R0_2
        R0 = nucleus%R0
        M0 = nucleus%mass
        M0R0_2 = M0 * R0**2
        
        eps2_2 = eps2**2
        eps2_3 = eps2**3
        B_eps2 =  2.0_dp/15_dp * (1.0_dp + 2.0_dp/9_dp * eps2_2) &
                * (1.0_dp - 1.0_dp/3_dp * eps2_2 - 2.0_dp/27_dp * eps2_3)**(-4.0_dp/3.0_dp) &
                * M0R0_2 / frK**2
        C_eps2 =  seven_points_derivative(E_gs_lambda, h) 

        E_0_lambda = 0.5_dp * sqrt(C_eps2 / B_eps2)*hbar_c

    end subroutine calculate_E_0_eps2

    subroutine calculate_E_0_gamma(this, nucleus, eps2, h, E_gs_lambda, E_0_lambda)
        class(zero_gamma), intent(in) :: this
        type(nucleus_property), intent(in) :: nucleus
        real(dp), intent(in) :: eps2
        real(dp), intent(in) :: h
        real(dp), intent(in) :: E_gs_lambda(7)
        real(dp), intent(out) :: E_0_lambda
        real(dp) :: R0, M0, B_gamma, C_gamma, M0R0_2
        R0 = nucleus%R0
        M0 = nucleus%mass
        M0R0_2 = M0 * R0**2
        B_gamma =  2.0_dp/15_dp * ((1.0_dp - 2.0_dp/3_dp * eps2)/(1.0_dp + 2.0_dp/3_dp * eps2))**(2.0_dp/3.0_dp) &
                * log((1.0_dp + 1.0_dp/3.0_dp * eps2) / (1.0_dp - 2.0_dp/3.0_dp * eps2) )**2 &
                * M0R0_2 / frK**2
        C_gamma =  seven_points_derivative(E_gs_lambda, h)

        E_0_lambda = 0.5_dp * sqrt(C_gamma / B_gamma)*hbar_c

    end subroutine calculate_E_0_gamma

    function seven_points_derivative(func,h)
        real(dp), intent(in) :: h
        real(dp), intent(in) :: func(7)
        real(dp) :: seven_points_derivative
        seven_points_derivative = 2.0_dp * func(1) - 27.0_dp * func(2) + 270.0_dp * func(3) - 490.0_dp * func(4) &
                                 + 270.0_dp * func(5) - 27.0_dp * func(6) + 2.0_dp * func(7)
        seven_points_derivative = seven_points_derivative / (180.0_dp * h**2)
    end function seven_points_derivative

end module zero_point_mod
