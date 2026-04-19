module common_function_mod
    use iso_fortran_env, only: real64
    use constant_mod
    use nucleus_mod, only: nucleus_property
    use grid_mod, only: grid_type
    use frldm_constant_mod, only: a_den, a_Yukawa
    use micro_constant_mod, only: a_pot
    implicit none
    private

    integer, parameter :: dp = real64

end module common_function_mod