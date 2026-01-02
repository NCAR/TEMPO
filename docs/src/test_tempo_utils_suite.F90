module test_tempo_utils_suite
  !! unit tests for module_mp_tempo_utils
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use module_mp_tempo_params, only : wp, sp, dp
  use module_mp_tempo_utils, only : calc_gamma_p
  implicit none
  private

  public :: collect_tempo_utils_suite

  contains

  ! collect all exported unit tests
  subroutine collect_tempo_utils_suite(testsuite)
    ! collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("checking that procedure calc_gamma_p(8.,8.) = 0.547", &
        CalcGammaP_a8x8) &
      ]

  end subroutine collect_tempo_utils_suite

  subroutine CalcGammaP_a8x8(error)
    !! test calc_gamma_p with a=8 and x=8
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: value_from_wolfram = 0.54703919051300551 

    call check(error, abs(calc_gamma_p(8.,8.) - value_from_wolfram) < 1.e-6_wp)
    if (allocated(error)) return
  end subroutine CalcGammaP_a8x8

end module test_tempo_utils_suite