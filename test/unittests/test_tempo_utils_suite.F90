module test_tempo_utils_suite
  !! unit tests for module_mp_tempo_utils
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use module_mp_tempo_params, only : wp, sp, dp, t_efrw, t_efsw, &
    initialize_array_efrw, initialize_array_efsw, initialize_bins_for_tables
  use module_mp_tempo_utils, only : calc_gamma_p, compute_efrw, compute_efsw, get_nuc
  use module_mp_tempo_cfgs, only : ty_tempo_cfgs, ty_tempo_spp
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
        CalcGammaP_a8x8), &
      new_unittest("checking initialization of rain-cloud water collection efficiency data after allocation resulting in all(t_efrw == 0.) being false", &
        TableEfrw_AllZeroFalse), &
      new_unittest("checking initialization of snow-cloud water collection efficiency data after allocation resulting in all(t_efsw == 0.) being false", &
        TableEfsw_AllZeroFalse), &
      new_unittest("checking that get_nuc applies the stochastic shape offset and keeps nu_c within 2-15", &
        GetNuc_ShapeOffset), &
      new_unittest("checking that set_spp_options decodes the spp_mp bitmask (0-7)", &
        SetSppOptions_Bitmask), &
      new_unittest("checking the per-column perturbations from spp_perturbations", &
        SppPerturbations_Values) &
      ]

  end subroutine collect_tempo_utils_suite

  subroutine CalcGammaP_a8x8(error)
    !! test calc_gamma_p with a=8 and x=8
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: value_from_wolfram = 0.54703919051300551_wp

    call check(error, abs(calc_gamma_p(8._wp, 8._wp) - value_from_wolfram) < 1.e-6_wp)
    if (allocated(error)) return
  end subroutine CalcGammaP_a8x8


  subroutine TableEfrw_allZeroFalse(error)
    !! test for proper initialization of t_efrw
    type(error_type), allocatable, intent(out) :: error
    call initialize_bins_for_tables() 
    call initialize_array_efrw()
    write(*,*) shape(t_efrw)
    call compute_efrw()
    call check(error, all(t_efrw==0.0_dp), .false.)
    if (allocated(error)) return
  end subroutine TableEfrw_AllZeroFalse


  subroutine TableEfsw_allZeroFalse(error)
    !! test for proper initialization of t_efsw
    type(error_type), allocatable, intent(out) :: error
    call initialize_bins_for_tables() 
    call initialize_array_efsw()
    call compute_efsw()
    call check(error, all(t_efsw==0.0_dp), .false.)
    if (allocated(error)) return
  end subroutine TableEfsw_AllZeroFalse


  subroutine GetNuc_ShapeOffset(error)
    !! nc = 100 cm-3 gives nu_c = 12; offsets shift it within 2-15; low nc stays at 15
    type(error_type), allocatable, intent(out) :: error

    call check(error, get_nuc(100.e6_wp), 12)
    if (allocated(error)) return
    call check(error, get_nuc(100.e6_wp, 0), 12)
    if (allocated(error)) return
    call check(error, get_nuc(100.e6_wp, 2), 14)
    if (allocated(error)) return
    call check(error, get_nuc(100.e6_wp, 5), 15)
    if (allocated(error)) return
    call check(error, get_nuc(100.e6_wp, -20), 2)
    if (allocated(error)) return
    call check(error, get_nuc(50.e6_wp, -5), 15)
    if (allocated(error)) return
  end subroutine GetNuc_ShapeOffset


  subroutine SetSppOptions_Bitmask(error)
    !! bit 1 = graupel intercept, bit 2 = cloud shape, bit 4 = nucleation
    type(error_type), allocatable, intent(out) :: error
    type(ty_tempo_cfgs) :: cfgs
    integer :: spp_mp

    do spp_mp = 0, 7
      call cfgs%set_spp_options(spp_mp)
      call check(error, cfgs%spp_graupel_intercept_flag .eqv. (mod(spp_mp, 2) == 1))
      if (allocated(error)) return
      call check(error, cfgs%spp_cloud_shape_flag .eqv. (mod(spp_mp/2, 2) == 1))
      if (allocated(error)) return
      call check(error, cfgs%spp_nucleation_flag .eqv. (mod(spp_mp/4, 2) == 1))
      if (allocated(error)) return
    enddo
    call cfgs%set_spp_options(7, spp_amplitude=0.5_wp, spp_stddev_cutoff=3._wp)
    call check(error, abs(cfgs%spp_pattern_max - 1.5_wp) < 1.e-12_wp)
    if (allocated(error)) return
  end subroutine SetSppOptions_Bitmask


  subroutine SppPerturbations_Values(error)
    !! pattern 0.4 with a maximum of 1: intercept +0.4, shape +nint(0.8) = 1, nucleation 0.25*(0.4+1)
    type(error_type), allocatable, intent(out) :: error
    type(ty_tempo_cfgs) :: cfgs
    type(ty_tempo_spp) :: spp

    ! all options off: no perturbation whatever the pattern
    spp = cfgs%spp_perturbations(0.4_wp)
    call check(error, spp%graupel_log10_n0_offset, 0._wp)
    if (allocated(error)) return
    call check(error, spp%cloud_shape_offset, 0)
    if (allocated(error)) return
    call check(error, spp%nucleation_offset, 0._wp)
    if (allocated(error)) return

    call cfgs%set_spp_options(7, spp_amplitude=1._wp, spp_stddev_cutoff=1._wp)
    spp = cfgs%spp_perturbations(0.4_wp)
    call check(error, abs(spp%graupel_log10_n0_offset - 0.4_wp) < 1.e-12_wp)
    if (allocated(error)) return
    call check(error, spp%cloud_shape_offset, 1)
    if (allocated(error)) return
    call check(error, abs(spp%nucleation_offset - 0.35_wp) < 1.e-12_wp)
    if (allocated(error)) return

    ! a pattern beyond the cutoff cannot make the nucleation perturbation negative
    spp = cfgs%spp_perturbations(-2._wp)
    call check(error, spp%nucleation_offset, 0._wp)
    if (allocated(error)) return
    call check(error, spp%cloud_shape_offset, -4)
    if (allocated(error)) return
  end subroutine SppPerturbations_Values

end module test_tempo_utils_suite