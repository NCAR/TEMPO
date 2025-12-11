module test_tempo_init_suite
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use module_mp_tempo_params, only : wp, sp, dp, &
    initialize_graupel_vars, dim_nrhg, nrhg, nrhg1, rho_g, am_g, pi, &
    initialize_parameters, ccg, &
    initialize_bins_for_tables, dc, nbc, t_efrw, t_efsw, &
    initialize_array_efrw, initialize_array_efsw
  use module_mp_tempo_init, only : compute_efrw, compute_efsw
  implicit none
  private

  public :: collect_tempo_init_suite

  contains

  ! collect all exported unit tests
  subroutine collect_tempo_init_suite(testsuite)
    ! collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("dim_nrhg=9 with hail aware = true", &
        InitializeGraupelVars_dimNrhg9_when_HailAwareTrue), &
      new_unittest("dim_nrhg=1 with hail aware = false", &
        InitializeGraupelVars_dimNrhg1_when_HailAwareFalse), &
      new_unittest("ccg(2,4) = 7!, or 5040", &
        InitializeParameters_ccgValue5040), &
      new_unittest("dc(nbc) = 100 microns (last dc bin)", &
        InitializeBinsForTables_lastDcBin100microns), &
      new_unittest("all values of t_efrw = 0 is false", &
        TableEfrw_AllZeroFalse), &
      new_unittest("all values of t_efsw = 0 is false", &
        TableEfsw_AllZeroFalse) &
      ]
  end subroutine collect_tempo_init_suite

    
  subroutine InitializeGraupelVars_dimNrhg9_when_HailAwareTrue(error)
    type(error_type), allocatable, intent(out) :: error

    call initialize_graupel_vars(hail_flag=.true.)
    call check(error, dim_nrhg, nrhg, "dim_nrhg with hail_aware = true")
    if (allocated(error)) return
  end subroutine InitializeGraupelVars_dimNrhg9_when_HailAwareTrue


  subroutine InitializeGraupelVars_dimNrhg1_when_HailAwareFalse(error)
    type(error_type), allocatable, intent(out) :: error

    call initialize_graupel_vars(hail_flag=.false.)
    call check(error, dim_nrhg, nrhg1, "dim_nrhg with hail_aware = false")
    if (allocated(error)) return
  end subroutine InitializeGraupelVars_dimNrhg1_when_HailAwareFalse


  subroutine InitializeParameters_ccgValue5040(error)
    type(error_type), allocatable, intent(out) :: error

    call initialize_parameters()
    call check(error, ccg(2,4), 5040._wp, "ccg(2,4) = 7!, or 5040")
    if (allocated(error)) return
  end subroutine InitializeParameters_ccgValue5040


  subroutine InitializeBinsForTables_lastDcBin100microns(error)
    type(error_type), allocatable, intent(out) :: error

    call initialize_bins_for_tables()
    call check(error, dc(nbc), 100.e-6_wp, "dc last bin = 100 microns")
    if (allocated(error)) return
  end subroutine InitializeBinsForTables_lastDcBin100microns


  subroutine TableEfrw_allZeroFalse(error)
    type(error_type), allocatable, intent(out) :: error
    call initialize_array_efrw()
    call compute_efrw()
    call check(error, all(t_efrw==0.0_dp), .false., "all values of t_efrw = 0 is false")
    if (allocated(error)) return
  end subroutine TableEfrw_AllZeroFalse


  subroutine TableEfsw_allZeroFalse(error)
    type(error_type), allocatable, intent(out) :: error
    call initialize_array_efsw()
    call compute_efsw()
    call check(error, all(t_efsw==0.0_dp), .false., "all values of t_efsw = 0 is false")
    if (allocated(error)) return
  end subroutine TableEfsw_AllZeroFalse

end module test_tempo_init_suite