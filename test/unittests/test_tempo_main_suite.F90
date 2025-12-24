module test_tempo_main_suite
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use module_mp_tempo_params, only : wp, sp, dp, r_c, t_nc, r_s, tc, t0, &  
    r_r, n0r_exp, crg, org2, org1, bm_r, org1, am_r, cre, mu_r, &
    idx_bg1, r_g, n0g_exp, cgg, ogg2, ogg1, bm_g, am_g, cge, r_i, nt_i
  use module_mp_tempo_main, only : get_cloud_table_index, get_snow_table_index, &
    get_temperature_table_index, get_rain_table_index, get_graupel_table_index, &
    get_ice_table_index
  implicit none
  private

  public :: collect_tempo_main_suite

  contains

  ! collect all exported unit tests
  subroutine collect_tempo_main_suite(testsuite)
    ! collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("get_cloud_table_index(rc=1.2345e-5) returns index where rc is between r_c(index-1) and r_c(index+1)", &
        getCloudTableIndex_massIndexCorrect), &
      new_unittest("get_cloud_table_index(nc=54.321e6) returns index where nc is between t_nc(index-1) and t_nc(index+1)", &
        getCloudTableIndex_numberIndexCorrect), &
      new_unittest("get_snow_table_index(rs=2.3456e-4) returns index where rs is between r_(index-1) and r_s(index+1)", &
        getSnowTableIndex_massIndexCorrect), &
      new_unittest("get_temperature_table_index(temp=-23.4) returns index where temp is between tc(index-1) and tc(index+1)", &
        getTemperatureTableIndex_IndexCorrect), &
      new_unittest("get_rain_table_index(rr=1.2345e-4) returns index where rr is between r_r(index-1) and r_r(index+1)", &
        getRainTableIndex_massIndexCorrect), &
      new_unittest("get_rain_table_index(ilam=3.4567e8) returns index where n0exp is between n0exp(index-1) and n0exp(index+1)", &
        getRainTableIndex_n0expIndexCorrectMvd), &
      new_unittest("get_rain_table_index(rr=1.2345e-4) returns index where rr is between r_r(index-1) and r_r(index+1)", &
        getGraupelTableIndex_massIndexCorrect), &
      new_unittest("get_rain_table_index(ilam=3.4567e8) returns index where n0exp is between n0exp(index-1) and n0exp(index+1)", &
        getGraupelTableIndex_n0expIndexCorrectMvd), &
      new_unittest("get_cloud_table_index(rc=1.2345e-5) returns index where rc is between r_c(index-1) and r_c(index+1)", &
        getIceTableIndex_massIndexCorrect), &
      new_unittest("get_cloud_table_index(nc=54.321e6) returns index where nc is between t_nc(index-1) and t_nc(index+1)", &
        getIceTableIndex_numberIndexCorrect) &
      ]

  end subroutine collect_tempo_main_suite

  subroutine getCloudTableIndex_massIndexCorrect(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_c, idx_n
    real(wp) :: mass, number
    mass = 1.2345e-4
    number = 54.321e6
    call get_cloud_table_index(mass, number, idx_c, idx_n)

    call check(error, (mass > r_c(idx_c-1) .and. mass < r_c(idx_c+1)), &
      "get_cloud_table_index(rc=1.2345e-5) returns index where rc is between r_c(index-1) and r_c(index+1)")
    if (allocated(error)) return
  end subroutine getCloudTableIndex_massIndexCorrect


  subroutine getCloudTableIndex_numberIndexCorrect(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_c, idx_n
    real(wp) :: mass, number
    mass = 1.2345e-4
    number = 54.321e6
    call get_cloud_table_index(mass, number, idx_c, idx_n)

    call check(error, (number > t_nc(idx_n-1) .and. number < t_nc(idx_n+1)), &
      "get_cloud_table_index(nc=54.321e6) returns index where nc is between t_nc(index-1) and t_nc(index+1)")
    if (allocated(error)) return
  end subroutine getCloudTableIndex_numberIndexCorrect


  subroutine getSnowTableIndex_massIndexCorrect(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_s
    real(wp) :: mass
    mass = 2.3456e-4
    call get_snow_table_index(mass, idx_s)

    call check(error, (mass > r_s(idx_s-1) .and. mass < r_s(idx_s+1)), &
      "get_snow_table_index(rc=2.3456e-4) returns index where rs is between r_s(index-1) and r_s(index+1)")
    if (allocated(error)) return
  end subroutine getSnowTableIndex_massIndexCorrect


  subroutine getTemperatureTableIndex_IndexCorrect(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_t
    real(wp) :: temp
    temp = t0 - 23.4_wp
    call get_temperature_table_index(temp, idx_t)

    call check(error, (temp-t0 < tc(idx_t-1) .and. temp-t0 > tc(idx_t+1)), &
      "get_temperture_table_index(temp=-23.4) returns index where temp is between tc(index-1) and tc(index+1)")
    if (allocated(error)) return
  end subroutine getTemperatureTableIndex_IndexCorrect


  subroutine getRainTableIndex_massIndexCorrect(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_r, idx_r1
    real(wp) :: mass, ilam
    mass = 1.2345e-4_wp
    ilam = 3.4567e8_wp
    call get_rain_table_index(mass, ilam, idx_r, idx_r1)

    call check(error, (mass > r_c(idx_r-1) .and. mass < r_r(idx_r+1)), &
      "get_rain_table_index(rr=1.2345e-4) returns index where rr is between r_r(index-1) and r_r(index+1)")
    if (allocated(error)) return
  end subroutine getRainTableIndex_massIndexCorrect


  subroutine getRainTableIndex_n0expIndexCorrectMvd(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_r, idx_r1
    real(wp) :: mass, ilam, lam, mvd, lam_exp_lower, lam_exp_upper, &
      lam_lower, lam_upper, mvd_lower, mvd_upper
    mass = 1.2345e-3_wp
    mvd = 1.234e-3_wp ! 1.234 mm
    lam = (3.0_wp + mu_r + 0.672_wp) / mvd
    ilam = 1._wp/lam
    call get_rain_table_index(mass, ilam, idx_r, idx_r1)
    lam_exp_upper = (n0r_exp(idx_r1-1)/(org1*mass/am_r))**(1._wp/cre(1))
    lam_exp_lower = (n0r_exp(idx_r1+1)/(org1*mass/am_r))**(1._wp/cre(1))
    lam_upper = lam_exp_upper / ((crg(3)*org2*org1)**bm_r)
    lam_lower = lam_exp_lower / ((crg(3)*org2*org1)**bm_r)
    mvd_upper = (3.0_wp + mu_r + 0.672_wp) / lam_upper
    mvd_lower = (3.0_wp + mu_r + 0.672_wp) / lam_lower

    call check(error, (mvd > mvd_lower .and. mvd < mvd_upper), &
      "get_rain_table_index(ilamr) returns index where mvd is between mvd(index-1) and mvd(index+1)")
    if (allocated(error)) return
  end subroutine getRainTableIndex_n0expIndexCorrectMvd


  subroutine getGraupelTableIndex_massIndexCorrect(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_g, idx_g1
    real(wp) :: mass, ilam
    integer :: idx
    mass = 1.2345e-4_wp
    ilam = 3.4567e8_wp
    idx = idx_bg1
    call get_graupel_table_index(mass, ilam, idx, idx_g, idx_g1)

    call check(error, (mass > r_g(idx_g-1) .and. mass < r_g(idx_g+1)), &
      "get_graupel_table_index(rr=1.2345e-4) returns index where rg is between r_g(index-1) and r_g(index+1)")
    if (allocated(error)) return
  end subroutine getGraupelTableIndex_massIndexCorrect


  subroutine getGraupelTableIndex_n0expIndexCorrectMvd(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_g, idx_g1
    integer :: idx
    real(wp) :: mass, ilam, lam, mvd, lam_exp_lower, lam_exp_upper, &
      lam_lower, lam_upper, mvd_lower, mvd_upper
    mass = 1.2345e-3_wp
    mvd = 5.6789e-3_wp ! 5.6789 mm
    idx = idx_bg1
    lam = (3.0_wp + mu_r + 0.672_wp) / mvd
    ilam = 1._wp/lam
    call get_graupel_table_index(mass, ilam, idx, idx_g, idx_g1)
    lam_exp_upper = (n0g_exp(idx_g1-1)/(ogg1*mass/am_g(idx)))**(1._wp/cge(1,1))
    lam_exp_lower = (n0g_exp(idx_g1+1)/(ogg1*mass/am_g(idx)))**(1._wp/cge(1,1))
    lam_upper = lam_exp_upper / ((cgg(3,1)*ogg2*ogg1)**bm_g)
    lam_lower = lam_exp_lower / ((cgg(3,1)*ogg2*ogg1)**bm_g)
    mvd_upper = (3.0_wp + mu_r + 0.672_wp) / lam_upper
    mvd_lower = (3.0_wp + mu_r + 0.672_wp) / lam_lower

    call check(error, (mvd > mvd_lower .and. mvd < mvd_upper), &
      "get_rain_table_index(ilamr) returns index where mvd is between mvd(index-1) and mvd(index+1)")
    if (allocated(error)) return
  end subroutine getGraupelTableIndex_n0expIndexCorrectMvd


  subroutine getIceTableIndex_massIndexCorrect(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_i, idx_i1
    real(wp) :: mass, number
    mass = 1.2345e-5
    number = 100000. ! 100 / L
    call get_ice_table_index(mass, number, idx_i, idx_i1)

    call check(error, (mass > r_i(idx_i-1) .and. mass < r_i(idx_i+1)), &
      "get_cloud_table_index(rc=1.2345e-5) returns index where rc is between r_c(index-1) and r_c(index+1)")
    if (allocated(error)) return
  end subroutine getIceTableIndex_massIndexCorrect


  subroutine getIceTableIndex_numberIndexCorrect(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: idx_i, idx_i1
    real(wp) :: mass, number
    mass = 1.2345e-5
    number = 100000. ! 100 / L
    call get_ice_table_index(mass, number, idx_i, idx_i1)

    call check(error, (number > nt_i(idx_i1-1) .and. number < nt_i(idx_i1+1)), &
      "get_cloud_table_index(nc=54.321e6) returns index where nc is between t_nc(index-1) and t_nc(index+1)")
    if (allocated(error)) return
  end subroutine getIceTableIndex_numberIndexCorrect

end module test_tempo_main_suite
