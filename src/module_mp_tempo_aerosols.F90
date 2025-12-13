module module_mp_tempo_aerosols
  !! initialize aerosols for tempo microphysics
  use module_mp_tempo_params, only : wp, sp, dp, &
    naccn0, naccn1, nain0, nain1, nwfa_default, aero_max

  implicit none
  private

  public :: init_water_friendly_aerosols, init_ice_friendly_aerosols, init_aerosol_emissions

  contains

  subroutine init_water_friendly_aerosols(hgt, nwfa)
    !! exponential profile of aerosols if nothing else is available
    real(wp), dimension(:), intent(in) :: hgt
    real(wp), dimension(:), intent(out) :: nwfa
    real(wp) :: h_01, niccn3
    integer :: k, nz
    
    nz = size(hgt)
    if(hgt(1) <= 1000.0_wp) then
      h_01 = 0.8_wp
    elseif(hgt(1) >= 2500.0_wp) then
      h_01 = 0.01_wp
    else
      h_01 = 0.8_wp*cos(hgt(1)*0.001_wp - 1.0_wp)
    endif
    niccn3 = -1.0_wp*log(naccn1/naccn0)/h_01
    nwfa(1) = naccn1+naccn0*exp(-((hgt(2)-hgt(1))/1000._wp)*niccn3)
    do k = 2, nz
      nwfa(k) = naccn1+naccn0*exp(-((hgt(k)-hgt(1))/1000._wp)*niccn3)
    enddo
  end subroutine init_water_friendly_aerosols


  subroutine init_ice_friendly_aerosols(hgt, nifa)
    !! exponential profile of aerosols if nothing else is available
    real(wp), dimension(:), intent(in) :: hgt
    real(wp), dimension(:), intent(out) :: nifa
    real(wp) :: h_01, niin3
    integer :: k, nz
    
    nz = size(hgt)
    if(hgt(1) <= 1000.0_wp) then
      h_01 = 0.8_wp
    elseif(hgt(1) >= 2500.0_wp) then
      h_01 = 0.01_wp
    else
      h_01 = 0.8_wp*cos(hgt(1)*0.001_wp - 1.0_wp)
    endif
    niin3 = -1.0_wp*log(nain1/nain0)/h_01
    nifa(1) = nain1+nain0*exp(-((hgt(2)-hgt(1))/1000._wp)*niin3)
    do k = 2, nz
      nifa(k) = nain1+nain0*exp(-((hgt(k)-hgt(1))/1000._wp)*niin3)
    enddo
  end subroutine init_ice_friendly_aerosols


  subroutine init_aerosol_emissions(deltaz, area, nwfa, nifa, nwfa2d, nifa2d)
    !! aerosol emissions calculated from lowest-model aerosol values
    ! scale the lowest level aerosol data into an emissions rate
    ! where: Nwfa=50 per cc, emit 0.875E4 aerosols per second per grid box unit
    ! that was tested as ~(20kmx20kmx50m = 2.e10 m**3) 
    real(wp), intent(in) :: deltaz, area, nwfa, nifa
    real(wp), intent(out) :: nwfa2d, nifa2d
    real(wp) :: deltaz_
    real(wp), parameter :: min_deltaz = 9._wp

    deltaz_ = deltaz
    if (deltaz < min_deltaz) deltaz_ = min_deltaz
    nwfa2d = max(nwfa_default, min(aero_max, nwfa)) * &
      0.000196_wp * (5._wp / deltaz_) * (area / 9.e6_wp)
    nifa2d = 0._wp
  end subroutine init_aerosol_emissions

end module module_mp_tempo_aerosols