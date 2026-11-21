! ------------------------------------------------------------------
! TEMPO_STRIDE1 compile-time specialization (per-column CPU fast build)
! When -DTEMPO_STRIDE1 is set, the horizontal block iterators its/ite/jts/jte
! (i.e. bi/bj) are hardcoded to a compile-time extent of 1 so the compiler can
! drop the horizontal dimension, remove 3D address arithmetic and specialize the
! vertical (k) loops. The driver forces block_stride=1 in this build, so every
! invocation runs the per-column form regardless of the requested -s stride.
! Without the macro the tokens fall back to the original runtime bounds, so the
! default build is byte-for-byte the original hybrid code.
! ------------------------------------------------------------------
#ifdef TEMPO_STRIDE1
#define TEMPO_ITS 1
#define TEMPO_ITE 1
#define TEMPO_JTS 1
#define TEMPO_JTE 1
#else
#define TEMPO_ITS its
#define TEMPO_ITE ite
#define TEMPO_JTS jts
#define TEMPO_JTE jte
#endif
module module_mp_tempo_aerosols
  !! contains produces used when aerosol-aware = true
  use module_mp_tempo_params, only : wp, sp, dp, eps, &
    naccn0, naccn1, nain0, nain1, nwfa_default, aero_max

  implicit none
  private

  public :: init_water_friendly_aerosols, init_ice_friendly_aerosols, &
    aerosol_collection_efficiency

  contains

  subroutine init_water_friendly_aerosols(kts, kte, its, ite, jts, jte, tempo_first_main, dz3d, nwfa, nwfa3d)
    !! column_j1b / column_i1b: initialize water-friendly aerosols over the horizontal tile
    !! sets water-friendly aerosols to an exponential profile when aerosol-aware = true
    !! and no initial condition is provided by the host model
    !! @note requires kte > kts (at least two vertical mass levels)
    integer, intent(in) :: kts, kte, its, ite, jts, jte
    logical, intent(in) :: tempo_first_main
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: dz3d
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: nwfa
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: nwfa3d
    integer :: i, j, k
    real(wp) :: hgt(kts:kte)
    real(wp) :: h_01, niccn3

    if (present(nwfa3d)) then
      !$acc parallel default(present) copyin(dz3d) copy(nwfa) copyin(nwfa3d)
      !$acc loop gang collapse(2) private(hgt, h_01, niccn3)
      column_j1b_p: do j = TEMPO_JTS, TEMPO_JTE
        column_i1b_p: do i = TEMPO_ITS, TEMPO_ITE
          if (tempo_first_main .and. j == jts .and. i == its .and. sum(nwfa3d(:,i,j)) < eps) then
            hgt = 0._wp
            !$acc loop seq
            do k = kts+1, kte
              hgt(k) = hgt(k-1) + dz3d(k,i,j)
            enddo
            if(hgt(kts) <= 1000.0_wp) then
              h_01 = 0.8_wp
            elseif(hgt(kts) >= 2500.0_wp) then
              h_01 = 0.01_wp
            else
              h_01 = 0.8_wp*cos(hgt(kts)*0.001_wp - 1.0_wp)
            endif
            niccn3 = -1.0_wp*log(naccn1/naccn0)/h_01
            nwfa(kts,i,j) = naccn1+naccn0*exp(-((hgt(kts+1)-hgt(kts))/1000._wp)*niccn3)
            !$acc loop vector
            do k = kts+1, kte
              nwfa(k,i,j) = naccn1+naccn0*exp(-((hgt(k)-hgt(kts))/1000._wp)*niccn3)
            enddo
          endif
        enddo column_i1b_p
      enddo column_j1b_p
      !$acc end parallel
    else
      !$acc parallel default(present) copyin(dz3d) copy(nwfa)
      !$acc loop gang collapse(2) private(hgt, h_01, niccn3)
      column_j1b: do j = TEMPO_JTS, TEMPO_JTE
        column_i1b: do i = TEMPO_ITS, TEMPO_ITE
          hgt = 0._wp
          !$acc loop seq
          do k = kts+1, kte
            hgt(k) = hgt(k-1) + dz3d(k,i,j)
          enddo
          if(hgt(kts) <= 1000.0_wp) then
            h_01 = 0.8_wp
          elseif(hgt(kts) >= 2500.0_wp) then
            h_01 = 0.01_wp
          else
            h_01 = 0.8_wp*cos(hgt(kts)*0.001_wp - 1.0_wp)
          endif
          niccn3 = -1.0_wp*log(naccn1/naccn0)/h_01
          nwfa(kts,i,j) = naccn1+naccn0*exp(-((hgt(kts+1)-hgt(kts))/1000._wp)*niccn3)
          !$acc loop vector
          do k = kts+1, kte
            nwfa(k,i,j) = naccn1+naccn0*exp(-((hgt(k)-hgt(kts))/1000._wp)*niccn3)
          enddo
        enddo column_i1b
      enddo column_j1b
      !$acc end parallel
    endif
  end subroutine init_water_friendly_aerosols


  subroutine init_ice_friendly_aerosols(kts, kte, its, ite, jts, jte, tempo_first_main, dz3d, nifa, nifa3d)
    !! column_j1b2 / column_i1b2: initialize ice-friendly aerosols over the horizontal tile
    !! sets ice-friendly aerosols to an exponential profile when aerosol-aware = true
    !! and no initial condition is provided by the host model
    !! @note requires kte > kts (at least two vertical mass levels)
    integer, intent(in) :: kts, kte, its, ite, jts, jte
    logical, intent(in) :: tempo_first_main
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in) :: dz3d
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(inout) :: nifa
    real(wp), dimension(kts:kte, TEMPO_ITS:TEMPO_ITE, TEMPO_JTS:TEMPO_JTE), intent(in), optional :: nifa3d
    integer :: i, j, k
    real(wp) :: hgt(kts:kte)
    real(wp) :: h_01, niin3

    if (present(nifa3d)) then
      !$acc parallel default(present) copyin(dz3d) copy(nifa) copyin(nifa3d)
      !$acc loop gang collapse(2) private(hgt, h_01, niin3)
      column_j1b2_p: do j = TEMPO_JTS, TEMPO_JTE
        column_i1b2_p: do i = TEMPO_ITS, TEMPO_ITE
          if (tempo_first_main .and. j == jts .and. i == its .and. sum(nifa3d(:,i,j)) < eps) then
            hgt = 0._wp
            !$acc loop seq
            do k = kts+1, kte
              hgt(k) = hgt(k-1) + dz3d(k,i,j)
            enddo
            if(hgt(kts) <= 1000.0_wp) then
              h_01 = 0.8_wp
            elseif(hgt(kts) >= 2500.0_wp) then
              h_01 = 0.01_wp
            else
              h_01 = 0.8_wp*cos(hgt(kts)*0.001_wp - 1.0_wp)
            endif
            niin3 = -1.0_wp*log(nain1/nain0)/h_01
            nifa(kts,i,j) = nain1+nain0*exp(-((hgt(kts+1)-hgt(kts))/1000._wp)*niin3)
            !$acc loop vector
            do k = kts+1, kte
              nifa(k,i,j) = nain1+nain0*exp(-((hgt(k)-hgt(kts))/1000._wp)*niin3)
            enddo
          endif
        enddo column_i1b2_p
      enddo column_j1b2_p
      !$acc end parallel
    else
      !$acc parallel default(present) copyin(dz3d) copy(nifa)
      !$acc loop gang collapse(2) private(hgt, h_01, niin3)
      column_j1b2: do j = TEMPO_JTS, TEMPO_JTE
        column_i1b2: do i = TEMPO_ITS, TEMPO_ITE
          hgt = 0._wp
          !$acc loop seq
          do k = kts+1, kte
            hgt(k) = hgt(k-1) + dz3d(k,i,j)
          enddo
          if(hgt(kts) <= 1000.0_wp) then
            h_01 = 0.8_wp
          elseif(hgt(kts) >= 2500.0_wp) then
            h_01 = 0.01_wp
          else
            h_01 = 0.8_wp*cos(hgt(kts)*0.001_wp - 1.0_wp)
          endif
          niin3 = -1.0_wp*log(nain1/nain0)/h_01
          nifa(kts,i,j) = nain1+nain0*exp(-((hgt(kts+1)-hgt(kts))/1000._wp)*niin3)
          !$acc loop vector
          do k = kts+1, kte
            nifa(k,i,j) = nain1+nain0*exp(-((hgt(k)-hgt(kts))/1000._wp)*niin3)
          enddo
        enddo column_i1b2
      enddo column_j1b2
      !$acc end parallel
    endif
  end subroutine init_ice_friendly_aerosols


  function aerosol_collection_efficiency(d, da, visc, rhoa, temp, species) result(eff_a)
  !$acc routine seq
    !! computes aerosol collection efficiency for precipitation scavenging
    !! from [Wang et al. (2010)](https://doi.org/10.5194/acp-10-5685-2010)
    use module_mp_tempo_params, only : rho_w, rho_s, av_s, bv_s, &
      idx_bg1, pi, av_g, bv_g, rho_g

    real(dp), intent(in) :: d
    real(wp), intent(in) :: da, visc, rhoa, temp
    character(len=1), intent(in) :: species
    real(wp) :: aval, cc, diff, re, sc, st, st2, vt, eff, rho_p
    real(wp), parameter :: boltzman = 1.3806503e-23_wp
    real(wp), parameter :: meanpath = 0.0256e-6_wp
    real(wp) :: eff_a

    vt = 1._wp
    rho_p = rho_w
    ! rain
    if (species == 'r') then
      vt = -0.1021_wp + 4.932e3_wp*d - 0.9551e6_wp*d*d + &
        0.07934e9_wp*d*d*d - 0.002362e12_wp*d*d*d*d
      rho_p = rho_w
    ! snow
    elseif (species == 's') then
      vt = av_s*d**bv_s
      rho_p = rho_s
    ! graupel
    elseif (species .eq. 'g') then
      vt = av_g(idx_bg1)*d**bv_g(idx_bg1)
      rho_p = rho_g(idx_bg1)
    endif
    cc = 1._wp + 2._wp*meanpath/da *(1.257_wp+0.4_wp*exp(-0.55_wp*da/meanpath))
    diff = boltzman*temp*cc/(3._wp*pi*visc*da)
    re = 0.5_wp*rhoa*d*vt/visc
    sc = visc/(rhoa*diff)
    st = (rho_p-rhoa)*da*da*vt*cc/(9._wp*visc*d)
    aval = log(1._wp + re)
    st2 = (1.2_wp + 1._wp/12._wp*aval)/(1._wp+aval)
    eff = 4._wp/(re*sc) * (1._wp + 0.4_wp*sqrt(re)*sc**0.3333_wp + &
      0.16_wp*sqrt(re)*sqrt(sc)) + 4._wp*da/d * &
      (0.02_wp + da/d*(1._wp+2._wp*sqrt(re)))
    if (st > st2) eff = eff + ((st-st2)/(st-st2+0.666667_wp))**1.5_wp
    eff_a = max(1.e-5_wp, min(eff, 1._wp))
  end function aerosol_collection_efficiency

end module module_mp_tempo_aerosols
